#!/usr/bin/env python3

import argparse
import os
import sys
import re
from ete3 import Tree
import csv
import pandas as pd
from collections import defaultdict

########################################################################################################################
# HELPER FUNCTIONS
########################################################################################################################
def reroot_by_tips(t, tip_names):
    """Root at the LCA of listed tips; raises if a tip is missing.
    Returns the list of matched tip nodes."""
    nodes = []
    for nm in tip_names:
        hits = t.search_nodes(name=nm)
        if not hits:
            raise RuntimeError(f"Tip not found for --root-by: {nm}")
        nodes.append(hits[0])
    lca = t.get_common_ancestor(nodes)
    t.set_outgroup(lca)
    return nodes

def force_outgroup_first(t, outgroup):
    if not outgroup:
        return
    root = t.get_tree_root()
    try:
        i = root.children.index(outgroup)
    except Exception as e:
        sys.stderr.write(f"[WARN] Cannot place outgroup first at root: {e}\n")
        return
    if i != 0:
        root.children.insert(0, root.children.pop(i))

def load_aliases(arg: str) -> dict:
    """Load species aliases from file or inline string."""
    if arg is None:
        return {}
    if os.path.isfile(arg):
        alias_dict = {}
        with open(arg, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) == 2:
                    species, alias = parts
                    alias_dict[species] = alias
        return alias_dict
    alias_dict = {}
    tokens = arg.strip().split()
    if all(":" in t for t in tokens):
        for t in tokens:
            species, alias = t.split(":", 1)
            alias_dict[species] = alias
        return alias_dict
    sys.stderr.write(f"[WARNING] --aliases not recognized as file or inline mapping: {arg}\n")
    return {}

def is_dup(node):
    """Return True for *any* duplication event (single‑ or multi‑species)."""
    return str(getattr(node, "D", "")).upper() == "Y"

def get_support(node):
    if not node.is_leaf() and node.name:
        support = node.name.strip()
        if support.replace('.', '', 1).isdigit():
            return float(support)
    return None

def is_low_support(node, tau):
    support = get_support(node)
    return (support is None) or (support < tau)

def is_boundary(node, tau):
    fb = getattr(node, "ForceBoundary", None)
    if fb is not None:          # explicit override wins
        return fb
    if is_dup(node) or (not node.is_leaf() and is_low_support(node, tau)):
        if IGNORE_SCO:
            species = [get_species_alias(leaf.S) for leaf in node.iter_leaves() if hasattr(leaf, "S")]
            if len(species) == len(set(species)):
                return False
        return True
    return False

def propagate_down_until_boundary(start, attr, value, tau):
    for child in start.children:
        if is_boundary(child, tau):
            continue
        if not hasattr(child, attr):
            setattr(child, attr, value)
        propagate_down_until_boundary(child, attr, value, tau)

def get_species_alias(name: str) -> str:
    return ALIASES.get(name, name)

def assign_node_ids(t):
    """Deterministic NodeID in preorder after rooting/sorting."""
    for i, n in enumerate(t.traverse("preorder"), start=1):
        n.NodeID = i

def _parse_pg_ignore_nodes(tokens):
    """
    tokens: list[str] from --PGignore-nodes.
    Accepts: "1 2 3", "1,2,3", mix, or file via "@file" or plain path.
    Returns: set[int]
    """
    if not tokens:
        return set()
    acc = []
    for tok in tokens:
        path = tok[1:] if tok.startswith("@") else tok
        if os.path.isfile(path):
            with open(path) as f:
                for line in f:
                    for part in line.replace(",", " ").split():
                        if part.strip().isdigit():
                            acc.append(int(part.strip()))
            continue
        for part in tok.replace(",", " ").split():
            if part.strip().isdigit():
                acc.append(int(part.strip()))
    return set(acc)

def annotate_leaf_counts(t):
    """Annotate each node with _nleaves = number of descendant leaves."""
    for n in t.traverse("postorder"):
        if n.is_leaf():
            n._nleaves = 1
        else:
            n._nleaves = sum(getattr(c, "_nleaves", 1 if c.is_leaf() else 0) for c in n.children)

def parse_id_list(s: str) -> set[int]:
    s = (s or "").strip()
    if not s:
        return set()
    parts = re.split(r"[,\s]+", s)
    out = set()
    for p in parts:
        if not p:
            continue
        try:
            out.add(int(p))
        except ValueError:
            sys.stderr.write(f"[WARN] Ignoring non-integer NodeID '{p}' in --force-pg-root\n")
    return out

def PG_mark_parent_groups(t, ignore_ids=None, max_pg_size=0, tau=0.0):
    """
    First non-dup in preorder opens a PG unless:
      - NodeID in ignore_ids, or
      - subtree leaf count > max_pg_size (when >0).
    Entire subtree gets ParentGroupNr.
    """
    ignore_ids = ignore_ids or set()
    pg = 0
    for n in t.traverse("preorder"):
        if getattr(n, "ParentGroupNr", None) is not None:
            continue
        if any(getattr(d, "ForcePGRoot", False) for d in n.iter_descendants()):
            continue
        if getattr(n, "ForcePGRoot", False):
            pg += 1
            for d in n.traverse():
                d.ParentGroupNr = pg
            continue
        if is_boundary(n, tau):
            continue
        if getattr(n, "NodeID", None) in ignore_ids:
            continue
        if max_pg_size and max_pg_size > 0:
            size = getattr(n, "_nleaves", None)
            if size is None:
                size = sum(1 for _ in n.iter_leaves())
            if size > max_pg_size:
                continue
        pg += 1
        for d in n.traverse():
            d.ParentGroupNr = pg
    return pg

def PG_assign_paranums_from_localnum(t):
    """
    Inside each PG, use distinct LocalNum values (OrthoCore proxies) as paralog groups.
    ParaNr order = ascending LocalNum within each PG.
    Returns dict: pg -> bool(has_multiple_paras).
    """
    pg_vals = sorted({getattr(leaf, "ParentGroupNr", None)
                      for leaf in t.iter_leaves() if getattr(leaf, "ParentGroupNr", None) is not None})
    has_multi = {}
    for pg in pg_vals:
        lns = sorted({getattr(leaf, "LocalNum", None)
                      for leaf in t.iter_leaves()
                      if getattr(leaf, "ParentGroupNr", None) == pg and getattr(leaf, "LocalNum", None) is not None})
        ln_to_para = {ln: i+1 for i, ln in enumerate(lns)}
        has_multi[pg] = len(lns) > 1
        for leaf in t.iter_leaves():
            if getattr(leaf, "ParentGroupNr", None) != pg:
                continue
            ln = getattr(leaf, "LocalNum", None)
            pn = ln_to_para.get(ln)
            if pn is not None:
                leaf.ParaNr = pn
    return has_multi

def _parse_pgpara_cell(x):
    """Parse canonical or display-style PG cells.
    Accepts interchangeable bare/suffixed singleton notation:
      - '25a'   -> ('25a', 1)
      - '25a-3' -> ('25a', 3)
      - '25'    -> (25, 1)
      - '25-3'  -> (25, 3)
    Legacy PGs are atomic locus names; trailing letters are part of PG, not Para.
    Return None on NaN/empty/malformed.
    """
    if pd.isna(x):
        return None
    s = str(x).strip()
    if not s:
        return None
    m = re.match(r"^\s*(\d+[A-Za-z]+)\s*-\s*(\d+)\s*$", s)
    if m:
        return m.group(1), int(m.group(2))
    m = re.match(r"^\s*(\d+[A-Za-z]+)\s*$", s)
    if m:
        return m.group(1), 1
    m = re.match(r"^\s*(\d+)\s*-\s*(\d+)\s*$", s)
    if m:
        return int(m.group(1)), int(m.group(2))
    m = re.match(r"^\s*(\d+)\s*$", s)
    if m:
        return int(m.group(1)), 1
    return None

def _load_revlog(path: str) -> pd.DataFrame:
    """Load CSV/Parquet revision table. Empty DF if file missing."""
    if not os.path.exists(path):
        return pd.DataFrame()
    if path.lower().endswith((".parquet", ".pq")):
        return pd.read_parquet(path)
    return pd.read_csv(path, index_col=0)

def _save_revlog(df: pd.DataFrame, path: str) -> None:
    if path.lower().endswith((".parquet", ".pq")):
        df.to_parquet(path, index=True)
    else:
        df.to_csv(path, index=True)

def _legacy_numeric_block_max(df: pd.DataFrame) -> int:
    """Largest pure numeric part seen in any PG token, including legacy strings like '25a'."""
    mx = 0
    for col in df.columns:
        ser = df[col].dropna().astype(str)
        for s in ser:
            parsed = _parse_pgpara_cell(s)
            if not parsed:
                continue
            pg = parsed[0]
            if isinstance(pg, int):
                mx = max(mx, pg)
            elif isinstance(pg, str):
                m = re.match(r"^(\d+)[A-Za-z]+$", pg)
                if m:
                    mx = max(mx, int(m.group(1)))
    return mx

def _legacy_string_override(pg_final: dict, A_cur: dict,
                             R: pd.DataFrame, cur_uid_set: set) -> dict:
    """Post-revival pass: if a temp PG was assigned integer I but the newest
    revlog column has a legacy string PG whose numeric prefix == I for any uid
    in that temp PG, the legacy string overrides the integer.

    Fixes two cases:
    (a) New ortholog labelled '25a' lands in the same tree clade as genes that
        had integer PG 25 -> revival picks integer 25 (oldest wins); override
        corrects to '25a'.
    (b) Genes manually re-labelled from '25' to '25a' in a newer revision;
        oldest-wins revival would keep integer 25; override corrects to '25a'.

    R.columns[0] is always the newest column (inserted via R.insert(0,...)).
    """
    if R.empty or len(R.columns) == 0:
        return pg_final
    newest_col = R.columns[0]
    uid_to_legacy: dict = {}
    for uid, val in R[newest_col].items():
        parsed = _parse_pgpara_cell(val)
        if parsed is not None and isinstance(parsed[0], str):
            uid_to_legacy[uid] = parsed[0]
    result = dict(pg_final)
    for tmp, uid_set in A_cur.items():
        cur_assign = result.get(tmp)
        if not isinstance(cur_assign, int):
            continue
        candidates: dict = defaultdict(set)
        for uid in uid_set:
            if uid in uid_to_legacy:
                candidates[uid_to_legacy[uid]].add(uid)
        for legacy_key in candidates:
            m = re.match(r"^(\d+)", legacy_key)
            if m and int(m.group(1)) == cur_assign:
                result[tmp] = legacy_key
                break
    return result

def _max_pg_and_para(df: pd.DataFrame):
    """Scan all columns to get global max *numeric* PG and per-PG max Para.
    Legacy string PGs do not contribute to pg_max, but do keep their own para_max entry.
    """
    pg_max = 0
    para_max = defaultdict(int)
    for col in df.columns:
        ser = df[col].dropna().astype(str)
        for s in ser:
            parsed = _parse_pgpara_cell(s)
            if not parsed:
                continue
            pg, pa = parsed
            if isinstance(pg, int) and pg > pg_max:
                pg_max = pg
            if pa is not None and pa > para_max[pg]:
                para_max[pg] = pa
    return pg_max, para_max

def _current_oc_sets_for_pg(t: Tree, pg_tmp: int) -> dict:
    """LocalNum groups inside a current temporary PG."""
    d = defaultdict(set)
    for leaf in t.iter_leaves():
        if getattr(leaf, "ParentGroupNr", None) == pg_tmp:
            ln = getattr(leaf, "LocalNum", None)
            if ln is not None:
                d[ln].add(leaf.name)
    return d

def _para_sets_from_prev_for_pg(R: pd.DataFrame, prev_col: str, pg_old) -> dict:
    """Map old Para -> set(UID) for a given old PG number using one historical column.
    Bare singleton forms like '25' or '25a' are normalized to Para 1.
    """
    d = defaultdict(set)
    if prev_col is None:
        return d
    ser = R[prev_col] if prev_col in R.columns else pd.Series(dtype=object)
    for uid, val in ser.items():
        parsed = _parse_pgpara_cell(val)
        if parsed and parsed[0] == pg_old:
            d[parsed[1]].add(uid)
    return d

def _old_pg_sets_from_col(R: pd.DataFrame, col: str, cur_uid_set: set[str]):
    """Return dict old_PG -> set(UID) for one historical column, intersected with current UIDs."""
    B = defaultdict(set)
    for uid, val in R[col].items():
        parsed = _parse_pgpara_cell(val)
        if parsed:
            b = parsed[0]
            if uid in cur_uid_set:
                B[b].add(uid)
    return B

def _history_columns_for_pg(R: pd.DataFrame, pg) -> list[str]:
    """All columns where this PG appears, ordered oldest→newest."""
    cols = []
    for col in list(R.columns)[::-1]:  # oldest first
        found = False
        for val in R[col].dropna():
            parsed = _parse_pgpara_cell(val)
            if parsed and parsed[0] == pg:
                found = True
                break
        if found:
            cols.append(col)
    return cols

def _revive_units(history_cols: list[str],
                  old_units_fn,
                  new_units: dict) -> dict:
    """
    Generic 1↔1 revival: oldest→newest; returns mapping new_unit -> old_label.
    old_units_fn(col) -> dict[old_label] = set(UIDs)
    new_units = dict[new_label] = set(UIDs)
    """
    assigned = {}
    taken_old = set()
    for col in history_cols:
        U_old = old_units_fn(col)
        recipients = {lab: set() for lab in U_old}
        for lab, S_old in U_old.items():
            for new_lab, S_new in new_units.items():
                if S_old & S_new:
                    recipients[lab].add(new_lab)
        oldcontrib = {new_lab: set() for new_lab in new_units}
        for new_lab, S_new in new_units.items():
            for lab, S_old in U_old.items():
                if S_old & S_new:
                    oldcontrib[new_lab].add(lab)
        for lab, rec in recipients.items():
            if lab in taken_old or len(rec) != 1:
                continue
            new_lab = next(iter(rec))
            if new_lab in assigned:
                continue
            if oldcontrib.get(new_lab, set()) == {lab}:
                assigned[new_lab] = lab
                taken_old.add(lab)
    return assigned

########################################################################################################################
# Parse arguments
########################################################################################################################
ap = argparse.ArgumentParser()
ap.add_argument("-i", "--infile", default="-")
ap.add_argument("-o", "--outfile", default="-")
ap.add_argument("-p", "--prefix", default="_",
                help="Gene family prefix used in ID construction (e.g. 'IR', 'GR', 'OR'). Defaults to '_' if omitted.")
ap.add_argument("-t", "--threshold", type=float, default=70.0,
                help="Minimum node support for boundary detection; float for posteriors or percentages")
ap.add_argument("-a", "--aliases",
                help="Either path to TSV file (SpeciesAlias) or inline mapping 'Spec1:Alias1 Spec2:Alias2'.")
ap.add_argument("--ignore-sco", dest="ignore_sco", action="store_true", default=True,
                help="Ignore low support within single-copy orthologs (default: True)")
ap.add_argument("--no-ignore-sco", dest="ignore_sco", action="store_false",
                help="Disable SCO-ignore behavior")
ap.add_argument("--minimal-id", type=int, default=1,
                help="Lowest numeric ID allowed in fresh assignment mode")
ap.add_argument("--include-outgroup", action="store_false", dest="skip_outgroup", default=True,
                help="Process the outgroup as well (default: skipped)")
ap.add_argument("--no-sort-tree", dest="sort_tree", action="store_false", default=True,
                help="Skip all pre-processing sorts (ladderize/core-order)")
ap.add_argument("--root-by", nargs="+", metavar="TIP",
                help="Root at the LCA of the listed tips (applies before outgroup handling)")
ap.add_argument("--PGignore-nodes", nargs="+",
                help="NodeIDs to ignore as PG roots (space- or comma-separated). "
                     "Also supports files prefixed with @ or plain paths.")
ap.add_argument("--max-pg-size", type=int, default=0,
                help="Max leaves allowed for a PG subtree. "
                     "If a candidate PG root exceeds this, it is skipped. 0 = no limit.")
ap.add_argument("--revlog", required=True, type=str,
                help="Path to the revision log table (CSV or Parquet). Index = immutable gene UID; columns = revisions. "
                     "The new column named by --revtag is inserted at position 0.")
ap.add_argument("--revtag", required=True, type=str,
                help="Name/label for the current revision column (e.g., 'rev_20251031').")
ap.add_argument("--no-allow-revive-pg", dest="allow_revive_pg", action="store_false", default=True,
                help="Disable revival of deprecated PG numbers if their exact sets reappear. Default is to allow revival.")
ap.add_argument("--temp-partition-out", type=str,
                help="Optional path to write the temporary '{PGtemp}-{LocalNum}' mapping for this run. "
                     "If omitted, no temp file is written.")
ap.add_argument("--force-pg-root", default="",
                help="Comma-separated NodeID list to force as PG roots, e.g. '56,86,184,490'.")
ap.add_argument("--show-single-para", dest="omit_single_para",
                action="store_false", default=True,
                help="Always append -ParaNr suffix, even for PG singletons. "
                     "Default: omit suffix when PG has only one paralog.")
ap.add_argument("--force-boundary-true", default="",
                help="Comma-separated NodeIDs to force as boundary=TRUE  regardless of support/dup status.")
ap.add_argument("--force-boundary-false", default="",
                help="Comma-separated NodeIDs to force as boundary=FALSE regardless of support/dup status.")

args = ap.parse_args()
if args.outfile == "-" and args.revtag:
    args.outfile = f"{args.revtag}.nwk"

########################################################################################################################
# Part 1 : Parse input
########################################################################################################################
TAU = args.threshold
PREFIX = args.prefix
IGNORE_SCO = args.ignore_sco
MINIMAL_ID_NUM = args.minimal_id
ALIASES = load_aliases(args.aliases)

if args.infile == "-":
    raw = sys.stdin.read()
else:
    try:
        with open(args.infile, "r", encoding="utf-8") as f:
            raw = f.read()
    except Exception as e:
        sys.stderr.write(f"[WARN] Failed to read input tree from '{args.infile}': {e}\n")
        raise
tree = Tree(raw, format=1)

if args.root_by:
    try:
        rooted_nodes = reroot_by_tips(tree, args.root_by)
    except Exception as e:
        sys.stderr.write(f"[ERROR] --root-by failed: {e}\n")
        sys.exit(2)

root = tree.get_tree_root()
outgroup = None
maingroup = tree

if args.skip_outgroup:
    if args.root_by:
        nodes = []
        lca = tree.get_common_ancestor(rooted_nodes)
        root = tree.get_tree_root()
        if len(root.children) == 2:
            c1, c2 = root.children
            outgroup = c1 if (lca is c1 or lca in c1.iter_descendants()) else c2
            maingroup = c2 if outgroup is c1 else c1
        else:
            outgroup = lca
            maingroup = root
        print("The following tips were selected as outgroup and will be skipped:", file=sys.stderr)
        for leaf in outgroup.iter_leaves():
            print(f"  {leaf.name}", file=sys.stderr)
        print("The following tips were selected as outgroup and will be skipped:", file=sys.stderr)
        for leaf in outgroup.iter_leaves():
            print(f"  {leaf.name}", file=sys.stderr)
    else:
        if len(root.children) == 2:
            c1, c2 = root.children
            n1 = sum(1 for _ in c1.iter_leaves())
            n2 = sum(1 for _ in c2.iter_leaves())
            outgroup = c1 if n1 < n2 else c2
            maingroup = c2 if outgroup is c1 else c1
            print("The following tips were selected as outgroup and will be skipped:", file=sys.stderr)
            for leaf in outgroup.iter_leaves():
                print(f"  {leaf.name}", file=sys.stderr)
        else:
            sys.stderr.write("[WARN] Root is not bifurcating; cannot determine outgroup reliably. Skipping outgroup removal.\n")
            outgroup = None

if args.skip_outgroup and outgroup is None:
    sys.stderr.write("[WARN] Outgroup skipping requested, but no outgroup was identified; using full tree.\n")
target_tree = maingroup if (args.skip_outgroup and outgroup is not None) else tree
if args.sort_tree:
    tree.ladderize()
    try:
        force_outgroup_first(tree, outgroup)
    except Exception as e:
        sys.stderr.write(f"[WARN] Failed to enforce outgroup-first ordering: {e}\n")

try:
    annotate_leaf_counts(target_tree)
except Exception as e:
    sys.stderr.write(f"[WARN] annotate_leaf_counts failed: {e}\n")

########################################################################################################################
# PART 2 : unified OrthoCore components (no SubGroups)
########################################################################################################################
assign_node_ids(tree)

PG_IGNORE_IDS = _parse_pg_ignore_nodes(args.PGignore_nodes) if hasattr(args, "PGignore_nodes") else set()

FORCE_BOUNDARY_TRUE  = parse_id_list(args.force_boundary_true)
FORCE_BOUNDARY_FALSE = parse_id_list(args.force_boundary_false)

for n in tree.traverse():
    nid = getattr(n, "NodeID", None)
    if nid in FORCE_BOUNDARY_TRUE:
        n.ForceBoundary = True
    elif nid in FORCE_BOUNDARY_FALSE:
        n.ForceBoundary = False

oc_counter = 1
for n in target_tree.traverse("preorder"):
    if is_boundary(n, TAU):
        continue
    if n.up is None or is_boundary(n.up, TAU):
        n.OrthoCoreNr = oc_counter
        oc_counter += 1

FORCE_IDS = parse_id_list(args.force_pg_root)
if FORCE_IDS:
    for n in target_tree.traverse():
        if getattr(n, "NodeID", None) in FORCE_IDS:
            n.ForcePGRoot = True
PG_mark_parent_groups(target_tree, ignore_ids=PG_IGNORE_IDS, max_pg_size=max(0, args.max_pg_size), tau=TAU)

########################################################################################################################
# Part 3 : LocalNum from boundary distances (single tree traversal)
########################################################################################################################
counter = 0
LocalNums = []

for node in target_tree.traverse("preorder"):
    if hasattr(node, "LocalNum"):
        continue
    if not hasattr(node, "OrthoCoreNr"):
        continue
    counter += 1
    node.LocalNum = counter
    propagate_down_until_boundary(node, "LocalNum", counter, TAU)
    LocalNums.append(counter)

PG_assign_paranums_from_localnum(target_tree)

if args.temp_partition_out:
    rows = []
    for leaf in target_tree.iter_leaves():
        rows.append({
            "Name": leaf.name,
            "PGtemp": getattr(leaf, "ParentGroupNr", None),
            "LocalNum": getattr(leaf, "LocalNum", None),
            "PGtemp-LocalNum": f"{getattr(leaf,'ParentGroupNr', '')}-{getattr(leaf,'LocalNum','')}"
        })
    pd.DataFrame(rows).to_csv(args.temp_partition_out, index=False)

########################################################################################################################
# Part 4 : REVISION RESOLUTION (revival oldest→newest; strict sets; per-PG Para revival)
########################################################################################################################

R = _load_revlog(args.revlog)

current_uids = [leaf.name for leaf in target_tree.iter_leaves()]
cur_uid_set = set(current_uids)

if R.empty:
    R = pd.DataFrame(index=current_uids)
else:
    missing = [u for u in current_uids if u not in R.index]
    if missing:
        R = pd.concat([R, pd.DataFrame(index=missing)], axis=0)

PG_MAX_EVER, PARA_MAX_EVER = _max_pg_and_para(R)
LEGACY_NUMERIC_BLOCK_MAX = _legacy_numeric_block_max(R)
NEXT_PG = max(MINIMAL_ID_NUM, LEGACY_NUMERIC_BLOCK_MAX + 1, PG_MAX_EVER + 1)

A_cur = defaultdict(set)
for leaf in target_tree.iter_leaves():
    pg_tmp = getattr(leaf, "ParentGroupNr", None)
    if pg_tmp is not None:
        A_cur[pg_tmp].add(leaf.name)

cols_pg = list(R.columns)[::-1] if (args.allow_revive_pg and len(R.columns) > 0) else ([R.columns[0]] if len(R.columns) > 0 else [])

def _old_units_fn_pg(col):
    return _old_pg_sets_from_col(R, col, cur_uid_set)

pg_final = {}
pg_assigned = _revive_units(cols_pg, _old_units_fn_pg, A_cur)
for pg_tmp, old_pg in pg_assigned.items():
    pg_final[pg_tmp] = old_pg

for pg_tmp in A_cur.keys():
    if pg_tmp not in pg_final:
        pg_final[pg_tmp] = NEXT_PG
        NEXT_PG += 1

# Legacy-string override: newest revlog col's legacy PG supersedes any integer
# assignment that shares the same numeric prefix (e.g. '25a' overrides int 25).
pg_final = _legacy_string_override(pg_final, A_cur, R, cur_uid_set)

# Intra-PG Para resolution
for pg_tmp, final_pg in pg_final.items():
    oc_new = _current_oc_sets_for_pg(target_tree, pg_tmp)
    kept = pg_assigned.get(pg_tmp) == final_pg

    if kept and len(R.columns) > 0:
        hist_cols = _history_columns_for_pg(R, final_pg)

        def _old_units_fn_para(col):
            return _para_sets_from_prev_for_pg(R, col, final_pg)

        reuse_map_ln_to_para = _revive_units(hist_cols, _old_units_fn_para, oc_new)

        para_next = PARA_MAX_EVER.get(final_pg, 0) + 1
        for u in oc_new.keys():
            if u not in reuse_map_ln_to_para:
                reuse_map_ln_to_para[u] = para_next
                para_next += 1

        for leaf in target_tree.iter_leaves():
            if getattr(leaf, "ParentGroupNr", None) == pg_tmp:
                ln = getattr(leaf, "LocalNum", None)
                if ln is not None:
                    leaf.ParaNr = reuse_map_ln_to_para.get(ln)

        if reuse_map_ln_to_para:
            PARA_MAX_EVER[final_pg] = max(PARA_MAX_EVER.get(final_pg, 0),
                                          max(reuse_map_ln_to_para.values()))

for leaf in target_tree.iter_leaves():
    pg_tmp = getattr(leaf, "ParentGroupNr", None)
    if pg_tmp is not None:
        leaf.ParentGroupNr = pg_final[pg_tmp]

# Compose current revision column in canonical storage form.
# Bare singleton forms in historical revlogs are accepted on input, but all newly
# written assignments are normalized to explicit PG-Para strings.
cur_series = pd.Series(index=R.index, dtype=object)
for uid in R.index:
    cur_series.at[uid] = pd.NA
for leaf in target_tree.iter_leaves():
    uid = leaf.name
    pg = getattr(leaf, "ParentGroupNr", None)
    pa = getattr(leaf, "ParaNr", None)
    if pg is None:
        continue
    pa = 1 if pa is None else int(pa)
    if isinstance(pg, str):
        cur_series.at[uid] = f"{pg}-{pa}"
    else:
        cur_series.at[uid] = f"{int(pg)}-{pa}"

R.insert(0, args.revtag, cur_series)
revlog_out = f"revlog_{args.revtag}.csv"
_save_revlog(R, revlog_out)
sys.stderr.write(f"Revlog written to {revlog_out}\n")

########################################################################################################################
# PART 5 : Create and apply gene-IDs
########################################################################################################################
# Find PGs where all members have ParaNr == 1 (no within-PG paralogs)
pg_max_para = defaultdict(int)
for leaf in target_tree.iter_leaves():
    pg = getattr(leaf, "ParentGroupNr", None)
    pa = getattr(leaf, "ParaNr", None)
    if pg is not None and pa is not None:
        pg_max_para[pg] = max(pg_max_para[pg], int(pa))

for leaf in target_tree.iter_leaves():
    pg = getattr(leaf, "ParentGroupNr", None)
    if pg is None:
        sys.stderr.write(f"[WARNING] Missing ParentGroupNr for {leaf.name}; placeholder ID set.\n")
        leaf.add_features(geneID="XXXXXXXXXXX")
        continue
    alias = get_species_alias(getattr(leaf, "S", ""))
    base = f"{alias}{PREFIX}{pg}"
    pa = getattr(leaf, "ParaNr", None)
    pa = 1 if pa is None else int(pa)
    if args.omit_single_para and pg_max_para.get(pg, 1) == 1:
        leaf.geneID = base
    else:
        leaf.geneID = f"{base}-{pa}"
    leaf.IDNr = f"{pg}-{pa}"

########################################################################################################################
# Part 6 : SAVE FINAL TREE
########################################################################################################################
features_to_include = ["NodeID","S","D","H","B","OrthoCoreNr","ParentGroupNr","ParaNr","geneID","LocalNum","IDNr"]

out=tree.write(format=1, features=features_to_include)
if args.outfile == "-": sys.stdout.write(out)
else:
    with open(args.outfile,"w") as fh: fh.write(out)
print("✅ Named tree written to", args.outfile, file=sys.stderr)
if args.outfile != "-":
    csv_path = os.path.splitext(args.outfile)[0] + ".csv"
    with open(csv_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Name", "geneID", "ParentGroupNr", "ParaNr"])
        for leaf in tree.iter_leaves():
            w.writerow([leaf.name,
                        getattr(leaf, "geneID", ""),
                        getattr(leaf, "ParentGroupNr", ""),
                        getattr(leaf, "ParaNr", "")])
    print(f"Wrote CSV: {csv_path}", file=sys.stderr)