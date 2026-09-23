# Reconaming naming controls. Leave a value blank to use the core-script default.
# This file is intentionally Maker-like: only non-empty values take effect.

# Gene-family label inserted between species alias and ParentGroup number.
prefix =

# Minimum internal-node support required to avoid a naming boundary.
threshold =

# Path to a two-column, tab-separated species-to-alias file; alternatively use an inline mapping.
aliases =

# true: do not split low-support single-copy ortholog nodes; false: treat them as boundaries.
ignore_sco =

# Lowest numeric ParentGroup identifier for a newly created registry.
minimal_id =

# true: name the rooted outgroup as well; false/default: exclude it from naming.
include_outgroup =

# true/default: ladderize and core-order the tree; false: preserve original ordering.
sort_tree =

# Optional whitespace-separated reconciled-tree leaf IDs used to root the tree.
# When blank, non-comment entries in --outgroup_genes are used automatically.
root_by =

# Node IDs excluded as ParentGroup roots. Accepts a whitespace/comma list or a file path.
pgignore_nodes =

# Largest permitted ParentGroup subtree; 0/default means no limit.
max_pg_size =

# true/default: restore unambiguous historical ParentGroup assignments; false: disable restoration.
allow_revive_pg =

# Comma-separated Node IDs forced to be ParentGroup roots.
force_pg_root =

# true: append '-1' to singleton paralog groups; false/default: omit singleton suffixes.
show_single_para =

# Comma-separated Node IDs forced to be, respectively, boundaries or non-boundaries.
force_boundary_true =
force_boundary_false =
