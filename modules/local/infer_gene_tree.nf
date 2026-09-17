process INFER_GENE_TREE {
    label 'gene_tree'
    tag 'IQ-TREE gene-tree inference'

    input:
    path trimmed_alignment
    val iqtree_args

    output:
    path 'gene_tree.treefile', emit: gene_tree
    path 'gene_tree.iqtree', emit: report, optional: true
    path 'gene_tree.log', emit: log, optional: true
    path 'gene_tree.model.gz', emit: model, optional: true

    script:
    def args = iqtree_args ?: ''
    """
python3 - '${trimmed_alignment}' '${args}' <<'PY'
import shlex
import subprocess
import sys

alignment, user_arg_string = sys.argv[1:]
user_tokens = shlex.split(user_arg_string)
defaults = [
    ('-mset', 'raxml'),
    ('-mrate', 'G'),
    ('-B', '1000'),
    ('-T', 'AUTO'),
]
default_flags = {flag for flag, _ in defaults}
overrides = {}
extras = []
i = 0
while i < len(user_tokens):
    token = user_tokens[i]
    if token in default_flags:
        if i + 1 >= len(user_tokens):
            raise SystemExit(f'Missing value for IQ-TREE option: {token}')
        overrides[token] = user_tokens[i + 1]
        i += 2
    else:
        extras.append(token)
        i += 1

effective_args = []
for flag, value in defaults:
    effective_args.extend([flag, overrides.get(flag, value)])
effective_args.extend(extras)

subprocess.run(
    ['iqtree', '-s', alignment, '-pre', 'gene_tree', *effective_args],
    check=True,
)
PY
"""
}
