process INFER_GENE_TREE {
    label 'gene_tree'
    tag 'IQ-TREE gene-tree inference'
    cpus { params.process_cpus.INFER_GENE_TREE }
    input:
    path trimmed_alignment
    val iqtree_args
    output:
    path 'gene_tree.treefile', emit: gene_tree
    path 'gene_tree.iqtree', emit: report, optional: true
    path 'gene_tree.log', emit: log, optional: true
    path 'gene_tree.model.gz', emit: model, optional: true
    script:
    def encoded_args = (iqtree_args ?: '').toString().getBytes('UTF-8').encodeBase64().toString()
    """
    python3 - '${trimmed_alignment}' '${encoded_args}' '${params.iqtree_threads}' '${task.cpus}' <<'PY'
import base64, re, shlex, subprocess, sys
alignment, encoded, configured_threads, cpu = sys.argv[1:]; cpus=int(cpu)
tokens=shlex.split(base64.b64decode(encoded).decode()); overrides={}; extra=[]; i=0
aliases={'-mset':'-mset','-mrate':'-mrate','-B':'-B'}
forbidden={'-T','-nt','--threads','-ntmax','--threads-max'}
while i < len(tokens):
    token=tokens[i]; flag, equals, value=token.partition('=')
    if flag in forbidden: raise SystemExit('Thread options belong in params.iqtree_threads; -ntmax is set from task.cpus')
    if flag in aliases:
        if not equals:
            i+=1
            if i==len(tokens): raise SystemExit(f'Missing value for {flag}')
            value=tokens[i]
        overrides[aliases[flag]]=value
    else:
        if re.match(r'^-(?:T|nt|ntmax)[0-9]',token): raise SystemExit('Use the dedicated iqtree_threads setting')
        extra.append(token)
    i+=1
threads='AUTO' if configured_threads.upper()=='AUTO' else configured_threads
if threads!='AUTO' and (not threads.isdecimal() or not 1<=int(threads)<=cpus): raise SystemExit(f'iqtree_threads must be AUTO or 1..{cpus}')
command=['iqtree','-s',alignment,'-pre','gene_tree','-mset',overrides.get('-mset','raxml'),'-mrate',overrides.get('-mrate','G'),'-B',overrides.get('-B','1000'),'-T',threads,'-ntmax',str(cpus),*extra]
subprocess.run(command,check=True)
PY
    """
}
