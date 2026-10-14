#!/bin/bash -euo pipefail
python3 - 'representative_proteins.faa' '' '64' <<'PY'
import base64, shlex, subprocess, sys
input_file, encoded, cpu = sys.argv[1:]; cpus = int(cpu)
tokens = shlex.split(base64.b64decode(encoded).decode())
threads = None; extra = []; i = 0
while i < len(tokens):
    token = tokens[i]; flag, equals, value = token.partition('=')
    if flag == '-t':
        if threads is not None: raise SystemExit('Duplicate FAMSA -t option')
        if not equals:
            i += 1
            if i == len(tokens): raise SystemExit('Missing FAMSA -t value')
            value = tokens[i]
        if not value.isdecimal() or not 1 <= int(value) <= cpus: raise SystemExit(f'FAMSA -t must be 1..{cpus}')
        threads = int(value)
    else:
        if token.startswith('-t') and token[2:].isdecimal(): raise SystemExit('Use -t N or -t=N')
        extra.append(token)
    i += 1
subprocess.run(['famsa','-t',str(threads or cpus),*extra,input_file,'representative_proteins.famsa.faa'],check=True)
PY
