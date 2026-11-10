"""Audit x86 FMA instructions in benchmark-owned objects and executables."""
import json, os, re, shutil, subprocess
from pathlib import Path
FMA=re.compile(r'\bv(?:f|fn)m(?:add|sub)[a-z0-9]*\b')
SYMBOL=re.compile(r'^\s*[0-9a-f]+ <([^>]+)>:')
def preflight():
    if not shutil.which('objdump'):raise RuntimeError('Missing objdump: apt install binutils')
def instructions(path):
    preflight()
    cp=subprocess.run(['objdump','-d',str(path)],stdout=subprocess.PIPE,stderr=subprocess.PIPE,env=dict(os.environ,LANG='C',LC_ALL='C'))
    if cp.returncode:raise RuntimeError('objdump failed: '+str(path)+'\n'+cp.stderr.decode(errors='replace'))
    symbol='';found=[]
    for line in cp.stdout.decode(errors='replace').splitlines():
        m=SYMBOL.match(line)
        if m:symbol=m[1]
        if FMA.search(line):found.append(dict(symbol=symbol,instruction=line.strip()))
    return found

def audit(executable, objects, report, explicit_objects=()):
    executable=Path(executable);allowed={Path(x).name for x in explicit_objects}
    records=[];explicit_symbols=set();unexpected=[]
    for path in objects:
        path=Path(path);hits=instructions(path);explicit=path.name in allowed
        records.append(dict(file=path.name,fma_count=len(hits),explicit_intrinsics_allowed=explicit,instructions=hits))
        if explicit:explicit_symbols.update(h['symbol'] for h in hits)
        elif hits:unexpected.append(path.name)
    linked=instructions(executable)
    unexpected_linked=[h for h in linked if h['symbol'] not in explicit_symbols]
    status='PASS' if not unexpected and not unexpected_linked else 'FAIL'
    data=dict(status=status,executable=executable.name,executable_fma_count=len(linked),
              explicit_intrinsics_allowed_objects=sorted(allowed),object_fma_counts=records,
              executable_instructions=linked,unexpected_objects=unexpected,
              unexpected_executable_instructions=unexpected_linked,
              note='Explicit upstream intrinsics are preserved; compiler contraction is disabled for floating-point builds. Shared-library code is outside this executable audit.')
    report=Path(report);report.parent.mkdir(parents=True,exist_ok=True);report.write_text(json.dumps(data,indent=2)+'\n')
    if status!='PASS':raise RuntimeError('Unexpected FMA instructions: '+str(report))
    return data
