"""Compare executable sections of trusted locally built ELF64 objects."""
import hashlib,json,struct
from pathlib import Path

def fingerprint(path):
    data=Path(path).read_bytes()
    if data[:6]!=b'\x7fELF\x02\x01':raise ValueError('expected little-endian ELF64')
    shoff=struct.unpack_from('<Q',data,40)[0]
    entsize,num,names_index=struct.unpack_from('<HHH',data,58)
    headers=[struct.unpack_from('<IIQQQQIIQQ',data,shoff+i*entsize) for i in range(num)]
    names_header=headers[names_index];names=data[names_header[4]:names_header[4]+names_header[5]]
    sections=[]
    for h in headers:
        if not h[2]&4:continue
        name=names[h[0]:].split(b'\0',1)[0].decode()
        sections.append((name,hashlib.sha256(data[h[4]:h[4]+h[5]]).hexdigest()))
    return sections

def write_report(out,sets):
    rows=[]
    for row in sets:
        for backend in ('ref','avx2'):
            selected=[next(v for v in row['variants'] if v['name']==name) for name in (backend,'vendored_'+backend,'vendored_'+backend+'_instrumented')]
            dirs=[Path(out)/'work'/(row['name']+'-'+v['name']) for v in selected]
            for i,source in enumerate(selected[0]['sources']):
                hashes=[fingerprint(d/(str(i)+'.o')) for d in dirs]
                rows.append(dict(set=row['name'],backend=backend,source=source,upstream_equals_vendor=hashes[0]==hashes[1],vendor_equals_instrumented=hashes[1]==hashes[2],executable_sections=hashes))
            for source in ('driver','cpucycles','speed_print'):
                hashes=[fingerprint(d/(source+'.o')) for d in dirs]
                rows.append(dict(set=row['name'],backend=backend,source=source+'.o',upstream_equals_vendor=hashes[0]==hashes[1],vendor_equals_instrumented=hashes[1]==hashes[2],executable_sections=hashes))
    (Path(out)/'source-object-comparison.json').write_text(json.dumps(rows,indent=2)+'\n')
    return rows
