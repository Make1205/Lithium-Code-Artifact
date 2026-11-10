#!/usr/bin/env python3
import csv,datetime,hashlib,io,json,os,platform,re,shutil,subprocess,sys,tarfile
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
BENCH=ROOT/'bench'
CFG=json.loads((BENCH/'config.json').read_text())
ENV=dict(os.environ,LANG='C',LC_ALL='C',GIT_TERMINAL_PROMPT='0')
for k in ('CPATH','C_INCLUDE_PATH','CPLUS_INCLUDE_PATH','LIBRARY_PATH','LD_LIBRARY_PATH','GCC_EXEC_PREFIX','COMPILER_PATH'):ENV.pop(k,None)
CC='/usr/bin/gcc'
BASE=['-O3','-march=native','-std=gnu11']
def call(args,cwd=ROOT,log=None):
 p=subprocess.run([str(x) for x in args],cwd=cwd,env=ENV,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
 s=p.stdout.decode(errors='replace')
 if log:
  with log.open('a') as f:f.write('$ '+' '.join(map(str,args))+'\n'+s)
 if p.returncode:raise RuntimeError('exit %d: %s\n%s'%(p.returncode,args,s[-6000:]))
 return s

def archive(commit,path):
 data=subprocess.check_output(['git','archive',commit],cwd=ROOT,env=ENV)
 path.mkdir(parents=True)
 with tarfile.open(fileobj=io.BytesIO(data)) as tf:
  for m in tf.getmembers():
   target=(path/m.name).resolve()
   if os.path.commonpath([str(path.resolve()),str(target)])!=str(path.resolve()):raise RuntimeError('archive path escape')
   if m.issym() and os.path.commonpath([str(path.resolve()),str((target.parent/m.linkname).resolve())])!=str(path.resolve()):raise RuntimeError('archive symlink escape')
  tf.extractall(path)
 shutil.copytree(BENCH,path/'bench',ignore=shutil.ignore_patterns('results','__pycache__'),dirs_exist_ok=True)

def main():
 if sys.version_info<(3,8):raise RuntimeError('Python >=3.8 required: apt install python3')
 if platform.machine()!='x86_64':raise RuntimeError('RDTSC/AVX2 benchmark requires x86-64')
 if not Path(CC).is_file():raise RuntimeError('Missing GCC: apt install build-essential')
 if not shutil.which('git'):raise RuntimeError('Missing Git: apt install git')
 probe=Path('/tmp')/('bench-toolchain-'+str(os.getpid()))
 cp=subprocess.run([CC,'-x','c','-','-o',str(probe),'-lm'],input=b'#include <stdint.h>\n#include <stddef.h>\n#include <stdio.h>\n#include <stdlib.h>\nint main(void){void *p=malloc(8);free(p);return 0;}\n',env=ENV,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
 if cp.returncode:raise RuntimeError('System C toolchain/headers/libm missing: apt install build-essential (includes libc6-dev)\n'+cp.stderr.decode())
 probe.unlink()
 for dep in CFG.get('dependencies',[]):
  p=subprocess.run([CC,'-x','c','-','-o','/tmp/bench-dep-'+str(os.getpid())]+dep['libs'],input=dep['source'].encode(),stdout=subprocess.PIPE,stderr=subprocess.PIPE,env=ENV)
  if p.returncode:raise RuntimeError('Missing dependency: '+dep['package']+'; apt install '+dep['package']+'\n'+p.stderr.decode())
  Path('/tmp/bench-dep-'+str(os.getpid())).unlink()
 quick=os.environ.get('QUICK')=='1';nt=11 if quick else 10001;ns=100 if quick else 10000
 if hasattr(os,'sched_getaffinity'):cpu=min(os.sched_getaffinity(0));os.sched_setaffinity(0,{cpu})
 else:cpu=None
 stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S.%fZ');out=BENCH/'results'/stamp;out.mkdir(parents=True)
 print('Results:',out,flush=True)
 meta={'scale':'QUICK' if quick else 'FULL','timing_iterations':nt,'signatures':ns,'cpu':platform.processor(),'affinity':cpu,'gcc':call([CC,'--version']).splitlines()[0],'benchmark_commit':call(['git','rev-parse','HEAD']).strip(),'tested_commit':CFG['tested_commit'],'initial_commit':CFG['initial_commit'],'initial_status':call(['git','status','--porcelain']),'flags':BASE,'configuration':CFG,'benchmark_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in BENCH.iterdir() if p.is_file()}}
 if Path('/proc/cpuinfo').exists():meta['cpu']=re.search(r'^model name\s*:\s*(.*)$',Path('/proc/cpuinfo').read_text(),re.M)[1]
 (out/'metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
 stage=out/'work'/'source';archive(CFG['tested_commit'],stage)
 results=[];lengths=[];correct=[];builds={}
 for row in CFG['sets']:
  for variant in row['variants']:
   v=variant['name'];key=row['name']+'-'+v;print('Build/KAT:',key,flush=True)
   folder=stage/variant['cwd'];dest=out/'work'/key;dest.mkdir(parents=True)
   flags=BASE+variant.get('flags',[])+['-I.', '-I'+str(stage/'bench')]
   log=out/(key+'-build.log');objs=[]
   for j,s in enumerate(variant['sources']):
    obj=dest/(str(j)+'.o');call([CC]+flags+['-c',s,'-o',obj],folder,log);objs.append(obj)
   extra=[]
   for j,spec in enumerate(variant.get('extra_sources',[])):
    obj=dest/('extra'+str(j)+'.o');call([CC]+flags+spec.get('flags',[])+['-c',spec['path'],'-o',obj],folder,log);extra.append(obj)
   driver=dest/'driver';call([CC]+flags+['-include',str(stage/'bench'/variant.get('adapter','adapter.h')),'-c',str(stage/'bench/driver.c'),'-o',dest/'driver.o'],folder,log)
   timers=[]
   for name in ('cpucycles','speed_print'):
    obj=dest/(name+'.o');call([CC]+BASE+['-c',str(stage/'bench'/(name+'.c')),'-o',obj],folder,log);timers.append(obj)
   call([CC,str(dest/'driver.o')]+objs+extra+timers+variant.get('libs',[])+['-o',driver],folder,log)
   for kat in variant.get('kat',[]):
    if 'command' in kat:
     text=call(kat['command'],folder,out/'kat.log')
    else:
     exe=dest/kat['name'];call([CC]+flags+kat.get('flags',[])+kat['sources']+objs+kat.get('extra_objects',[])+variant.get('libs',[])+['-o',exe],folder,log)
     text=call([exe]+kat.get('args',[]),folder,out/'kat.log')
    if kat.get('expect') and kat['expect'] not in text:raise RuntimeError('KAT expected output missing: '+key)
   builds[row['name'],v]=driver
 # Cross checks and lengths precede all timing.
 for row in CFG['sets']:
  files={};variants=[v['name'] for v in row['variants']]
  for v in variants:
   key=row['name']+'-'+v;f=out/'work'/(key+'.vectors');files[v]=f
   text=call([builds[row['name'],v],'generate',ns,f],log=out/'correctness.log')
   n,q,pk,sk,sig=map(int,re.search(r'META (\d+) (-?\d+) (\d+) (\d+) (\d+)',text).groups())
   vals=list(map(int,re.search(r'LENGTH (\d+) (\d+) (\d+) (\d+)',text).groups()))
   total,count,lo,hi=vals
   lengths.append({'set':row['name'],'build':v,'n':n,'q':q,'pk':pk,'sk':sk,'sig':sig,'count':count,'mean':total/count,'min':lo,'max':hi})
  equal=None
  if len(variants)==2:
   equal=files[variants[0]].read_bytes()==files[variants[1]].read_bytes()
   if row.get('require_equal',False) and not equal:raise RuntimeError('Expected byte equality failed: '+row['name'])
  for v in variants:
   other=variants[-1] if v==variants[0] else variants[0]
   call([builds[row['name'],v],'check',ns,files[other]],log=out/'correctness.log')
  correct.append({'set':row['name'],'count_each':ns,'cross_builds':len(variants)==2,'byte_equal':equal})
  print('Correctness:',row['name'],correct[-1],flush=True)
 (out/'correctness.json').write_text(json.dumps(correct,indent=2)+'\n');(out/'lengths.json').write_text(json.dumps(lengths,indent=2)+'\n')
 for row in CFG['sets']:
  for variant in row['variants']:
   v=variant['name'];key=row['name']+'-'+v;print('Timing:',key,nt,flush=True)
   text=call([builds[row['name'],v],'time',nt],log=out/(key+'-timing.log'))
   for op,value in re.findall(r'^(KeyGen|Sign|Verify|ExpandedSign)\nmedian: (\d+) cycles/ticks',text,re.M):
    size=next(x for x in lengths if x['set']==row['name'] and x['build']==v)
    results.append(dict(size,operation=op,iterations=nt,median_cycles=int(value),kcycles=int(value)/1000))
   for op in ('KeyGen','Sign','Verify'):
    if not any(x['set']==row['name'] and x['build']==v and x['operation']==op for x in results):raise RuntimeError('Missing timing: '+key+' '+op)
 (out/'timings.json').write_text(json.dumps(results,indent=2)+'\n')
 with (out/'bench.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(results[0]));w.writeheader();w.writerows(results)
 lookup={(x['set'],x['build'],x['operation']):x['kcycles'] for x in results}
 s='# '+CFG['scheme']+' — '+meta['scale']+'\n\n'
 s+='QUICK is a workflow check only. Full-scale figures are measured locally.\n\n' if quick else 'Full-scale local benchmark.\n\n'
 s+='CPU: '+meta['cpu']+'; affinity: '+str(cpu)+'. GCC: '+meta['gcc']+'.\n\nTested commit: `'+CFG['tested_commit']+'`. Benchmark commit: `'+meta['benchmark_commit']+'`.\n\n'
 s+='Base flags: `'+ ' '.join(BASE)+'`. Each operation: '+str(nt)+' iterations; correctness and lengths: '+str(ns)+' signatures per build/set. Message: 32 bytes. Median units: thousand RDTSC ticks.\n\n'
 s+='Exact per-build options (in addition to base flags and local include paths):\n\n'
 for rr in CFG['sets']:
  for vv in rr['variants']:
   s+='* '+rr['name']+' '+vv['name']+': `'+ ' '.join(vv.get('flags',[]))+'`; link `'+ ' '.join(vv.get('libs',[]))+'`.\n'
 s+='\n'
 s+='| 参数集 | n | q | pk | sig（仓库实际） | 论文 sig | ref KeyGen/Sign/Verify | AVX2 KeyGen/Sign/Verify |'
 if CFG.get('calibration'):s+=' ratio to previous ref / AVX2 |'
 s+='\n|---|---:|---:|---:|---|---:|---|---|'+('---|' if CFG.get('calibration') else '')+'\n'
 mismatches=[]
 for row in CFG['sets']:
  name=row['name'];sizes=next(x for x in lengths if x['set']==name);trip=[];ratios=[]
  for v in ('ref','avx2'):
   trip.append('/'.join('%.3f'%lookup[name,v,op] for op in ('KeyGen','Sign','Verify')) if (name,v,'Sign') in lookup else '--')
   if CFG.get('calibration'):
    ratios.append('/'.join('%.3f'%(lookup[name,v,op]/base) for op,base in zip(('KeyGen','Sign','Verify'),CFG['calibration'][name][v])) if (name,v,'Sign') in lookup else '--')
  sigstr=str(sizes['sig'])+'; unpadded avg %.4f, max %d'%(sizes['mean'],sizes['max'])
  s+='| %s | %d | %s | %d | %s | %s | %s | %s |'%(name,sizes['n'],'--' if sizes['q']<0 else sizes['q'],sizes['pk'],sigstr,row.get('paper_sig','--'),*trip)
  if ratios:s+=' '+' ; '.join(ratios)+' |'
  s+='\n'
  if sizes['pk']!=row.get('paper_pk',sizes['pk']) or sizes['q']!=row.get('paper_q',sizes['q']) or sizes['sig']!=row.get('paper_sig',sizes['sig']):
   mismatches.append('Mismatch for '+name+': repository (q,pk,sig max/fixed)='+str((sizes['q'],sizes['pk'],sizes['sig']))+'; paper='+str((row.get('paper_q'),row.get('paper_pk'),row.get('paper_sig')))+'.')
 s+='\n'+'\n\n'.join(mismatches)+'\n'
 s+='\nExpanded-key Sign (thousand ticks; `--` means no supported expanded-key API):\n\n| 参数集 | ref | AVX2 |\n|---|---:|---:|\n'
 for row in CFG['sets']:
  name=row['name'];s+='| '+name+' | '+' | '.join('%.3f'%lookup[name,v,'ExpandedSign'] if (name,v,'ExpandedSign') in lookup else '--' for v in ('ref','avx2'))+' |\n'
 s+='\nCorrectness: all configured upstream KAT/self-tests passed; all generated signatures verified; changed messages rejected. Cross-build/equality details are in correctness.json.\n\n'
 s+=CFG.get('notes','')+'\n\nBuild flags and commands per variant are in config.json, metadata.json and build logs. All source is extracted locally from the tested commit; algorithm files are unchanged. See bench/NOTES.md.\n'
 (out/'SUMMARY.md').write_text(s)
 print('PASS:',out/'SUMMARY.md',flush=True)
if __name__=='__main__':
 try:main()
 except Exception as e:print('FAILED:',e,file=sys.stderr,flush=True);sys.exit(1)
