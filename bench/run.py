#!/usr/bin/env python3
import csv,datetime,hashlib,io,itertools,json,os,platform,re,shutil,subprocess,sys,tarfile
from pathlib import Path
import fma_audit
from source_text_audit import write_report
ROOT=Path(__file__).resolve().parent.parent
BENCH=ROOT/'bench'
CFG=json.loads((BENCH/'config.json').read_text())
ENV=dict(os.environ,LANG='C',LC_ALL='C',GIT_TERMINAL_PROMPT='0')
for k in ('CPATH','C_INCLUDE_PATH','CPLUS_INCLUDE_PATH','LIBRARY_PATH','LD_LIBRARY_PATH','GCC_EXEC_PREFIX','COMPILER_PATH'):ENV.pop(k,None)
CC='/usr/bin/gcc'
if CFG.get('timing_seed'):ENV.setdefault('BENCH_SEED',CFG['timing_seed'])
if CFG.get('default_modes'):
 selected=[int(x) for x in os.environ.get('COPA_MODES',','.join(map(str,CFG['default_modes']))).split(',')]
 if not selected or len(set(selected))!=len(selected) or any(x not in [r['mode'] for r in CFG['sets']] for x in selected):raise RuntimeError('Unsupported COPA_MODES')
 CFG['sets']=[r for mode in selected for r in CFG['sets'] if r['mode']==mode]
BASE=['-O3','-march=native','-std=c99','-ffp-contract=off']
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
 fma_audit.preflight()
 if CFG.get('source_attribution'):
  vendor=BENCH/'vendor/copa-3c676bd'
  manifest=json.loads((vendor/'MANIFEST.json').read_text())
  for name,item in manifest['files'].items():
   if hashlib.sha256((vendor/name).read_bytes()).hexdigest()!=item['sha256']:raise RuntimeError('Vendored source checksum mismatch: '+name)
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
 meta={'scale':'QUICK' if quick else 'FULL','timing_iterations':nt,'signatures':ns,'cpu':platform.processor(),'affinity':cpu,'gcc':call([CC,'--version']).splitlines()[0],'benchmark_commit':call(['git','rev-parse','HEAD']).strip(),'tested_commit':CFG['tested_commit'],'initial_commit':CFG['initial_commit'],'initial_status':call(['git','status','--porcelain']),'timing_warmup':20,'timing_key_population':'fresh key/signature per triple','timing_entropy':'SHAKE(seed,counter)' if ENV.get('BENCH_SEED') is not None else 'native OS','calibration_requested':os.environ.get('CALIBRATION')=='1','benchmark_seed':ENV.get('BENCH_SEED'),'flags':BASE,'configuration':CFG,'benchmark_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in BENCH.iterdir() if p.is_file()}}
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
   fma_audit.audit(driver,[dest/'driver.o']+objs+extra+timers,out/(key+'-fma-audit.json'))
   for kat in variant.get('kat',[]):
    if 'command' in kat:
     text=call(kat['command'],folder,out/'kat.log')
    else:
     exe=dest/kat['name'];call([CC]+flags+kat.get('flags',[])+kat['sources']+objs+kat.get('extra_objects',[])+variant.get('libs',[])+['-o',exe],folder,log)
     text=call([exe]+kat.get('args',[]),folder,out/'kat.log')
    if kat.get('expect') and kat['expect'] not in text:raise RuntimeError('KAT expected output missing: '+key)
    if kat.get('expected_file'):
     actual=text.encode();expected=(stage/kat['expected_file']).read_bytes()
     (out/(key+'-kat.rsp')).write_bytes(actual)
     proof=dict(actual_sha256=hashlib.sha256(actual).hexdigest(),expected_sha256=hashlib.sha256(expected).hexdigest(),fixture=kat['expected_file'],byte_equal=actual==expected)
     (out/(key+'-kat-equality.json')).write_text(json.dumps(proof,indent=2)+'\n')
     if actual!=expected:raise RuntimeError('Native KAT byte mismatch: '+key)

   builds[row['name'],v]=driver
 if CFG.get('source_attribution'):write_report(out,CFG['sets'])
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
  equal=all(files[v].read_bytes()==files[variants[0]].read_bytes() for v in variants[1:])
  for verifier,signer in itertools.permutations(variants,2):
   call([builds[row['name'],verifier],'check',ns,files[signer]],log=out/'correctness.log')
  correct.append({'set':row['name'],'count_each':ns,'cross_builds':True,'byte_equal':equal,'ordered_cross_pairs':len(variants)*(len(variants)-1)})
  (out/'correctness.json').write_text(json.dumps(correct,indent=2)+'\n')
  if row.get('require_equal',False) and not equal:raise RuntimeError('Expected byte equality failed after cross verification: '+row['name'])
  print('Correctness:',row['name'],correct[-1],flush=True)
 (out/'correctness.json').write_text(json.dumps(correct,indent=2)+'\n');(out/'lengths.json').write_text(json.dumps(lengths,indent=2)+'\n')
 for row in CFG['sets']:
  for variant in row['variants']:
   v=variant['name'];key=row['name']+'-'+v;print('Timing:',key,nt,flush=True)
   text=call([builds[row['name'],v],'time',nt],log=out/(key+'-timing.log'))
   for op,value in re.findall(r'^(KeyGen|Sign|Verify|ExpandedSign)\nmedian: (\d+) cycles/ticks',text,re.M):
    size=next(x for x in lengths if x['set']==row['name'] and x['build']==v)
    results.append(dict(size,operation=op,iterations=nt,median_cycles=int(value),kcycles=int(value)/1000,source_commit=variant.get('source_commit',CFG['tested_commit']),native_instrumentation=variant.get('instrumentation',False)))
   for op in ('KeyGen','Sign','Verify'):
    if not any(x['set']==row['name'] and x['build']==v and x['operation']==op for x in results):raise RuntimeError('Missing timing: '+key+' '+op)
 (out/'timings.json').write_text(json.dumps(results,indent=2)+'\n')
 with (out/'bench.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(results[0]));w.writeheader();w.writerows(results)
 lookup={(x['set'],x['build'],x['operation']):x['kcycles'] for x in results}
 s='# '+CFG['scheme']+' — '+meta['scale']+'\n\n'
 s+='QUICK is a workflow check only. Full-scale figures are measured locally.\n\n' if quick else 'Full-scale local benchmark.\n\n'
 s+='Calibration status: '+CFG.get('calibration_status','not evaluated')+'\n\n'
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
 if CFG.get('calibration'):
  ratios=[]
  for rr in CFG['sets']:
   for vv in ('ref','avx2'):
    for op,target in zip(('KeyGen','Sign','Verify'),CFG['calibration'][rr['name']][vv]):
     ratios.append(dict(set=rr['name'],build=vv,operation=op,ratio=lookup[rr['name'],vv,op]/target))
  (out/'calibration.json').write_text(json.dumps(dict(status='REPORT_ONLY',historical_environment=CFG['calibration_environment'],measurements=ratios),indent=2)+'\n')
  s+='Historical COPA ratios are report-only; no 3% gate or machine check can reject this run. Historical messages: 59 bytes; current protocol: 32 bytes.\n\n'
 if CFG.get('source_attribution'):
  s+='\nSource attribution with the SAME new driver (thousand cycles):\n\n| Set | Build | Upstream 81389d1 | COPA vendor, instrumentation off | COPA vendor, native instrumentation on | Vendor/upstream | Instrumented/vendor |\n|---|---|---|---|---|---|---|\n'
  for rr in CFG['sets']:
   name=rr['name']
   for backend in ('ref','avx2'):
    names=(backend,'vendored_'+backend,'vendored_'+backend+'_instrumented')
    trips=['/'.join('%.3f'%lookup[name,v,op] for op in ('KeyGen','Sign','Verify')) for v in names]
    ratios=['/'.join('%.3f'%(lookup[name,names[i],op]/lookup[name,names[i-1],op]) for op in ('KeyGen','Sign','Verify')) for i in (1,2)]
    s+='| %s | %s | %s | %s | %s | %s | %s |\n'%(name,backend,*trips,*ratios)
  s+='\nAll three source conditions use the identical driver, message, RNG provider, warmup, counter and flags (only the documented native instrumentation macro differs). This controls driver differences; FULL local results are needed to attribute the user-machine discrepancy quantitatively.\n'
 if CFG.get('source_attribution'):
  decomposition=[]
  s+='\nDriver/protocol residual against historical COPA (ratios; valid for attribution only on the same local environment):\n\n| Set | Build | Vendored instrumented / upstream (code) | Vendored new-driver / historical vendor old-driver | Upstream new-driver / historical vendor old-driver |\n|---|---|---|---|---|\n'
  for rr in CFG['sets']:
   for v in ('ref','avx2'):
    code=[];driver=[];overall=[]
    for op,target in zip(('KeyGen','Sign','Verify'),CFG['calibration'][rr['name']][v]):
     u=lookup[rr['name'],v,op];w=lookup[rr['name'],'vendored_'+v+'_instrumented',op]
     code.append(w/u);driver.append(w/target);overall.append(u/target)
     decomposition.append(dict(set=rr['name'],build=v,operation=op,code_factor=w/u,driver_protocol_factor=w/target,overall_factor=u/target))
    trips=['/'.join('%.3f'%x for x in values) for values in (code,driver,overall)]
    s+='| %s | %s | %s | %s | %s |\n'%(rr['name'],v,*trips)
  reference=CFG['calibration_environment']
  comparable=(not quick and reference['cpu'] in meta['cpu'] and reference['gcc'] in meta['gcc'] and cpu==0 and ENV.get('BENCH_SEED')==reference['seed'])
  (out/'attribution.json').write_text(json.dumps(dict(local_environment_comparable=comparable,formula='overall_factor = driver_protocol_factor / code_factor',measurements=decomposition),indent=2)+'\n')
  s+='\nDecomposition: overall = driver/protocol factor divided by code factor. Historical messages were 59 bytes; the message change is part of the driver/protocol factor. Comparable FULL local environment: '+str(comparable)+'. This label is informational and never a pass/fail gate.\n'
 s+='\nFMA audit: selected integer/fixed-point paths; zero FMA instructions in all benchmark-owned objects and linked drivers. See *-fma-audit.json.\n'
 s+=CFG.get('notes','')+'\n\nBuild flags and commands per variant are in config.json, metadata.json and build logs. Algorithm sources use the pinned tested commit and, when configured, the verified vendored snapshot. Algorithm source text is unchanged. See bench/NOTES.md.\n'
 (out/'SUMMARY.md').write_text(s)
 print('PASS (historical ratios are report-only):',out/'SUMMARY.md',flush=True)
if __name__=='__main__':
 try:main()
 except Exception as e:print('FAILED:',e,file=sys.stderr,flush=True);sys.exit(1)
