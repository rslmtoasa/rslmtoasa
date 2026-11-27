from pathlib import Path
import argparse
import csv
import os
import sys
import time
import numpy as np
root = Path(__file__).resolve().parents[3]
sys.path.insert(0,str(root/'tests/regression'))
from test_local_projection_slab import slab
parser = argparse.ArgumentParser(description='Measure the native magnetic slab at 20, 50 and 250 orders.')
parser.add_argument('driver', type=Path)
parser.add_argument('output', type=Path)
args = parser.parse_args()
exe = args.driver.resolve()
d = args.output.resolve()
d.mkdir(parents=True, exist_ok=True)
slab(d,exe,3,True,1,5,True,0,'local')
record=bytearray((d/'projection.bin').read_bytes())
env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'}
with (d/'slab_performance.csv').open('w') as report:
 w=csv.writer(report);w.writerow(['order','observable','projection','repeat','seconds','peak_rss_KiB'])
 for order in [20,50,250]:
  for repeat in [0,1]:
   for name,proj in [('spin_torque','seed'),('spin_torque_onsite','seed'),('spin_torque_onsite','local')]:
    record[12:16]=np.int32(order).tobytes();record[-56:-26]=name.encode().ljust(30);record[-26:-10]=proj.encode().ljust(16)
    (d/'projection.bin').write_bytes(record)
    start=time.monotonic()
    pid=os.fork()
    if pid==0:
     os.chdir(d)
     with open('log','w') as log:
      os.dup2(log.fileno(),1);os.dup2(log.fileno(),2);os.execve(exe,[str(exe),'projection'],env)
    _,status,usage=os.wait4(pid,0)
    assert status==0
    w.writerow([order,name,proj,repeat,time.monotonic()-start,usage.ru_maxrss]);report.flush()
