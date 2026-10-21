"""Compile and exercise the actual optional Fortran input hook when FC is available."""
import os
from pathlib import Path
import shlex
import shutil
import struct
import subprocess
import tempfile
import unittest
import numpy as np


FC=os.environ.get('FC') or shutil.which('gfortran')


@unittest.skipUnless(FC,'Set FC to run Fortran saved-Sbar hook integration test')
class HookTests(unittest.TestCase):
    def test_loader_and_rejections(self):
        source=Path(__file__).resolve().parents[2]/'source'
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            driver=p/'driver.f90'
            driver.write_text('''program driver
use precision_mod
use gf_audit_io_mod
implicit none
complex(rp) :: s(9,9,1,1)
s=0
open(16,file='view.sbar',status='old',position='append')
call gf_audit_load_sbar(s)
write(*,*) real(sum(s),rp)
end program
''')
            cmd=shlex.split(FC)+shlex.split(os.environ.get('FFLAGS',''))+[
                '-ffree-line-length-0',str(source/'precision.f90'),str(source/'gf_audit_io.f90'),str(driver),'-o',str(p/'check')]
            subprocess.run(cmd,cwd=p,check=True,capture_output=True)
            label=' VECTOR= 0.000000 0.000000 0.000000 ICLUS= 1\n'
            (p/'view.sbar').write_text(label)
            saved=p/'view.sbar.';saved.write_text(label)
            marker=struct.pack('<I',72)
            data=b''.join(marker+np.ones(9,dtype='<f8').tobytes()+marker for _ in range(9))
            binary=p/'sbar.reciprocal';binary.write_bytes(data)
            env=os.environ.copy()
            for k in ('RSLMTO_AUDIT_SBAR','RSLMTO_AUDIT_SBAR_VIEW','RSLMTO_GF_AUDIT'): env.pop(k,None)
            def run(): return subprocess.run([str(p/'check')],cwd=p,env=env,capture_output=True,text=True)
            r=run();self.assertEqual(r.returncode,0);self.assertEqual(float(r.stdout),0)
            env.update(RSLMTO_AUDIT_SBAR=str(binary),RSLMTO_AUDIT_SBAR_VIEW=str(saved),RSLMTO_GF_AUDIT='1')
            r=run();self.assertEqual(r.returncode,0,r.stderr)
            self.assertIn('DIAGNOSTIC frozen Sbar override:',r.stdout)
            self.assertEqual(float(r.stdout.splitlines()[-1]),81)
            saved.write_text(label.replace('0.000000','1.000000',1))
            self.assertNotEqual(run().returncode,0)
            saved.write_text(label)
            binary.write_bytes(data[:-1]);self.assertNotEqual(run().returncode,0)
            binary.write_bytes(data+data[:80]);self.assertNotEqual(run().returncode,0)
            binary.write_bytes(data);env.pop('RSLMTO_GF_AUDIT')
            self.assertNotEqual(run().returncode,0)

    def test_exhaustive_graph_checker(self):
        source=Path(__file__).resolve().parents[2]/'source'
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);driver=p/'graph.f90'
            driver.write_text('''program graph
use precision_mod
use gf_audit_io_mod
implicit none
integer :: nn(2,2), types(2), i
complex(rp) :: h(18,18,2,2), ho(18,18,2,2), o(18,18,2), e(18,18,2)
character :: mode
call get_command_argument(1,mode)
nn(1,:)=[2,2];nn(2,:)=[2,1];types=[1,2]
h=0;o=0;e=0
do i=1,18
o(i,i,:)=0.5_rp
end do
h(1,2,2,1)=0.1_rp;h(2,1,2,2)=0.1_rp
if (mode == 'p') h(2,1,2,2)=0.2_rp
if (mode == 'm') nn(2,2)=0
ho=0.5_rp*h
call gf_audit_graph(nn,types,h,ho,o,e,e)
end program
''')
            cmd=shlex.split(FC)+shlex.split(os.environ.get('FFLAGS',''))+[
                '-ffree-line-length-0',str(source/'precision.f90'),str(source/'gf_audit_io.f90'),str(driver),'-o',str(p/'check')]
            subprocess.run(cmd,cwd=p,check=True,capture_output=True)
            env=os.environ.copy();env['RSLMTO_GF_AUDIT_H']='1'
            def run(mode): return subprocess.run([str(p/'check'),mode],cwd=p,env=env,capture_output=True)
            self.assertEqual(run('r').returncode,0)
            r=dict(line.split() for line in (p/'gf-audit-graph.dat').read_text().splitlines())
            self.assertEqual(int(r['directed_blocks']),4);self.assertEqual(float(r['bare_h_max_Ry']),0)
            self.assertEqual(run('p').returncode,0)
            r=dict(line.split() for line in (p/'gf-audit-graph.dat').read_text().splitlines())
            self.assertAlmostEqual(float(r['bare_h_max_Ry']),0.1)
            self.assertNotEqual(run('m').returncode,0)


if __name__=='__main__': unittest.main()
