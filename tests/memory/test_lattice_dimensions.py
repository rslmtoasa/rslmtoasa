"""Compile the actual lattice dimension scanner and exercise namelist ordering.

Run: python3 tests/memory/test_lattice_dimensions.py
Requires gfortran (or set FC to a gfortran-compatible compiler).
"""
import os
from pathlib import Path
import re
import subprocess
import tempfile

root = Path(__file__).resolve().parents[2]
source = (root / "source/lattice.f90").read_text()
scanner = re.search(
    r"   subroutine read_lattice_dimensions\(.*?end subroutine read_lattice_dimensions",
    source, re.S | re.I
).group()
# Only replace the logging dependency; compile the production scanner unchanged.
prefix = """
module scanner_test
implicit none
type logger
contains
procedure :: fatal
end type
type(logger) :: g_logger
contains
subroutine fatal(this, message, file, line)
class(logger) :: this
character(*) :: message, file
integer :: line
print *, message
error stop 1
end subroutine
"""
driver = """
end module
program test
use scanner_test
implicit none
integer :: ndim=20, ntype=0, nclu=0, njij=0, njijk=0, u, ios
real, allocatable :: ct(:)
character(10000) :: tag=''
namelist /lattice/ ndim, ntype, nclu, njij, njijk, ct, tag
call read_lattice_dimensions('input.nml', ndim, ntype, nclu, njij, njijk)
if (ntype /= 2 .or. ndim /= 20) error stop 2
allocate(ct(ntype))
ct = 0
open(newunit=u, file='input.nml', status='old')
read(u, nml=lattice, iostat=ios)
close(u)
if (ios /= 0) error stop 3
if (any(abs(ct - [3.0,4.0]) > 1.e-6)) error stop 4
end program
"""
cases = [
    ("ordered", "&lattice ndim=20, ntype=2, ct=3,4 /", True),
    ("array_first", "&lattice ct=3,4, ntype=2, ndim=20 /", True),
    ("comments_quotes",
     '&control ntype=999 /\n&LATTICE tag="it\'s ! ntype=999 /", '
     'ct=3,4, ! ntype=500\n NTYPE=2, ndim=20 /', True),
    ("doubled_quote", "&lattice tag='it''s ntype=999', ct=3,4, ntype=2 /", True),
    ("long_record", "&lattice tag='" + "x"*5000 + "', ct=3,4, ntype=2 /", True),
    ("repeat", "&lattice ct=3,4, ntype=1*2 /", True),
    ("null", "&lattice ct=3,4, ntype=2, ndim=, /", True),
    ("negative", "&lattice ntype=-1, ct=3,4 /", False),
    ("missing", "&control ntype=2 /", False),
]
with tempfile.TemporaryDirectory() as tmp:
    work = Path(tmp)
    (work / "scanner.F90").write_text(prefix + scanner + driver)
    subprocess.run(
        [os.environ.get("FC", "gfortran"), "-cpp", "-ffree-line-length-none",
         "-fcheck=all", "-O0", "-g", "scanner.F90", "-o", "test"],
        cwd=work, check=True,
    )
    for name, text, success in cases:
        (work / "input.nml").write_text(text + "\n")
        run = subprocess.run([str(work / "test")], cwd=work, capture_output=True, text=True)
        assert (run.returncode == 0) == success, (name, run.stdout, run.stderr)
        print("PASS", name)
