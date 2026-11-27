"""post_processing = 'bdg' on the bcc-Fe density_of_states deck (6x6x6 full-BZ mesh, 300 K, E_F solved for 8 electrons).

Usage: route.py <rslmto.x> <deck directory> <scratch directory> <mpiexec>
The deck directory is only read. Checks:
  (a) post_processing = 'density_of_states' with auto_find_fermi solves E_F; that E_F goes into &energy fermi of the BdG
      deck. N_normal of the BdG route (Delta = 0 pass) against the N printed by the density_of_states route:
      tolerance 1e-7 electrons, the printed precision (EF to 1e-10 Ry, N to 9 digits) plus eigensolver round-off of two
      independent codes.
  (b) consistency: Delta after 15 iterations at lambda_dd = 2 eV from 1 and 2 MPI ranks, max|diff|/max|Delta| <= 1e-12.
  (c) the run with reciprocal_mode = 'generalized_overlap_proxy' exits non-zero with the guard message.
  (d) lambda = 0: bdg_dos.dat total DOS (electron-projected, |u|^2 = 1 on electron states and 0 on hole states at
      Delta = 0) against dos_kspace.dat of the density_of_states run (gaussian, same sigma, same mesh, same E_F). Both grids
      are 1001 points on E_F +- 0.5 Ry (the DOS deck range is absolute, so it is centred on E_F); dos_kspace.dat prints 6
      decimals (states/Ry), so the tolerance is 5e-6 per point (rounding 5e-7 plus 1e-10 Ry of E_F truncation times the
      DOS slope of about 1e4 states/Ry^2; here mu_dos - E_F is below 3e-8 Ry).
  (e) restart: a max_iter = 0 run with delta_file = the bdg_delta.dat of the lambda = 2 eV run. Delta of the new file
      equals the old to 1e-14 relative (I/O round trip, plumbing) and bdg_dos.dat has n_dos rows.
"""
import pathlib
import re
import shutil
import subprocess
import sys

binary, deck, scratch, mpiexec = sys.argv[1:5]
RY2EV, N_DOS, SIGMA, WINDOW = 13.605703976, 1001, 0.005, 0.5  # sigma and window in Ry
DOS_KEYS = f"dos_window = {WINDOW * RY2EV!r}\ndos_sigma = {SIGMA * RY2EV!r}\nn_dos = {N_DOS}\n"
BDG = "\n&bdg\nlambda(3,3,1) = {lam}\nmix = 0.5\ntol = 1.0e-12\nmax_iter = {iters}\n" + DOS_KEYS + "/\n"
TOL_N, TOL_RANK, TOL_DOS, TOL_RESTART = 1.0e-7, 1.0e-12, 5.0e-6, 1.0e-14

shutil.rmtree(scratch, ignore_errors=True)
shutil.copytree(deck, scratch)
text = pathlib.Path(deck, "input.nml").read_text()
text = text.replace("post_processing = 'density_of_states'", "post_processing = 'bdg'")
text = re.sub(r"use_symmetry_reduction = \.true\.", "use_symmetry_reduction = .false.", text)
text = re.sub(r"^(\s*nk[123]) = 12$", r"\1 = 6", text, flags=re.M)
for key in ("post_processing = 'bdg'", "use_symmetry_reduction = .false.", "nk3 = 6"):
    if key not in text:
        sys.exit("route.py: " + key + " not set in " + deck + "/input.nml")


def run(name, deck_text, ranks=1):
    pathlib.Path(scratch, name + ".nml").write_text(deck_text)
    command = [binary, name + ".nml"]
    if ranks > 1:
        command = [mpiexec, "-np", str(ranks)] + command
    result = subprocess.run(command, cwd=scratch, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    pathlib.Path(scratch, name + ".log").write_text(result.stdout)
    return result


def delta_file():
    rows = [line.split() for line in pathlib.Path(scratch, "bdg_delta.dat").read_text().splitlines()[1:]]
    return [float(x) for row in rows for x in row[2:4]]


failures = []

# (a) lambda = 0 against the normal-state k-space occupation sum
dos_text = text.replace("post_processing = 'bdg'", "post_processing = 'density_of_states'")
dos_text = dos_text.replace("auto_find_fermi = .false.", "auto_find_fermi = .true.")
dos = run("dos", dos_text)
n_ref = re.search(r"Canonical k-space occupations: EF=\s*(\S+) Ry, N=\s*(\S+),", dos.stdout)
if not n_ref:
    sys.exit("route.py: (a) density_of_states run failed, see " + scratch + "/dos.log")
text = re.sub(r"^fermi = .*$", "fermi = " + n_ref.group(1).replace("E", "d"), text, flags=re.M)
zero = run("lambda0", text + BDG.format(lam=0.0, iters=0))
n_bdg = re.search(r"bdg: N_normal =\s*(\S+) electrons", zero.stdout)
if dos.returncode or zero.returncode or not n_bdg:
    sys.exit("route.py: (a) runs failed, see " + scratch + "/dos.log and lambda0.log")
ef, n_ref, n_bdg = n_ref.group(1), float(n_ref.group(2)), float(n_bdg.group(1))
# same run on the BdG grid: the energy range of the DOS deck is absolute and is read in single precision, so centre it on
# mu_dos, E_F rounded to a multiple of 2**-24 (mu_dos +- 0.5 are then exact in single precision)
mu_dos = round(float(ef) * 2**24) / 2**24
for key, shift in (("dos_energy_min", -WINDOW), ("dos_energy_max", WINDOW)):
    dos_text = re.sub(r"^(\s*" + key + r") = .*$", r"\1 = " + repr(mu_dos + shift), dos_text, flags=re.M)
dos_text = re.sub(r"^(\s*n_energy_points) = .*$", r"\1 = " + str(N_DOS), dos_text, flags=re.M)
if run("dos_grid", dos_text).returncode:
    sys.exit("route.py: (d) density_of_states run on the BdG grid failed, see " + scratch + "/dos_grid.log")
print(f"(a) E_F {ef} Ry  N_normal bdg {n_bdg:.10f}  density_of_states {n_ref:.10f}  |diff| {abs(n_bdg - n_ref):.2e}  tol {TOL_N:.0e}")
if abs(n_bdg - n_ref) > TOL_N:
    failures.append("a")

# (d) lambda = 0 DOS of the BdG route (mu = mu_dos) against the normal-state DOS
if run("lambda0_dos", re.sub(r"^fermi = .*$", "fermi = " + repr(mu_dos).replace("e", "d"), text, flags=re.M)
       + BDG.format(lam=0.0, iters=0)).returncode:
    sys.exit("route.py: (d) lambda = 0 run at mu_dos failed, see " + scratch + "/lambda0_dos.log")
ref = [[float(x) for x in line.split()] for line in pathlib.Path(scratch, "dos_kspace.dat").read_text().splitlines()
       if not line.startswith("#")]
mine = [[float(x) for x in line.split()] for line in pathlib.Path(scratch, "bdg_dos.dat").read_text().splitlines()[1:]]
if len(ref) != N_DOS or len(mine) != N_DOS:
    sys.exit("route.py: (d) grid sizes " + str(len(ref)) + " and " + str(len(mine)) + ", expected " + str(N_DOS))
grid = max(abs(r[0] - m[0] / RY2EV - (mu_dos - float(ef))) for r, m in zip(ref, mine))
diff = max(abs(r[1] - m[1] * RY2EV) for r, m in zip(ref, mine))
print(f"(d) lambda = 0 DOS: max|E grid diff| {grid:.2e} Ry (6 printed decimals)  max|DOS diff| {diff:.2e} states/Ry (peak {max(r[1] for r in ref):.3f})  tol {TOL_DOS:.0e}")
if grid > 1.0e-6 or diff > TOL_DOS:
    failures.append("d")

# (b) 1 rank against 2 ranks
nonzero = text + BDG.format(lam=2.0, iters=15)
if run("rank1", nonzero).returncode:
    sys.exit("route.py: (b) 1-rank run failed, see " + scratch + "/rank1.log")
d1 = delta_file()
if run("rank2", nonzero, ranks=2).returncode:
    sys.exit("route.py: (b) 2-rank run failed, see " + scratch + "/rank2.log")
d2 = delta_file()
scale = max(abs(x) for x in d1)
diff = max(abs(x - y) for x, y in zip(d1, d2)) / scale
print(f"(b) max|Delta| {scale:.6e} Ry  max|diff|/max|Delta| {diff:.2e}  tol {TOL_RANK:.0e}")
if scale < 1.0e-5 or diff > TOL_RANK:
    failures.append("b")

# (e) restart from the bdg_delta.dat of the 2-rank run
first = pathlib.Path(scratch, "bdg_delta.dat").read_text()
pathlib.Path(scratch, "restart_in.dat").write_text(first)
restart = run("restart", text + BDG.format(lam=2.0, iters=0).replace("max_iter = 0", "max_iter = 0\ndelta_file = 'restart_in.dat'"))
if restart.returncode:
    sys.exit("route.py: (e) restart run failed, see " + scratch + "/restart.log")
d3 = delta_file()
n_rows = len(pathlib.Path(scratch, "bdg_dos.dat").read_text().splitlines()) - 1
rt = max(abs(x - y) for x, y in zip(d2, d3)) / scale
print(f"(e) restart: max|Delta diff|/max|Delta| {rt:.2e}  tol {TOL_RESTART:.0e}  bdg_dos.dat rows {n_rows}  expected {N_DOS}")
if rt > TOL_RESTART or n_rows != N_DOS:
    failures.append("e")

# (c) generalized S(k) must stop before any BdG work
guard = run("guard", nonzero.replace("&reciprocal", "&reciprocal\n   reciprocal_mode = 'generalized_overlap_proxy'", 1))
ok = guard.returncode != 0 and "reciprocal_mode must be 'ham_only'" in guard.stdout
print(f"(c) reciprocal_mode = generalized_overlap_proxy: exit code {guard.returncode}, guard message {'found' if ok else 'missing'}")
if not ok:
    failures.append("c")

sys.exit("route.py: failed " + ",".join(failures) if failures else 0)
