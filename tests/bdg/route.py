"""post_processing = 'bdg' on the bcc-Fe density_of_states deck (6x6x6 full-BZ mesh, 300 K, E_F solved for 8 electrons).

Usage: route.py <rslmto.x> <deck directory> <scratch directory> <mpiexec>
The deck directory is only read. Checks:
  (a) post_processing = 'density_of_states' with auto_find_fermi solves E_F; that E_F goes into &energy fermi of the BdG
      deck. N_normal of the BdG route (Delta = 0 pass) against the N printed by the density_of_states route:
      tolerance 1e-7 electrons, the printed precision (EF to 1e-10 Ry, N to 9 digits) plus eigensolver round-off of two
      independent codes.
  (b) consistency: Delta after 15 iterations at lambda_dd = 2 eV from 1 and 2 MPI ranks, max|diff|/max|Delta| <= 1e-12.
  (c) the run with reciprocal_mode = 'generalized_overlap_proxy' exits non-zero with the guard message.
"""
import pathlib
import re
import shutil
import subprocess
import sys

binary, deck, scratch, mpiexec = sys.argv[1:5]
BDG = "\n&bdg\nlambda(3,3,1) = {lam}\nmix = 0.5\ntol = 1.0e-12\nmax_iter = {iters}\n/\n"
TOL_N, TOL_RANK = 1.0e-7, 1.0e-12

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
if "auto_find_fermi = .true." not in dos_text:
    sys.exit("route.py: auto_find_fermi not set in " + deck + "/input.nml")
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
print(f"(a) E_F {ef} Ry  N_normal bdg {n_bdg:.10f}  density_of_states {n_ref:.10f}  |diff| {abs(n_bdg - n_ref):.2e}  tol {TOL_N:.0e}")
if abs(n_bdg - n_ref) > TOL_N:
    failures.append("a")

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

# (c) generalized S(k) must stop before any BdG work
guard = run("guard", nonzero.replace("&reciprocal", "&reciprocal\n   reciprocal_mode = 'generalized_overlap_proxy'", 1))
ok = guard.returncode != 0 and "reciprocal_mode must be 'ham_only'" in guard.stdout
print(f"(c) reciprocal_mode = generalized_overlap_proxy: exit code {guard.returncode}, guard message {'found' if ok else 'missing'}")
if not ok:
    failures.append("c")

sys.exit("route.py: failed " + ",".join(failures) if failures else 0)
