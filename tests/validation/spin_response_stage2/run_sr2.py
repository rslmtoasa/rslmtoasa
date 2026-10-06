"""run_sr2.py <N> <eta> <xi,xi,..> [--bin B] [--tag T] [--fermi E_F] [--uscale S] [--method M] [--env K=V] [--static]
Copy tests/spin_response/fe_bcc to ./runs/m<N>_e<eta><tag>, set the mesh to N^3, run post_processing=spin_response
at q direct (xi/2, -xi/2, xi/2), omega 0 to 0.05 Ry with spacing <= eta/4. Work directory: the current directory.
--fermi, --uscale and --env need a scratch binary (see coefficient_goldstone.patch for SR_COEF_GOLDSTONE)."""
import argparse, math, os, pathlib, re, shutil, subprocess, time
REPO = pathlib.Path(__file__).resolve().parents[3]
SP = pathlib.Path.cwd()
ap = argparse.ArgumentParser()
ap.add_argument("N", type=int); ap.add_argument("eta", type=float); ap.add_argument("xis")
ap.add_argument("--bin", default=str(REPO / "build/bin/rslmto.x")); ap.add_argument("--tag", default="")
ap.add_argument("--fermi", type=float); ap.add_argument("--uscale", type=float)
ap.add_argument("--method", default="juelich"); ap.add_argument("--env", action="append", default=[]); ap.add_argument("--static", action="store_true")
a = ap.parse_args()
xis = [float(x) for x in a.xis.split(",")]
wmax = 0.05
nw = math.ceil(wmax / (a.eta / 4) - 1e-9) + 1
d = SP / "runs" / f"m{a.N}_e{a.eta:g}{a.tag}"
shutil.rmtree(d, ignore_errors=True)
shutil.copytree(REPO / "tests/spin_response/fe_bcc", d)
text = (REPO / "tests/spin_response/fe_bcc/input.nml").read_text()
old = "pre_processing = 'bravais'"
assert old in text
for k in ("nk1", "nk2", "nk3"):
    text, n = re.subn(rf"^{k} = 12$", f"{k} = {a.N}", text, flags=re.M); assert n == 1
if a.fermi is not None:
    text, n = re.subn(r"^auto_find_fermi = \.true\.$", "auto_find_fermi = .false.", text, flags=re.M); assert n == 1
    text, n = re.subn(r"^fermi = -0\.067656$", f"fermi = {a.fermi:.10f}", text, flags=re.M); assert n == 1
text = text.replace(old, "pre_processing = 'none'\npost_processing = 'spin_response'")
q = "".join(f"q_list(:, {i+1}) = {x/2:.17g}, {-x/2:.17g}, {x/2:.17g}\n" for i, x in enumerate(xis))
win = ("omega_min = 0.0\nomega_max = 0.0\nn_omega = 1\neta = 0.0\n" if a.static else f"omega_min = 0.0\nomega_max = {wmax}\nn_omega = {nw}\neta = {a.eta:g}\n")
text += f"\n&spin_response\nmethod = '{a.method}'\nn_q = {len(xis)}\n{q}{win}/\n"
(d / "input_driver.nml").write_text(text)
env = dict(os.environ, OMP_NUM_THREADS="1")
for kv in a.env: k, v = kv.split("="); env[k] = v
if a.uscale is not None: env["SR_U_SCALE"] = repr(a.uscale)
t0 = time.time()
r = subprocess.run([a.bin, "input_driver.nml"], cwd=d, env=env, stdout=open(d / "driver.log", "w"), stderr=subprocess.STDOUT)
dt = time.time() - t0
line = f"m{a.N} eta={a.eta:g}{a.tag} nq={len(xis)} nw={nw} rc={r.returncode} wall={dt:.1f}s"
open(SP / "runs" / "times2.txt", "a").write(line + "\n"); print(line)
