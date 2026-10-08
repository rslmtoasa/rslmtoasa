"""run3a.py <N> <xi,xi,..> --tag T [--static | --window wmin,wmax,nw --eta E] [--temp K] [--bin B] [--method M]
Copy tests/spin_response/fe_bcc to ./runs/<tag>, set the mesh to N^3, run post_processing=spin_response at q direct
(xi/2, -xi/2, xi/2) (Gamma-H, xi in 2pi/a). --temp K replaces the deck temperature (300 K). --static: n_omega = 1, omega = 0, eta = 0. Work directory: the current directory."""
import argparse, os, pathlib, re, shutil, subprocess, time
REPO = pathlib.Path(__file__).resolve().parents[3]
ap = argparse.ArgumentParser()
ap.add_argument("N", type=int); ap.add_argument("xis"); ap.add_argument("--tag", required=True)
ap.add_argument("--static", action="store_true"); ap.add_argument("--window"); ap.add_argument("--eta", type=float, default=0.0)
ap.add_argument("--bin", default=str(REPO / "build-0c/bin/rslmto.x")); ap.add_argument("--method", default="juelich"); ap.add_argument("--temp", type=float)
a = ap.parse_args()
xis = [float(x) for x in a.xis.split(",")]
d = pathlib.Path.cwd() / "runs" / a.tag
shutil.rmtree(d, ignore_errors=True)
shutil.copytree(REPO / "tests/spin_response/fe_bcc", d)
text = (REPO / "tests/spin_response/fe_bcc/input.nml").read_text()
old = "pre_processing = 'bravais'"
assert old in text
for k in ("nk1", "nk2", "nk3"):
    text, n = re.subn(rf"^{k} = 12$", f"{k} = {a.N}", text, flags=re.M); assert n == 1
if a.temp is not None:
    text, n = re.subn(r"^temperature = 300.0", f"temperature = {a.temp:g}", text, flags=re.M); assert n == 1
text = text.replace(old, "pre_processing = 'none'\npost_processing = 'spin_response'")
q = "".join(f"q_list(:, {i+1}) = {x/2:.17g}, {-x/2:.17g}, {x/2:.17g}\n" for i, x in enumerate(xis))
if a.static:
    win = "omega_min = 0.0\nn_omega = 1\neta = 0.0\n"
else:
    wmin, wmax, nw = a.window.split(","); win = f"omega_min = {wmin}\nomega_max = {wmax}\nn_omega = {nw}\neta = {a.eta:g}\n"
text += f"\n&spin_response\nmethod = '{a.method}'\nn_q = {len(xis)}\n{q}{win}/\n"
(d / "input_driver.nml").write_text(text)
t0 = time.time()
r = subprocess.run([a.bin, "input_driver.nml"], cwd=d, env=dict(os.environ, OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1"),
                   stdout=open(d / "driver.log", "w"), stderr=subprocess.STDOUT)
for f in d.iterdir():  # keep only the response outputs, the log and the input; the SCF-database copies are not needed afterwards
    if not (f.name.startswith("spin_response_") or f.name in ("driver.log", "input_driver.nml")): f.unlink() if f.is_file() else shutil.rmtree(f)
line = f"{a.tag} N={a.N} nq={len(xis)} rc={r.returncode} wall={time.time()-t0:.1f}s"
open(d.parent / "times.txt", "a").write(line + "\n"); print(line)
