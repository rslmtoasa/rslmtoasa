# VAL-19: spin-response stiffness, Stage 3a measurements

bcc Fe (`tests/spin_response/fe_bcc`), Juelich amplitudes, q along Gamma-H, direct (xi/2, -xi/2, xi/2), xi in 2pi/a. Units: Ry unless stated; omega_s and D in meV and meV A^2 (|q| in 1/A from `spin_response_dispersion.dat`). Binary: Release build of `fable_v4b` at `9fbc954` plus the `--temp` option of `run3a.py`. Runs are single-thread processes.

Definitions (from the scripts): delta_q = 1 + U Re tr chi0(q, 0) (static, eta = 0); omega_s = U M delta_q; omega_s^W = U W delta_q; D = omega_s / |q|^2; D^W = omega_s^W / |q|^2. M = Juelich moment of the run, U = U_Juelich of the run. W(0) = pi eta tr L(omega = 0) at eta = 1e-3, q = 0 (`weight_q0.py` estimator (b)), measured per mesh and temperature. D is the raw curve at each q; no D0 fit is reported.

## Per-run U, M, W(0)

| mesh N^3, T | kT (Ry) | E_F used (Ry) | U_Juelich (Ry) | M_Juelich (mu_B) | W(0) | W(0)/M - 1 |
|---|---|---|---|---|---|---|
| 60, 300 K | 1.900104e-03 | -0.06879051 | 0.07618994 | 2.067677 | 2.54474 | +23.07% |
| 160, 300 K | 1.900104e-03 | -0.06873430 | 0.07618988 | 2.067585 | 2.54592 | +23.14% |
| 160, 100 K | 6.333681e-04 | -0.06873518 | 0.07618277 | 2.067708 | 2.54720 | +23.19% |
| 160, 600 K | 3.800209e-03 | -0.06869599 | 0.07621380 | 2.067013 | 2.54453 | +23.10% |

## Zeroth moment at q = 0 (60^3, 300 K)

Integral of -Im tr chi0 / pi over -0.05..0.05 Ry, eta = 1e-3, spacing eta/4 (`moment_q0.py`):

```
m0_N60_e1e-3: raw 0.00223 (omega>0: 0.00153, omega<0: 0.00070), tail-corrected 0.00225 (factor 1.01290); M_juelich 2.06768
```

## delta_q, omega_s, D (60^3, 300 K; W = 2.54474)

```
# N xi |q|(1/A) delta_q omega_s(meV) D(meV A^2) omega_s^W(meV) D^W(meV A^2), W = 2.544741
60 0.1000 0.219600 4.878564e-03 10.45666 216.835 12.86927 266.864
60 0.2000 0.439199 1.709414e-02 36.63940 189.944 45.09301 233.768
60 0.2333 0.512399 2.221793e-02 47.62167 181.379 58.60917 223.228
60 0.2667 0.585599 2.859037e-02 61.28028 178.698 75.41915 219.928
60 0.3000 0.658799 3.405708e-02 72.99758 168.191 89.83992 206.997
60 0.3333 0.731999 3.994542e-02 85.61858 159.789 105.37290 196.656
60 0.3667 0.805199 4.629436e-02 99.22683 153.046 122.12091 188.358
60 0.4000 0.878399 5.145068e-02 110.27883 142.925 135.72287 175.902
60 0.4333 0.951599 5.692294e-02 122.00801 134.735 150.15827 165.822
60 0.4667 1.024798 6.176557e-02 132.38765 126.058 162.93275 155.143
60 0.5000 1.097998 6.618704e-02 141.86457 117.671 174.59622 144.821
```

## delta_q, omega_s, D (160^3)

T = 100 K, W = 2.54720:
```
# N xi |q|(1/A) delta_q omega_s(meV) D(meV A^2) omega_s^W(meV) D^W(meV A^2), W = 2.547204
160 0.0250 0.054900 3.606356e-04 0.77292 256.444 0.95216 315.912
160 0.0500 0.109800 1.344334e-03 2.88120 238.985 3.54935 294.405
160 0.0750 0.164700 2.897191e-03 6.20932 228.906 7.64924 281.989
160 0.1000 0.219600 4.966912e-03 10.64518 220.744 13.11377 271.935
```
T = 300 K, W = 2.54592:
```
# N xi |q|(1/A) delta_q omega_s(meV) D(meV A^2) omega_s^W(meV) D^W(meV A^2), W = 2.545922
160 0.0250 0.054900 3.410614e-04 0.73099 242.533 0.90011 298.643
160 0.0500 0.109800 1.307230e-03 2.80177 232.397 3.44997 286.162
160 0.0750 0.164700 2.847827e-03 6.10372 225.014 7.51582 277.071
160 0.1000 0.219600 4.914200e-03 10.53256 218.409 12.96928 268.938
```
T = 600 K, W = 2.54453:
```
# N xi |q|(1/A) delta_q omega_s(meV) D(meV A^2) omega_s^W(meV) D^W(meV A^2), W = 2.54453
160 0.0250 0.054900 3.185940e-04 0.68287 226.564 0.84062 278.905
160 0.0500 0.109800 1.246883e-03 2.67253 221.677 3.28994 272.888
160 0.0750 0.164700 2.740291e-03 5.87346 216.525 7.23034 266.546
160 0.1000 0.219600 4.760853e-03 10.20428 211.602 12.56165 260.485
```

## Raw D^W(q) and D(q) per temperature (160^3), meV A^2

| xi | D^W 100 K | D^W 300 K | D^W 600 K | D 100 K | D 300 K | D 600 K |
|---|---|---|---|---|---|---|
| 0.0250 | 315.912 | 298.643 | 278.905 | 256.444 | 242.533 | 226.564 |
| 0.0500 | 294.405 | 286.162 | 272.888 | 238.985 | 232.397 | 221.677 |
| 0.0750 | 281.989 | 277.071 | 266.546 | 228.906 | 225.014 | 216.525 |
| 0.1000 | 271.935 | 268.938 | 260.485 | 220.744 | 218.409 | 211.602 |

Ratio D^W(xi = 0.025) / D^W(xi = 0.1):

| T (K) | D^W(0.025) | D^W(0.1) | D^W(0.025)/D^W(0.1) |
|---|---|---|---|
| 100 | 315.912 | 271.935 | 1.1617 |
| 300 | 298.643 | 268.938 | 1.1105 |
| 600 | 278.905 | 260.485 | 1.0707 |

D^W(xi = 0.025) is 315.9 (100 K), 298.6 (300 K), 278.9 (600 K); D^W(xi = 0.1) is 271.9, 268.9, 260.5. The small-q rise of D^W changes with temperature. Not interpreted.

## W(q) (60^3, 300 K)

Runs `wq_x<xi>_e<eta>`: window [0.5, 1.6] omega_s, spacing eta/8. Columns: eta, peak (meV), peak/omega_s(M), W = pi eta h, W/W(0). `W_pos/W(0)` is by construction: it re-expresses the peak position (M peak / omega_s(M) = peak / (U delta_q)) and is not an independent weight measurement.

```
# W(0) = 2.54474 (eta 1e-3); columns: eta(Ry) peak(meV) peak/omega_s(M) W=pi*eta*h W/W(0); fit; and W_pos = M peak/omega_s(M) = peak/(U delta_q), from the peak position only
xi=0.1 eta=0.0005    12.640  1.209   2.4425  0.960  W_pos/W(0) = 0.982
xi=0.1 eta=0.00025    12.643  1.209   2.4085  0.946  W_pos/W(0) = 0.982
xi=0.1 eta=0.000125    12.667  1.211   2.3634  0.929  W_pos/W(0) = 0.984
xi=0.1 fit over 3 eta: W = 2.4711 (W/W(0) = 0.971), Gamma = 0.082 meV
xi=0.2 eta=0.002    40.124  1.095   1.9602  0.770  W_pos/W(0) = 0.890
xi=0.2 eta=0.001    38.920  1.062   1.7721  0.696  W_pos/W(0) = 0.863
xi=0.2 eta=0.0005    37.746  1.030   1.5231  0.599  W_pos/W(0) = 0.837
xi=0.2 eta=0.00025    37.527  1.024   1.2863  0.505  W_pos/W(0) = 0.832
xi=0.2 eta=0.000125    37.954  1.036   1.0437  0.410  W_pos/W(0) = 0.842
xi=0.2 fit over 5 eta: W = 2.1008 (W/W(0) = 0.826), Gamma = 2.187 meV
xi=0.2333 eta=0.002    48.681  1.022   1.4397  0.566  W_pos/W(0) = 0.831
xi=0.2333 eta=0.001    47.673  1.001   1.1174  0.439  W_pos/W(0) = 0.813
xi=0.2333 eta=0.0005    49.640  1.042   0.7685  0.302  W_pos/W(0) = 0.847
xi=0.2333 eta=0.00025    50.433  1.059   0.5152  0.202  W_pos/W(0) = 0.860
xi=0.2333 fit over 4 eta: W = 1.9666 (W/W(0) = 0.773), Gamma = 10.123 meV
xi=0.3333 eta=0.002    91.013  1.063   1.2006  0.472  W_pos/W(0) = 0.864
xi=0.3333 eta=0.001    89.844  1.049   0.9965  0.392  W_pos/W(0) = 0.853
xi=0.3333 eta=0.0005    89.335  1.043   0.8197  0.322  W_pos/W(0) = 0.848
xi=0.3333 fit over 3 eta: W = 1.4328 (W/W(0) = 0.563), Gamma = 5.436 meV
```

## Classification of interior maxima of tr L (60^3, eta = 1e-3)

A maximum counts if its height is >= 10 % of the dominant one; it is the RPA branch if within 10 % of omega_s. Last column: deviation of the nearest maximum from omega_s.

With omega_s = U M delta_q:
```
# eta=1e-3  xi omega_s(meV) | maxima: omega(meV) height/dominant class | nearest-to-omega_s deviation
0.2000    36.64 | 38.9(1.00)RPA | +0.063
0.2333    47.62 | 47.6(1.00)RPA | -0.000
0.2667    61.28 | 38.5(0.30)other  92.4(1.00)other | -0.372
0.3000    73.00 | 96.5(1.00)other | +0.322
0.3333    85.62 | 89.8(1.00)RPA  163.0(0.22)other  200.9(0.19)other  257.5(0.15)other | +0.049
0.3667    99.23 | 145.2(1.00)other  153.9(1.00)other  280.1(0.37)other  508.0(0.11)other  548.7(0.11)other  560.6(0.11)other  618.2(0.12)other | +0.464
0.4000   110.28 | 167.8(1.00)other  382.6(0.18)other  396.6(0.18)other  446.0(0.21)other  497.9(0.19)other  556.2(0.16)other  617.5(0.12)other | +0.522
0.4333   122.01 | 135.3(1.00)other  280.5(0.45)other  335.0(0.47)other  374.9(0.47)other  440.8(0.41)other | +0.109
0.4667   132.39 | 120.5(0.63)RPA  185.5(0.68)other  231.7(0.81)other  296.6(0.99)other  346.3(1.00)other  410.3(0.77)other  580.9(0.30)other | -0.090
0.5000   141.86 | 78.2(0.21)other  151.0(0.37)RPA  210.2(0.60)other  274.7(0.90)other  310.6(0.97)other  342.0(1.00)other  416.6(0.79)other  601.7(0.24)other | +0.064
# maxima: 5 RPA branch, 38 other, 10 q
```
With omega_s^W = U W delta_q, W = 2.54474:
```
# eta=1e-3  xi omega_s(meV) | maxima: omega(meV) height/dominant class | nearest-to-omega_s deviation
0.2000    45.09 | 38.9(1.00)other | -0.137
0.2333    58.61 | 47.6(1.00)other | -0.188
0.2667    75.42 | 38.5(0.30)other  92.4(1.00)other | +0.225
0.3000    89.84 | 96.5(1.00)RPA | +0.074
0.3333   105.37 | 89.8(1.00)other  163.0(0.22)other  200.9(0.19)other  257.5(0.15)other | -0.147
0.3667   122.12 | 145.2(1.00)other  153.9(1.00)other  280.1(0.37)other  508.0(0.11)other  548.7(0.11)other  560.6(0.11)other  618.2(0.12)other | +0.189
0.4000   135.72 | 167.8(1.00)other  382.6(0.18)other  396.6(0.18)other  446.0(0.21)other  497.9(0.19)other  556.2(0.16)other  617.5(0.12)other | +0.237
0.4333   150.16 | 135.3(1.00)RPA  280.5(0.45)other  335.0(0.47)other  374.9(0.47)other  440.8(0.41)other | -0.099
0.4667   162.93 | 120.5(0.63)other  185.5(0.68)other  231.7(0.81)other  296.6(0.99)other  346.3(1.00)other  410.3(0.77)other  580.9(0.30)other | +0.138
0.5000   174.60 | 78.2(0.21)other  151.0(0.37)other  210.2(0.60)other  274.7(0.90)other  310.6(0.97)other  342.0(1.00)other  416.6(0.79)other  601.7(0.24)other | -0.135
# maxima: 2 RPA branch, 41 other, 10 q
```

## Commands

Work directory: an empty directory; `S=tests/validation/spin_response_stage3a`; `rslmto.x` from a Release build at `build-0c/bin` (run3a default). Wall times are in `runs/times.txt` of that directory: 60^3 q = 0 96 s, 60^3 static (11 q) 294 s, 60^3 scan (2 q, 201 omega) 244 s, 160^3 q = 0 1167 s (300 K), 1495 s (100 K), 1512 s (600 K), 160^3 static (4 q) 2429-2465 s; at most 5 to 9 processes at once on 16 cores; 1.9 GB resident per 160^3 process.

```
python3 $S/run3a.py 60 0 --tag q0_e1e-3 --eta 1e-3 --window=-0.01,0.01,81
python3 $S/run3a.py 60 0 --tag m0_N60_e1e-3 --eta 1e-3 --window=-0.05,0.05,401
python3 $S/run3a.py 60 0.1,0.2,0.2333333333333333,0.26666666666666666,0.3,0.3333333333333333,0.36666666666666664,0.4,0.43333333333333335,0.4666666666666667,0.5 --tag st_N60_a --static
python3 $S/run3a.py 60 <xi1>,<xi2> --tag sc60_e1e-3_<i> --eta 1e-3 --window=0,0.05,201      # i = 0..4, pairs as in VAL-18
python3 $S/wq_runs.py runs
for T in 100 300 600:
  python3 $S/run3a.py 160 0 --tag q0_N160_T${T}_e1e-3 --eta 1e-3 --window=0,0.001,2 --temp $T
  python3 $S/run3a.py 160 0.025,0.05,0.075,0.1 --tag st_N160_T$T --static --temp $T      # T = 300 without --temp is the deck
python3 $S/weight_q0.py runs          # 60^3 q = 0
python3 $S/moment_q0.py runs
python3 $S/wq.py runs
python3 $S/classify.py runs [--weight W]
python3 $S/omega_s.py <dir with one st_N160 per temperature> --weight W(0)      # one directory per temperature
```
(the q = 0 160^3 runs at 100 and 600 K use `--temp`; the 300 K one is the deck value.)

## Holds by construction / not checked

- W_pos/W(0): by construction (see above).
- omega_s^W at q -> 0 equals U W delta_q; with W = W(0) this is not an independent check of D.
- 60^3 q = 0 W(0) = 2.54474 equals the value used as the script default; the 160^3 values are separate measurements.
- Not checked: other meshes (80^3, 120^3), W(q) at 160^3, convergence of D with mesh at fixed temperature, T dependence of the classification.
