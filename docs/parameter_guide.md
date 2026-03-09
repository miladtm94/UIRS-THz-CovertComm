# System Parameter Guide

This document explains every parameter defined in `src/core/SystemParams.m` and `src/core/Flightconstants.m`.

---

## Communication Parameters (`SystemParams.m`)

| Parameter | Symbol | Value | Unit | Description |
|-----------|--------|-------|------|-------------|
| `fc` | f_c | 300e9 | Hz | THz carrier frequency (0.3 THz) |
| `C` | c | 3e8 | m/s | Speed of light |
| `lambda_c` | λ_c | C/fc | m | Carrier wavelength |
| `PSD` | — | -250 dBm/Hz | W/Hz | Power spectral density of noise |
| `BW` | W | 10e9 → 100 | Hz → normalized | Transmission bandwidth (normalised) |
| `delta2` | σ² | PSD·BW | W | AWGN noise power at UEs |
| `hj0` | h_j0 | — | — | UCJ channel reference gain at 1 m |
| `hr0` | h_r0 | — | — | UIRS cascaded channel reference gain |
| `kf` | κ(f_c,μ) | 3.2094e-4 | 1/m | Molecular absorption coefficient |
| `Pa_max` | p_a^max | 1 | W | Maximum AP transmit power |
| `Pj_max` | p̂_j^max | 1 | W | Maximum UCJ AN power |
| `P_tot` | P_tot | 40 | W | Total network power budget |
| `pl` | ρ | 2.3 | — | Path-loss exponent |
| `epsilon` | ε | 0.01 | — | Covertness constraint (detection error rate) |

---

## IRS Parameters (`SystemParams.m`)

| Parameter | Symbol | Value | Unit | Description |
|-----------|--------|-------|------|-------------|
| `Lx` | L_x | 6 | — | Number of IRS elements along x-axis |
| `Ly` | L_y | 5 | — | Number of IRS elements along y-axis |
| `L` | L | Lx·Ly = 30 | — | Total number of IRS reflecting elements |
| `delta_x` | δ_x | 1e-3 | m | IRS element spacing along x |
| `delta_y` | δ_y | 1e-3 | m | IRS element spacing along y |

---

## Geometric / Network Topology Parameters

| Parameter | Symbol | Value | Unit | Description |
|-----------|--------|-------|------|-------------|
| `H` | H | 50 | m | UAV reference altitude |
| `Hr` | H_r | H = 50 | m | UIRS flight altitude |
| `Hj` | H_j | H = 50 | m | UCJ flight altitude |
| `R` | R_o | 2H = 100 | m | UAV circular flight radius |
| `R1` | R_1 | 2H-1 = 99 | m | Minimum user distribution radius |
| `R2` | R_2 | 4H = 200 | m | Maximum user distribution radius |
| `Ds` | D_s | 0.1·H = 5 | m | Minimum UAV-UAV safety distance |
| `qa` | q_a | [0;0;0] | m | AP location (at origin) |
| `qr_I` | q_r^I | [Rr;0;Hr] | m | UIRS initial position |
| `qj_I` | q_j^I | [Rj;0;Hj] | m | UCJ initial position |
| `numUsr` | K | 5 | — | Number of ground UEs |

---

## Mission / Discretisation Parameters

| Parameter | Symbol | Value | Unit | Description |
|-----------|--------|-------|------|-------------|
| `T` | T | 30 | s | Total mission duration |
| `N` | N | 30 | — | Number of time slots |
| `dt` | δ_t | T/N = 1 | s | Time slot duration |
| `Amax` | a^max | 6 | m/s² | Maximum UAV acceleration |
| `Vmax` | v^max | 25 | m/s | Maximum UAV speed |

---

## UAV Flight Power Constants (`Flightconstants.m`)

These constants correspond to a rotary-wing UAV with typical parameters from the literature [Zeng et al., 2019].

| Parameter | Symbol | Value | Unit | Description |
|-----------|--------|-------|------|-------------|
| `omega` | Ω | 300 | rad/s | Rotor angular velocity |
| `r` | r | 0.4 | m | Rotor blade radius |
| `rho` | ρ_air | 1.225 | kg/m³ | Air density |
| `s` | s | 0.05 | — | Rotor solidity |
| `A` | A | 0.503 | m² | Rotor disc area |
| `v0` | v_0 | 4.03 | m/s | Mean rotor induced velocity in hover |
| `d0` | d_0 | 0.6 | — | Fuselage drag ratio |
| `delta_coeff` | δ | 0.012 | — | Profile drag coefficient |
| `k_corr` | k | 0.1 | — | Incremental correction factor |
| `Wgt` | W | 20 | N | UAV weight |
| `Po` | P_o | ~79.86 | W | Blade profile power (hover) |
| `Pi` | P_i | ~88.63 | W | Induced power (hover) |
| `c0` | c_0 | 3/(Ω²r²) | s²/m² | Speed-dependent profile power coefficient |
| `c1` | c_1 | 0.5·d_0·ρ·s·A | — | Parasitic drag power coefficient |
| `c2` | c_2 | 1/(2v_0²) | s²/m² | Induced power coefficient |

### Flight Power Formula

The UAV propulsion power at speed v (in m/s) is:

```
P_f(v) = P_o(1 + c_0 ‖v‖²) + c_1 ‖v‖³ + P_i · √(√(1 + c_2²‖v‖⁴) − c_2‖v‖²)
```

The hover power is: `P̄_f = P_o + P_i ≈ 168.5 W`

---

## Algorithm Convergence Parameters (`SystemParams.m`)

| Parameter | Value | Description |
|-----------|-------|-------------|
| `eps_algo` | 1e-1 | Outer loop fractional convergence tolerance |
| `eps_frac_algo` | 1e-2 | Dinkelbach's algorithm convergence tolerance |
| `flag` | 'false' | Convergence flag (set to 'converged' when done) |

---

## Scaling Factors (for Figures)

When plotting results, the raw outputs are scaled for readability:

| Quantity | Raw Unit | Scaled Unit | Scale Factor |
|----------|----------|-------------|-------------|
| mAEE | bits/J | kbits/J | `1e8 / ((Po+Pi)*2) / 1e3` |
| mACT | — | Mbps | `1e8 / 1e6` |
| APC | — | kW | `(Po+Pi)*2 / 1e3` |

---

## Notes on `SystemParams_fc.m`

This is a variant of `SystemParams.m` used for the **carrier frequency sweep** experiment.
It accepts externally set `fc` and `kf` variables before the script is called, allowing the caller to override the default 300 GHz and 3.2094e-4 values for `κ(f_c, μ)`.

Usage:
```matlab
fc = 350e9;           % Set carrier frequency
kf = ka(fc, T, p);   % Compute absorption coefficient
SystemParams_fc       % Load all other parameters using these fc/kf values
```
