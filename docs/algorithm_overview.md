# BSCA Algorithm Overview

This document describes the **Block Successive Convex Approximation (BSCA)** algorithm used to solve the mAEE optimisation problem in the paper.

---

## Optimisation Problem

The original problem is:

```
(P)  maximise   min_{k ∈ K}  (1/N) Σ_{n=1}^{N}  α_k[n] · R_k^lb[n] / P_f[n]

     subject to:
       C1:  UIRS trajectory & velocity constraints
       C2:  UCJ trajectory & velocity constraints
       C3:  UAV-UAV safety distance  ‖q_r[n] − q_j[n]‖ ≥ D_s
       C4:  Binary user scheduling  α_k[n] ∈ {0,1},  Σ_k α_k[n] ≤ 1
       C5:  IRS phase/amplitude constraints  0 ≤ ρ_l[n] ≤ 1,  0 ≤ φ_l[n] ≤ 2π
       C6:  AP and UCJ power budgets
       C7:  Covertness constraint  ζ_m[n] ≥ 1 − ε  ∀m ∈ W, n ∈ N
```

where `R_k^lb[n]` is the lower-bound covert data rate from AP to UE k in slot n,
`P_f[n]` is the total propulsion power of UIRS and UCJ in slot n,
and `ζ_m[n]` is the minimum detection error rate for Willie m.

This problem is **non-convex** due to binary variables, coupled optimisation variables, and the fractional structure of the objective.

---

## Solution Strategy

### Block Coordinate Descent (BCD)

At each outer iteration, the variables are partitioned into **5 blocks**, each solved individually while the others are held fixed:

```
Block 1: User scheduling   {α_k[n]}
Block 2: IRS beamforming   {ρ_l[n], φ_l[n]}
Block 3: Power allocation  {p_a[n], p̂_j[n]}
Block 4: UIRS trajectory   {q_r[n], v_r[n]}
Block 5: UCJ trajectory    {q_j[n], v_j[n]}
```

### Successive Convex Approximation (SCA)

Within each block, non-convex terms are replaced by **first-order Taylor approximations** around the current iterate, yielding a tractable convex surrogate solved via **CVX + MOSEK**.

---

## Algorithm Flowchart

```
┌─────────────────────────────────────────────────────────┐
│  INITIALISE: circular trajectory (Feasible_Init.m)       │
│  pa = Pa_max/2,  pj = Pj_max/2,  α = uniform           │
└──────────────────────────┬──────────────────────────────┘
                           │
                    ┌──────▼──────┐
                    │   Compute   │
                    │  Channels   │  ← Channels.m
                    │  AEE_Calc   │  ← AEE_Calc.m
                    └──────┬──────┘
                           │
               ┌───────────▼──────────────┐
               │     OUTER LOOP           │
               │  while NOT converged     │
               └───────────┬──────────────┘
                           │
          ┌────────────────▼──────────────────┐
          │  SUB-PROBLEM 1: User Scheduling    │
          │  Solve via CVX (Sch_optim.m)       │
          │  → Update α_k[n]                  │
          └────────────────┬──────────────────┘
                           │
          ┌────────────────▼──────────────────┐
          │  SUB-PROBLEM 2: IRS Beamforming    │
          │  Solve via CVX SDP                 │
          │  (IRSbeamformingNEW_optim.m)       │
          │  → Update Φ[n] (phase shifts)      │
          └────────────────┬──────────────────┘
                           │
          ┌────────────────▼──────────────────┐
          │  SUB-PROBLEM 3: Power Allocation   │
          │  Solve via CVX                     │
          │  (JointPowerWoRecast_optim.m)      │
          │  → Update p_a[n], p̂_j[n]          │
          └────────────────┬──────────────────┘
                           │
          ┌────────────────▼──────────────────┐
          │  SUB-PROBLEM 4: UIRS Trajectory    │
          │  Solve via CVX (Trj_uavIRS_optim)  │
          │  → Update q_r[n], v_r[n]           │
          └────────────────┬──────────────────┘
                           │
          ┌────────────────▼──────────────────┐
          │  SUB-PROBLEM 5: UCJ Trajectory     │
          │  Solve via CVX (Trj_CJU_optim)     │
          │  → Update q_j[n], v_j[n]           │
          └────────────────┬──────────────────┘
                           │
               ┌───────────▼──────────────┐
               │  Compute new mAEE         │
               │  Check convergence:       │
               │  |mAEE_new - mAEE_old|/  │
               │  |mAEE_old| ≤ eps_algo   │
               └───────────┬──────────────┘
                           │
                    ┌──────▼──────┐
                    │  CONVERGED? │
                    └──────┬──────┘
                        Yes│
                    ┌──────▼──────────┐
                    │  Output results  │
                    │  Save Itr struct │
                    └─────────────────┘
```

---

## Sub-Problem Details

### Sub-Problem 1: User Scheduling (`subProb1_Alt.m` → `Sch_optim.m`)

- **Challenge:** Binary constraint C4 makes this NP-hard
- **Approach:** Continuous relaxation of α_k[n] ∈ [0,1]; penalty-based method to recover near-binary solution
- **CVX formulation:** Maximise min_{k} Σ_n α·A_kn subject to Σ_k α_k[n] ≤ 1, covertness C7

### Sub-Problem 2: IRS Beamforming (`subProb2_NEW.m` → `IRSbeamformingNEW_optim.m`)

- **Challenge:** Unit-modulus constraints on phase shifts
- **Approach:** Semidefinite relaxation (SDR); reformulate as SDP with lifted variable W = φφ†
- **CVX formulation:** Hermitian PSD variable W of size L×L; constrain diagonal elements

### Sub-Problem 3: Joint Power Allocation (`subProb3_Alt.m` → `JointPowerWoRecast_optim.m`)

- **Challenge:** Rate expression is a difference of log functions (not convex in p_a, p̂_j jointly)
- **Approach:** SCA linearisation of the non-convex interference term; inner penalty loop
- **CVX formulation:** Maximise auxiliary variable η subject to linearised rate constraints and power budgets

### Sub-Problem 4 & 5: UAV Trajectory (`subProb4.m`, `subProb5.m`)

- **Challenge:** Non-convex channel gains (distance in denominator and exponent) w.r.t. trajectory
- **Approach:** First-order Taylor expansion of distance-dependent terms at current trajectory iterate
- **CVX formulation:** Quadratic and SOCP constraints for kinematics; linearised channel lower bounds

---

## Convergence and Complexity

- The BSCA algorithm is guaranteed to converge to a **stationary point** of the original problem.
- Each outer iteration increases the objective; the sequence is bounded, ensuring convergence.
- **Per-iteration complexity:** O(N³L³) for the SDP in Sub-problem 2 (dominant); O(N³) for trajectory sub-problems.
- **Typical convergence:** 20–50 outer iterations for the default parameters.

---

## Benchmark Descriptions

| Scheme | Script | Description |
|--------|--------|-------------|
| **Proposed (JTCD)** | `simulations/Proposed.m` | Full joint optimisation: trajectory + communication design |
| **Benchmark I (CD)** | `simulations/Benchmark1.m` | Communication design only; trajectory fixed to initial circular path |
| **Benchmark II (TD)** | `simulations/Benchmark2.m` | Trajectory design only; no user scheduling or power optimisation |
| **Benchmark III (IFTR)** | `experiments/Main.m` (Feasible_Init only) | Fixed initial feasible circular trajectory and equal resource allocation |

---

## Notes on `flag` Variable

The convergence flag is managed in `Proposed.m`:

```matlab
flag = 'false';
while (~strcmp(flag, 'converged'))
    % ... run sub-problems ...
    if abs(AEE_new - AEE_old) / abs(AEE_old) <= eps_algo
        flag = 'converged';
    end
end
```
