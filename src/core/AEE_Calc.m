%% AEE_Calc.m
% ==========================================================================
% OBJECTIVE FUNCTION CALCULATOR
% Computes the three key performance metrics at the current iterate:
%
%   ATR(ite_index)  = minimum Average Covert Throughput (mACT) [normalised]
%   APC(ite_index)  = Average Propulsion Power Consumption [normalised]
%   AEE(ite_index)  = minimum Average Energy Efficiency (mAEE) [normalised]
%
% The mAEE is defined as:
%   mAEE = min_{k} (1/N) Σ_n α_k[n] · R_k^lb[n] / P_f[n]
%
% where R_k^lb[n] is the lower-bound covert rate (bits/s) for UE k in
% slot n (Eq. 20 in paper), and P_f[n] = P_f_r[n] + P_f_j[n] is the
% total UAV propulsion power.
%
% Normalisation: outputs are normalised by (Po+Pi)*2 so that propulsion
% power is expressed relative to the hover baseline.
%
% Dependencies:
%   Channels.m          (computes G_auk, G_jk from current trajectories)
%   flight_pow.m        (computes PF_r, PF_j from current velocities)
%
% Inputs (from workspace):
%   alpha, pa, pj, qr, vr, qj, vj   : current iterate variables
%   BW, numUsr, N, Po, Pi            : system parameters
%
% Outputs (to workspace):
%   ATR(ite_index), APC(ite_index), AEE(ite_index) appended to arrays
% ==========================================================================

% Recompute channel gains at current trajectory iterate
Channels

% Instantaneous covert throughput for each UE in each slot [bits/s, normalised]
ITR = alpha * BW .* log2(1 + repmat(pa, numUsr, 1) .* G_auk ./ ...
      (0.5 * repmat(pj, numUsr, 1) .* G_jk + 1));

% mACT: minimum (over UEs) of average (over time slots) covert throughput
ATR(ite_index) = min(mean(ITR, 2));

% UAV propulsion power [W] for UIRS and UCJ at each time slot
[PF_r, PF_j] = flight_pow(vr, vj);

% Normalised average flight power (relative to 2×hover baseline)
IFP = (PF_r + PF_j) / (Po + Pi) / 2;
APC(ite_index) = mean(IFP);

% mAEE: minimum (over UEs) of average (over slots) rate-to-power ratio
AEE(ite_index) = min(mean(ITR ./ repmat(IFP, numUsr, 1), 2));
