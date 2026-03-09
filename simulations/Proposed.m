%% Proposed.m
% ==========================================================================
% PROPOSED ALGORITHM: Joint Trajectory and Communication Design (JTCD)
%
% Implements the Block Successive Convex Approximation (BSCA) algorithm for
% the mAEE maximisation problem in:
%
%   "Aerial Intelligent Reflecting Surface Enabled Terahertz Covert
%    Communications in Beyond-5G Internet of Things"
%   M. Tatar Mamaghani and Y. Hong, IEEE IoT Journal.
%
% Algorithm Description:
%   Jointly optimises user scheduling (alpha), AP transmit power (p_a),
%   UCJ AN power (p_j_hat), IRS beamforming (Phi), UIRS trajectory/velocity
%   (q_r, v_r), and UCJ trajectory/velocity (q_j, v_j) via 5 sequential
%   convex sub-problems per outer iteration.
%
% Inputs (from workspace, set by SystemParams.m):
%   All system parameters, epsilon (covertness), eps_algo (convergence tol.)
%
% Outputs (to workspace):
%   ATR, APC, AEE  : performance metrics across inner iterations
%   Itr            : struct array of per-outer-iteration snapshots
%   fig1           : trajectory figure handle
% ==========================================================================

ite_index = 1;

% ── STEP 1: Feasible initialisation (circular trajectory) ──────────────
Feasible_Init
print_results(ATR, APC, AEE);

% ── STEP 2: Visualise network topology ─────────────────────────────────
Visulization

% ── STEP 3: Record initial state ───────────────────────────────────────
Out_itr = 1;
uavTrj(qr, qj, fig1_ax, lgd, Out_itr);
Itr(Out_itr).alpha    = alpha;
Itr(Out_itr).pa       = pa;
Itr(Out_itr).pj       = pj;
Itr(Out_itr).Qr       = {qr, vr};
Itr(Out_itr).Qj       = {qj, vj};
Itr(Out_itr).mAEE     = AEE(end);
Itr(Out_itr).APC      = APC(end);
Itr(Out_itr).mACT     = ATR(end);
Itr(Out_itr).AEE_Outer = AEE(end);
AEE_Outer(Out_itr)    = AEE(end);

% ── MAIN BSCA OUTER LOOP ───────────────────────────────────────────────
while (~strcmp(flag, 'converged'))

    fprintf('\n Iteration #%d starts now ...\n', Out_itr);
    Out_itr = Out_itr + 1;

    % Sub-problem 1: User Scheduling
    fprintf('\n User scheduling is being optimized...\n')
    subProb1_Alt
    print_results(ATR, APC, AEE);

    % Sub-problem 2: IRS Beamforming (SDP)
    ite_index = ite_index + 1;
    subProb2_NEW
    print_results(ATR, APC, AEE);

    % Sub-problem 3: Joint Power Allocation
    ite_index = ite_index + 1;
    fprintf('\n Joint network power allocation is being optimized... \n')
    subProb3_Alt
    print_results(ATR, APC, AEE);

    % Sub-problem 4: UIRS Trajectory and Velocity
    fprintf('\n Joint UIRS velocity and trajectory is being optimized...\n')
    subProb4
    print_results(ATR, APC, AEE);

    % Sub-problem 5: UCJ Trajectory and Velocity
    fprintf('\n Joint UCJ velocity and trajectory is being optimized...\n')
    subProb5
    print_results(ATR, APC, AEE);

    % Record iteration state
    Itr(Out_itr).alpha     = alpha;
    Itr(Out_itr).pa        = pa;
    Itr(Out_itr).pj        = pj;
    Itr(Out_itr).Qr        = {qr, vr};
    Itr(Out_itr).Qj        = {qj, vj};
    Itr(Out_itr).mAEE      = AEE(end);
    Itr(Out_itr).APC       = APC(end);
    Itr(Out_itr).mACT      = ATR(end);
    AEE_Outer(Out_itr)     = AEE(end);
    Itr(Out_itr).AEE_Outer = AEE(end);
    uavTrj(qr, qj, fig1_ax, lgd, Out_itr);

    % Convergence check: relative improvement in mAEE
    if ((AEE_Outer(Out_itr) - AEE_Outer(Out_itr-1)) / AEE_Outer(Out_itr-1) <= eps_algo)
        fprintf('\nThe overall algorithm converged in %d outer iterations.\n', Out_itr);
        flag = 'converged';
    end

end

% Final convergence plot
figure;
plot(1:ite_index, AEE, '-ob', 'LineWidth', 1.5, 'MarkerSize', 4);
xlabel('Inner Iteration Index'); ylabel('mAEE');
title('Convergence — Proposed JTCD Algorithm'); grid on;
