%% Benchmark2.m
% ==========================================================================
% BENCHMARK II: Trajectory Design Only (TD) — No User Scheduling or Power Opt.
%
% This benchmark scheme optimises only the UIRS and UCJ trajectories while
% fixing user scheduling to a round-robin policy and using equal power
% allocation. IRS beamforming is still optimised.
%
% Variables optimised:
%   - UIRS trajectory and velocity (q_r, v_r)
%   - UCJ trajectory and velocity (q_j, v_j)
%   - IRS beamforming (Phi)
%
% Fixed (not optimised):
%   - User scheduling (alpha): round-robin
%   - AP power (p_a): uniform = Pa_max / 2
%   - UCJ AN power (p_j): uniform = Pj_max / 2
%
% Reference: Section IV-B, Benchmark II in the paper.
%
% Usage: Run SystemParams.m in the workspace before calling this script.
% ==========================================================================

ite_index = 1;

Feasible_Init
print_results(ATR, APC, AEE);
Visulization

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

while (~strcmp(flag, 'converged'))

    fprintf('\n Iteration #%d starts now ...\n', Out_itr);
    Out_itr = Out_itr + 1;

    % NOTE: User scheduling and power are NOT updated (Benchmark II — TD only)

    % Sub-problem 2: IRS Beamforming (still optimised)
    ite_index = ite_index + 1;
    subProb2_NEW
    print_results(ATR, APC, AEE);

    % Sub-problem 4: UIRS Trajectory and Velocity
    fprintf('\n Joint UIRS velocity and trajectory is being optimized...\n')
    subProb4
    print_results(ATR, APC, AEE);

    % Sub-problem 5: UCJ Trajectory and Velocity
    fprintf('\n Joint UCJ velocity and trajectory is being optimized...\n')
    subProb5
    print_results(ATR, APC, AEE);

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

    if ((AEE_Outer(Out_itr) - AEE_Outer(Out_itr-1)) / AEE_Outer(Out_itr-1) <= eps_algo)
        fprintf('\nBenchmark II (TD) converged in %d outer iterations.\n', Out_itr);
        flag = 'converged';
    end

end

figure;
plot(1:ite_index, AEE, '-^r', 'LineWidth', 1.5, 'MarkerSize', 4);
xlabel('Inner Iteration Index'); ylabel('mAEE');
title('Convergence — Benchmark II (TD, Fixed Scheduling/Power)'); grid on;
