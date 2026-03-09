%% quick_start.m
% ==========================================================================
% QUICK START EXAMPLE
% Aerial IRS-Enabled THz Covert Communications — UIRS-THz-CovertComm
%
% This script demonstrates how to run the proposed BSCA algorithm with
% default parameters and display key performance metrics.
%
% Prerequisites:
%   - MATLAB R2021a or later
%   - CVX 2.2+ installed and configured (http://cvxr.com/cvx/)
%   - MOSEK solver (recommended) or SDPT3/SeDuMi
%
% Usage:
%   Run this script from the repository root:
%       >> cd /path/to/UIRS-THz-CovertComm
%       >> quick_start
% ==========================================================================

clc; clear all; close all;

%% ── 1. Set up paths ────────────────────────────────────────────────────
fprintf('Setting up paths...\n');
addpath(genpath('src'));
addpath('simulations');

%% ── 2. Set simulation parameters ───────────────────────────────────────
T       = 30;    % Mission duration [s]
epsilon = 0.01;  % Covertness constraint (max allowed detection error rate)
N       = 30;    % Number of time slots

fprintf('Parameters: T=%ds, epsilon=%.3f, N=%d\n', T, epsilon, N);

%% ── 3. Load system parameters ───────────────────────────────────────────
fprintf('\nLoading system parameters...\n');
SystemParams  % Defined in src/core/SystemParams.m

%% ── 4. Run the proposed JTCD algorithm ─────────────────────────────────
fprintf('\n========================================\n');
fprintf('  Running Proposed Algorithm (JTCD)...\n');
fprintf('========================================\n');
Proposed      % Defined in simulations/Proposed.m

%% ── 5. Display final results ────────────────────────────────────────────
fprintf('\n========================================\n');
fprintf('  Final Results Summary\n');
fprintf('========================================\n');

Flightconstants  % Load Po, Pi for scaling

scale_AEE = 1e8 / ((Po + Pi) * 2);  % Normalisation factor

fprintf('  Minimum Average Energy Efficiency (mAEE): %.4f kbits/J\n', ...
        scale_AEE * AEE(end) / 1e3);
fprintf('  Min. Average Covert Throughput (mACT):   %.4f Mbps\n', ...
        1e8 * ATR(end) / 1e6);
fprintf('  Average Power Consumption (APC):         %.4f kW\n', ...
        (Po + Pi) * 2 * APC(end) / 1e3);
fprintf('  Total outer iterations to convergence:   %d\n', Out_itr);
fprintf('========================================\n\n');

%% ── 6. Quick convergence plot ───────────────────────────────────────────
figure('Name', 'Convergence - mAEE vs Iteration');
plot(0:length(AEE)-1, scale_AEE * AEE / 1e3, '-ob', ...
     'LineWidth', 2, 'MarkerSize', 5);
xlabel('Inner Iteration Index');
ylabel('mAEE [kbits/Joule]');
title('Convergence of Proposed JTCD Algorithm');
grid on;
fprintf('Convergence plot displayed.\n');
fprintf('To reproduce all paper figures, run: experiments/Main.m\n\n');
