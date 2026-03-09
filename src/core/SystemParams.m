%% SystemParams.m
% ==========================================================================
% SYSTEM PARAMETERS — UIRS-UCJ THz Covert Communication System
%
% Defines all system-level parameters for the simulation. This script must
% be run (not called as a function) before any simulation script.
%
% Categories:
%   1. Communication constants (THz channel, bandwidth, noise)
%   2. IRS array parameters
%   3. Network geometry and topology
%   4. UAV mission parameters (altitude, speed, duration)
%   5. Algorithm convergence settings
%
% Reference:
%   M. Tatar Mamaghani and Y. Hong, "Aerial Intelligent Reflecting Surface
%   Enabled Terahertz Covert Communications in Beyond-5G IoT," IEEE IoT J.
%
% See also: docs/parameter_guide.md for full parameter descriptions.
% ==========================================================================

%% ── Algorithm flags & convergence thresholds ───────────────────────────
flag           = 'false';     % Outer-loop convergence flag
eps_algo       = 1e-1;        % Outer-loop fractional convergence tolerance
eps_frac_algo  = 1e-2;        % Dinkelbach's algorithm termination criterion

%% ── THz Communication Constants ────────────────────────────────────────
fc       = 300e9;             % Carrier frequency: 0.3 THz [Hz]
C        = 3e8;               % Speed of light [m/s]
lambda_c = C / fc;            % Carrier wavelength [m]

PSD    = db2pow(-250);        % Noise PSD: -250 dBm/Hz → W/Hz
BW     = 10e9;                % Allocated bandwidth [Hz]
delta2 = PSD * BW;            % Noise power at UEs [W]
BW     = 1e2;                 % Normalised BW used in rate computations

% Reference channel gains (normalised by noise power)
hj0 = (lambda_c / 4 / pi) / sqrt(delta2);            % UCJ reference gain at 1 m
hr0 = (lambda_c / (8 * (sqrt(pi^3)))) / sqrt(delta2); % UIRS cascaded reference gain

% THz molecular absorption coefficient κ(fc, μ) at 0.3 THz
% (standard atmosphere: 296 K, 101325 Pa, ~50% humidity)
kf = 3.2094e-04;              % [1/m]

%% ── Power Constraints ───────────────────────────────────────────────────
Pa_max = 1;    % Maximum AP transmit power [W]
P_tot  = 40;   % Total network power budget [W]
Pj_max = 1;    % Maximum UCJ AN power [W]
pl     = 2.3;  % Path-loss exponent (2 ≤ ρ ≤ 4)

%% ── IRS Array Parameters ────────────────────────────────────────────────
delta_x = 1e-3;   % IRS element spacing along x-axis [m]
delta_y = 1e-3;   % IRS element spacing along y-axis [m]
Lx      = 6;      % Number of elements along x-axis
Ly      = 5;      % Number of elements along y-axis
L       = Lx * Ly; % Total IRS elements (30)

%% ── Network Geometry ────────────────────────────────────────────────────
H   = 50;          % Reference UAV altitude [m]
Ds  = 0.1 * H;    % Minimum UAV-UAV safety distance [m]
Hr  = H;           % UIRS flight altitude [m]
Hj  = H;           % UCJ flight altitude [m]
R   = 2 * H;       % Circular flight radius [m]
R1  = 2 * H - 1;  % Inner boundary of UE distribution [m]
R2  = 4 * H;       % Outer boundary of UE distribution [m]
Rr  = R;           % UIRS circular flight radius [m]
Rj  = R - 4 * Ds; % UCJ circular flight radius [m]

% AP location (at origin)
qa = zeros(3, 1);

% UAV initial positions (on x-axis at respective radii)
qr_I = [Rr; 0; Hr];
qj_I = [Rj; 0; Hj];

%% ── UAV Kinematics ──────────────────────────────────────────────────────
Amax = 6;    % Maximum acceleration [m/s²]
Vmax = 25;   % Maximum speed [m/s]
dt   = 1;    % Time slot duration δt = T/N [s]

% Compute minimum number of time slots for feasible circular trajectory
j = 3;
while true
    theta = linspace(0, 2*pi, j);
    x = R * cos(theta);
    y = R * sin(theta);
    y(end) = 0;
    q = [x; y];
    tmp = diff(q, [], 2) / dt;
    v = [tmp, tmp(:, end)];
    j = j + 1;
    if max(norms(v)) <= Vmax
        break;
    end
end
Nmin = j;    % Minimum feasible time slots

%% ── User Equipment (UE) Placement ──────────────────────────────────────
numUsr = 5;  % Number of ground UEs (K)

% Generate uniformly distributed UE positions in annular region [R1, R2]
rng default   % Fixed seed for reproducibility
[xk, yk] = UsrRandDist(0, 0, R1, R2, numUsr);
qk = [xk; yk; zeros(1, numUsr)];  % 3D UE coordinates (z=0, ground level)
