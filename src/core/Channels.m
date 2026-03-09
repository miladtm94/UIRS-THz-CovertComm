%% Channels.m
% ==========================================================================
% THz CHANNEL GAIN COMPUTATION
%
% Computes all channel power gains for the current UIRS and UCJ trajectories:
%
%   G_auk(k,n)     : AP-UIRS-UE_k cascaded channel gain in slot n
%                    (includes THz path loss + molecular absorption)
%   G_jk(k,n)      : UCJ-to-UE_k direct channel gain in slot n
%   G_auw(m,k,n)   : AP-UIRS-Willie_m gain when UE k is scheduled in slot n
%   G_jw(m,k,n)    : UCJ-to-Willie_m gain when UE k is scheduled in slot n
%
% Channel Model (THz, Eq. 7 and 9 in paper):
%   h_jk[n]  = (λ_c/4π) · d_jk^(-ρ/2) · exp(-jωd_jk/c) · exp(-κd_jk/2)
%   h̃_ark[n] = (λ_c/8π√π) · (d_ar·d_rk)^(-ρ/2) · exp(-j·...) · exp(-κd_ark/2)
%   h_ark[n] = e_k†[n] · Φ[n] · e_a[n] · h̃_ark[n]
%
% Note: In the current implementation, IRS beamforming is absorbed into the
% channel gain as G_auk = |h̃_ark|² · L² (maximum coherent gain with optimal
% phase alignment), consistent with the closed-form approach used in the
% paper for trajectory and power sub-problems.
%
% Inputs (from workspace):
%   qr, qj, qk, qa      : 3D positions of UIRS, UCJ, UEs, AP
%   hr0, hj0, kf, pl    : channel model parameters
%   lambda_c             : carrier wavelength
%   Lx, Ly, L, delta_x, delta_y : IRS geometry parameters
%   N, numUsr            : number of time slots and UEs
%
% Outputs (to workspace):
%   Ea(L,N)           : AP-to-IRS steering vectors
%   Ek(K,L,N)         : IRS-to-UE steering vectors
%   G_auk(K,N)        : AP-UIRS-UE cascaded channel gains
%   G_jk(K,N)         : UCJ-UE channel gains
%   G_auw(K-1,K,N)    : AP-UIRS-Willie channel gains
%   G_jw(K-1,K,N)     : UCJ-Willie channel gains
% ==========================================================================

%% ── IRS Steering Vectors ────────────────────────────────────────────────
% AP-to-IRS array steering vector (Eq. 11): e_a[n] = exp(-j*theta_l[n])
theta = zeros(L, N);
for lx = 1:Lx
    for ly = 1:Ly
        l = (ly-1)*Lx + lx;
        delta_lxly = [(lx-1)*delta_x; (ly-1)*delta_y; 0];
        temp = (2*pi*(qr-qa)' * delta_lxly) ./ lambda_c ./ reshape(norms(qr-qa), N, 1);
        theta(l, :) = wrapTo2Pi(temp');
    end
end
Ea = exp(-1i * theta);  % (L × N) steering vectors from AP to IRS

% IRS-to-UE array steering vectors (Eq. 13): e_k[n] = exp(-j*beta_l_k[n])
beta = zeros(numUsr, L, N);
for usr = 1:numUsr
    for lx = 1:Lx
        for ly = 1:Ly
            l = (ly-1)*Lx + lx;
            delta_lxly = [(lx-1)*delta_x; (ly-1)*delta_y; 0];
            temp = (2*pi*(qk(:,usr)-qr)' * delta_lxly) ./ lambda_c ./ ...
                   reshape(norms(qk(:,usr)-qr), N, 1);
            beta(usr, l, :) = wrapTo2Pi(temp');
        end
    end
end
Ek = exp(-1i * beta);   % (K × L × N) steering vectors from IRS to UEs

%% ── AP-UIRS-UE Cascaded Channel Gains (G_auk) ──────────────────────────
beam   = zeros(numUsr, N);
G_auk  = zeros(numUsr, N);
G_jk   = zeros(numUsr, N);

for usr = 1:numUsr
    % AP-to-UIRS and UIRS-to-UE distances
    d_ar  = norms(qa - qr);           % AP to UIRS (1×N)
    d_rk  = norms(qk(:,usr) - qr);   % UIRS to UE k (1×N)
    d_ark = d_ar + d_rk;              % Total cascaded distance

    % AP-UIRS-UE_k complex channel (Eq. 9)
    h_ark = hr0 .* exp(-1i*2*pi*d_ark/lambda_c) ...
            .* exp(-kf.*d_ark/2) ...
            .* (d_ar .* d_rk).^(-pl/2);

    % Channel gain with maximum IRS beamforming gain (L² coherent combining)
    for n = 1:N
        G_auk(usr, n) = (abs(h_ark(n))^2) * (L^2);
    end

    % UCJ-to-UE_k direct channel gain (Eq. 7)
    d_jk  = norms(qk(:,usr) - qj);
    h_jk  = hj0 .* exp(-1i*2*pi*d_jk/lambda_c) ...
            .* exp(-kf.*d_jk/2) .* (d_jk).^(-pl/2);
    G_jk(usr, :) = abs(h_jk .* conj(h_jk));
end

%% ── Willie Channel Gains (G_auw, G_jw) ─────────────────────────────────
% G_auw(m,k,n): gain from AP-UIRS to Willie m, when UE k is scheduled
% G_jw(m,k,n):  gain from UCJ to Willie m, when UE k is scheduled
G_auw = zeros(numUsr-1, numUsr, N);
G_jw  = zeros(numUsr-1, numUsr, N);

for n = 1:N
    for k = 1:numUsr
        j = 1;
        for w = 1:numUsr
            if (w ~= k)   % w is a Willie (unscheduled UE)
                d_ar  = norms(qa - qr(:,n));
                d_rm  = norms(qk(:,w) - qr(:,n));
                d_arm = d_ar + d_rm;
                h_arm = hr0 .* exp(-1i*2*pi*d_arm/lambda_c) ...
                        .* exp(-kf.*d_arm/2) .* (d_ar.*d_rm).^(-pl/2);
                G_auw(j, k, n) = abs(h_arm)^2 * (L^2);

                d_jm  = norms(qk(:,w) - qj(:,n));
                h_jm  = hj0 .* exp(-1i*2*pi*d_jm/lambda_c) ...
                        .* exp(-kf.*d_jm/2) .* (d_jm).^(-pl/2);
                G_jw(j, k, n)  = abs(h_jm .* conj(h_jm));
                j = j + 1;
            end
        end
    end
end
