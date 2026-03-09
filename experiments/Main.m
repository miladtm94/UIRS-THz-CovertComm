%% Main.m
% ==========================================================================
% MAIN EXPERIMENT RUNNER
% Runs all algorithms and generates all figures from the paper:
%
%   "Aerial Intelligent Reflecting Surface Enabled Terahertz Covert
%    Communications in Beyond-5G Internet of Things"
%   M. Tatar Mamaghani and Y. Hong, IEEE IoT Journal.
%
% SECTIONS (navigate using MATLAB's section runner or run sequentially):
%   Section 1 : Run Proposed, Benchmark I, Benchmark II (T=30, eps=0.01)
%   Section 2 : Generate convergence figures (Figs. 4-6)
%   Section 3 : Generate detection rate and velocity figures (Figs. 7-9)
%   Section 4 : mAEE vs. number of IRS elements (Fig. 10)
%   Section 5 : mAEE vs. carrier frequency (Fig. 11)
%
% RUNTIME WARNING:
%   Running this entire script may take several hours. Individual sections
%   can be run independently after generating the required .mat files.
%
% OUTPUT FILES:
%   All results saved to: data/results/
%   All figures saved to: figures/  (both .fig and .eps formats)
%
% Prerequisites:
%   CVX + MOSEK installed and configured.
%   All paths added via: addpath(genpath('../src')); addpath('../simulations');
% ==========================================================================

clc; clear all; close all;

% ── Path setup ─────────────────────────────────────────────────────────
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repoRoot, 'src')));
addpath(fullfile(repoRoot, 'simulations'));
resultsDir = fullfile(repoRoot, 'data', 'results');
figuresDir = fullfile(repoRoot, 'figures');
if ~exist(resultsDir,'dir'), mkdir(resultsDir); end
if ~exist(figuresDir,'dir'), mkdir(figuresDir); end

%% ═════════════════════════════════════════════════════════════════════════
%% SECTION 1: Run all algorithms (T=30, epsilon=0.01)
%% ═════════════════════════════════════════════════════════════════════════

T = 30; epsilon = 0.01; N = 30;
filename = fullfile(resultsDir, sprintf('myResults_T%deps%d.mat', T, 100*epsilon));

% ── Proposed JTCD ─────────────────────────────────────────────────────
clearvars -except filename T resultsDir figuresDir epsilon N
SystemParams; Proposed;
Itr_Prop = Itr; ATR_Prop = ATR; APC_Prop = APC; AEE_Prop = AEE;
save(filename, 'Itr_Prop', 'ATR_Prop', 'APC_Prop', 'AEE_Prop');
saveas(fig1, fullfile(figuresDir, sprintf('TrjProp_%deps%d', T, 100*epsilon)), 'fig');
saveas(fig1, fullfile(figuresDir, sprintf('TrjProp_%deps%d', T, 100*epsilon)), 'epsc');
close all;

% ── Benchmark I: Fixed Trajectory CD ─────────────────────────────────
clearvars -except filename T resultsDir figuresDir epsilon N
SystemParams; Benchmark1;
Itr_RAD = Itr; ATR_RAD = ATR; APC_RAD = APC; AEE_RAD = AEE;
save(filename, 'Itr_RAD', 'ATR_RAD', 'APC_RAD', 'AEE_RAD', '-append');
saveas(fig1, fullfile(figuresDir, sprintf('TrjBenchmarkI_%deps%d', T, 100*epsilon)), 'fig');
saveas(fig1, fullfile(figuresDir, sprintf('TrjBenchmarkI_%deps%d', T, 100*epsilon)), 'epsc');
close all;

% ── Benchmark II: Trajectory Design Only ─────────────────────────────
clearvars -except filename T resultsDir figuresDir epsilon N
SystemParams; Benchmark2;
Itr_JTD = Itr; ATR_JTD = ATR; APC_JTD = APC; AEE_JTD = AEE;
save(filename, 'Itr_JTD', 'ATR_JTD', 'APC_JTD', 'AEE_JTD', '-append');
saveas(fig1, fullfile(figuresDir, sprintf('TrjBenchmarkII_%deps%d', T, 100*epsilon)), 'fig');
saveas(fig1, fullfile(figuresDir, sprintf('TrjBenchmarkII_%deps%d', T, 100*epsilon)), 'epsc');
close all;

%% ═════════════════════════════════════════════════════════════════════════
%% SECTION 2: Convergence Figures (Figs. 4-6)
%% ═════════════════════════════════════════════════════════════════════════

clc; close all;
T = 30; epsilon = 0.01;
filename = fullfile(resultsDir, sprintf('myResults_T%deps%d.mat', T, 100*epsilon));
SystemParams; Flightconstants;
load(filename);

% Plotting style definitions
colours   = {'#0072BD','#77AC30','#A2142F','#7E2F8E','#D95319','#FF00FF','#EDB120'};
marker    = {'s','h','o','+','x','p','|'};
linestyle = {'-','--',':','-.'};
lgd_name  = {'Proposed JTCD', 'Benchmark I - CD', 'Benchmark II - TD'};

AEE_set = {AEE_Prop, AEE_RAD, AEE_JTD};
ATR_set = {ATR_Prop, ATR_RAD, ATR_JTD};
APC_set = {APC_Prop, APC_RAD, APC_JTD};
Itr_set = {Itr_Prop, Itr_RAD, Itr_JTD};

% Determine max iteration count for x-axis alignment
for i = 1:length(AEE_set), lengthArr(i) = numel(AEE_set{i}); end
l_offset = max(lengthArr);
scale    = 1e8 / ((Po+Pi)*2);

% --- Fig. 4: mAEE vs. iteration index ---
fig1 = figure(1); set(fig1,'defaultLegendAutoUpdate','off');
legend(lgd_name{1:length(AEE_set)}, 'AutoUpdate','off');
for i = 1:length(AEE_set)
    AEE_i = AEE_set{i};
    plt(i) = plot(0:l_offset-1, scale*[AEE_i, repmat(AEE_i(end),1,l_offset-length(AEE_i))]/1e3, '-or');
    plt(i).LineStyle = linestyle{mod(i,4)+1}; plt(i).Marker = marker{i};
    plt(i).MarkerSize = 5; plt(i).MarkerIndices = 1:3:l_offset;
    plt(i).Color = colours{i}; plt(i).LineWidth = 2; hold on;
end
legend(lgd_name); xlabel('Iteration index');
ylabel('minimum Average Energy Efficiency [kbits/Joule]');
xlim([0,l_offset-1]); grid on;
saveas(fig1, fullfile(figuresDir, sprintf('AEE_Itr_%deps%d', T, 100*epsilon)), 'fig');
saveas(fig1, fullfile(figuresDir, sprintf('AEE_Itr_%deps%d', T, 100*epsilon)), 'epsc');

% --- Fig. 5: mACT vs. iteration index ---
fig3 = figure(3);
for i = 1:3
    ATR_i = ATR_set{i};
    plt(i) = plot(0:l_offset-1, 1e8*[ATR_i, repmat(ATR_i(end),1,l_offset-length(ATR_i))]/1e6, '-or');
    plt(i).LineStyle = linestyle{mod(i,4)+1}; plt(i).Marker = marker{i};
    plt(i).MarkerSize = 5; plt(i).MarkerIndices = 1:3:l_offset;
    plt(i).Color = colours{i}; plt(i).LineWidth = 2; hold on;
end
xlabel('Iteration index'); ylabel('minimum Average Covert Throughput [Mbps]');
legend(lgd_name,'Location','northwest'); xlim([0,l_offset-1]); grid on;
saveas(fig3, fullfile(figuresDir, sprintf('mACT_Itr_%deps%d', T, 100*epsilon)), 'fig');
saveas(fig3, fullfile(figuresDir, sprintf('mACT_Itr_%deps%d', T, 100*epsilon)), 'epsc');

% --- Fig. 6: APC vs. iteration index ---
fig2 = figure(2);
for i = 1:3
    APC_i = APC_set{i};
    plt(i) = plot(0:l_offset-1, (Po+Pi)*2*[APC_i, repmat(APC_i(end),1,l_offset-length(APC_i))]/1e3, '-or');
    plt(i).LineStyle = linestyle{mod(i,4)+1}; plt(i).Marker = marker{i};
    plt(i).MarkerSize = 5; plt(i).MarkerIndices = 1:3:l_offset;
    plt(i).Color = colours{i}; plt(i).LineWidth = 2; hold on;
end
xlabel('Iteration index'); ylabel('Average Power Consumption [kW]');
legend(lgd_name,'Location','southwest'); xlim([0,l_offset-1]); grid on;
saveas(fig2, fullfile(figuresDir, sprintf('APC_Itr_%deps%d', T, 100*epsilon)), 'fig');
saveas(fig2, fullfile(figuresDir, sprintf('APC_Itr_%deps%d', T, 100*epsilon)), 'epsc');

%% ═════════════════════════════════════════════════════════════════════════
%% SECTION 3: Detection Rate, Velocity, and Power Figures (Figs. 7-9)
%% ═════════════════════════════════════════════════════════════════════════

clc; close all;
tt = linspace(0, T, T);

% --- Fig. 7: Min. error detection rate over time slots ---
fig4 = figure(4); N = T;
for S = 1:length(Itr_set)
    Itr = Itr_set{S};
    alpha = Itr(end).alpha; pa = Itr(end).pa; pj = Itr(end).pj;
    qr = Itr(end).Qr{1}; qj = Itr(end).Qj{1};
    Channels;
    idx = zeros(1, N);
    for n = 1:N, [~, idx(n)] = max(alpha(:,n)); end
    mYzeta = zeros(numUsr, N, S);
    for n = 1:N
        for k = 1:numUsr
            if (k == idx(n))
                zeta = 1 - repmat(pa(n), numUsr-1, 1) .* G_auw(:,k,n) ...
                           ./ repmat(pj(n), numUsr-1, 1) ./ G_jw(:,k,n);
                j = 1;
                for jj = 1:numUsr
                    if (jj ~= k)
                        mYzeta(jj, n, S) = zeta(j); j = j + 1;
                    end
                end
            end
        end
    end
    subplot(3,1,S); mycol = {'k', colours{1}, colours{2}, 'r', 'b'};
    for j = 1:numUsr
        yplot = reshape(mYzeta(j,:,S), 1, N); yplot(yplot==0) = nan;
        plt_h = plot(tt, yplot, 'LineStyle', linestyle{mod(j,4)+1}, ...
                     'Color', mycol{j}, 'Marker', marker{j}, 'LineWidth', 1.5, 'MarkerSize', 5);
        hold on; lgd = legend; lgd.String(j) = {sprintf('#%d', j)};
    end
    title(lgd_name{S}); grid on; hold on;
end
han = axes(fig4,'visible','off');
han.XLabel.Visible = 'on'; han.YLabel.Visible = 'on';
ylabel(han,'Minimum Error Detection Rate','FontSize',12,'Interpreter','latex');
xlabel(han,'Time slot [s]','FontSize',12,'Interpreter','latex');
saveas(fig4, fullfile(figuresDir, sprintf('minZeta_Time_%deps%d', T, 100*epsilon)), 'fig');
saveas(fig4, fullfile(figuresDir, sprintf('minZeta_Time_%deps%d', T, 100*epsilon)), 'epsc');

% --- Figs. 8-9: UIRS and UCJ velocity profiles ---
close all;
fig5 = figure(5); fig6 = figure(6);
for i = 1:length(Itr_set)
    Itr = Itr_set{i};
    vr = Itr(end).Qr{2}; vj = Itr(end).Qj{2};
    figure(5); plt(i) = plot(tt, norms(vr));
    plt(i).LineStyle = linestyle{i}; plt(i).Marker = marker{i};
    plt(i).Color = colours{i}; plt(i).LineWidth = 2; plt(i).MarkerSize = 5; hold on;
    figure(6); qlt(i) = plot(tt, norms(vj));
    qlt(i).LineStyle = linestyle{i}; qlt(i).Marker = marker{i};
    qlt(i).Color = colours{i}; qlt(i).LineWidth = 2; qlt(i).MarkerSize = 5; hold on;
end
figure(5); ylabel('UIRS velocity [m/s]','Interpreter','latex'); xlabel('Time [s]','Interpreter','latex');
legend(lgd_name,'Location','southeast'); xlim([0 T]); grid on;
saveas(fig5, fullfile(figuresDir, sprintf('UIRSvel_Time_%deps%d', T, 100*epsilon)), 'fig');
saveas(fig5, fullfile(figuresDir, sprintf('UIRSvel_Time_%deps%d', T, 100*epsilon)), 'epsc');
figure(6); ylabel('UCJ velocity [m/s]','Interpreter','latex'); xlabel('Time [s]','Interpreter','latex');
legend(lgd_name,'Location','southeast'); xlim([0 T]); grid on;
saveas(fig6, fullfile(figuresDir, sprintf('UCJvel_Time_%deps%d', T, 100*epsilon)), 'fig');
saveas(fig6, fullfile(figuresDir, sprintf('UCJvel_Time_%deps%d', T, 100*epsilon)), 'epsc');

%% ═════════════════════════════════════════════════════════════════════════
%% SECTION 4: mAEE vs. Number of IRS Elements (Fig. 10)
%% ═════════════════════════════════════════════════════════════════════════

clc; clear all; close all;
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repoRoot,'src'))); addpath(fullfile(repoRoot,'simulations'));
resultsDir = fullfile(repoRoot,'data','results');
figuresDir = fullfile(repoRoot,'figures');

T = 30; epsilon = 0.01; N = 30;
LL = [1 2 3 4 5 6 7 8 9; 2 3 4 5 6 7 8 9 10];  % Lx, Ly pairs
AEE_IRS = zeros(4, length(LL));

for ii = 1:length(LL)
    clearvars -except T epsilon N AEE_IRS LL ii resultsDir figuresDir repoRoot
    addpath(genpath(fullfile(repoRoot,'src'))); addpath(fullfile(repoRoot,'simulations'));
    SystemParams; Lx = LL(1,ii); Ly = LL(2,ii); L = Lx*Ly;
    Proposed; AEE_IRS(1,ii) = AEE(end); close all;

    clearvars -except T epsilon N AEE_IRS LL ii resultsDir figuresDir repoRoot
    addpath(genpath(fullfile(repoRoot,'src'))); addpath(fullfile(repoRoot,'simulations'));
    SystemParams; Lx = LL(1,ii); Ly = LL(2,ii); L = Lx*Ly;
    Benchmark1; AEE_IRS(2,ii) = AEE(end); close all;

    clearvars -except T epsilon N AEE_IRS LL ii resultsDir figuresDir repoRoot
    addpath(genpath(fullfile(repoRoot,'src'))); addpath(fullfile(repoRoot,'simulations'));
    SystemParams; Lx = LL(1,ii); Ly = LL(2,ii); L = Lx*Ly;
    Benchmark2; AEE_IRS(3,ii) = AEE(end); close all;

    clearvars -except T epsilon N AEE_IRS LL ii resultsDir figuresDir repoRoot
    addpath(genpath(fullfile(repoRoot,'src'))); addpath(fullfile(repoRoot,'simulations'));
    SystemParams; Lx = LL(1,ii); Ly = LL(2,ii); L = Lx*Ly;
    ite_index = 1; Feasible_Init; AEE_IRS(4,ii) = AEE(end); close all;
end

Flightconstants;
scale = 1e8 / ((Po+Pi)*2);
AEE_IRS_scaled = scale * AEE_IRS / 1e3;
save(fullfile(resultsDir,'AEE_IRS_scaled.mat'), 'AEE_IRS_scaled', 'LL');

colours   = {'#0072BD','#77AC30','#A2142F','#7E2F8E','#D95319','#FF00FF','#EDB120'};
marker    = {'s','h','o','+','x','p','|'};
linestyle = {'-','--',':','-.'};
IRSnum    = LL(1,:).*LL(2,:);
figX = figure(10);
lgd_name4 = {'Proposed JTCD','Benchmark I - CD','Benchmark II - TD','Benchmark III - IFTR'};
for i = 1:4
    qq(i) = plot(IRSnum, AEE_IRS_scaled(i,:));
    qq(i).LineStyle = linestyle{i}; qq(i).Marker = marker{i};
    qq(i).MarkerSize = 5; qq(i).Color = colours{i}; qq(i).LineWidth = 2; hold on;
end
xlabel('Number of IRS elements','Interpreter','latex');
ylabel('minimum Average Energy Efficiency [kbits/Joule]','Interpreter','latex');
legend(lgd_name4,'Location','northwest'); grid on;
saveas(figX, fullfile(figuresDir, sprintf('IRSnum_%deps%d', T, 100*epsilon)), 'fig');
saveas(figX, fullfile(figuresDir, sprintf('IRSnum_%deps%d', T, 100*epsilon)), 'epsc');

%% ═════════════════════════════════════════════════════════════════════════
%% SECTION 5: mAEE vs. Carrier Frequency (Fig. 11)
%% ═════════════════════════════════════════════════════════════════════════

clc; clear all; close all;
repoRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repoRoot,'src'))); addpath(fullfile(repoRoot,'simulations'));
resultsDir = fullfile(repoRoot,'data','results');
figuresDir = fullfile(repoRoot,'figures');

T = 30; epsilon = 0.01; N = 30;
fcVec = 1e9 * linspace(275, 400, 6); fcVec = fcVec(2:end-1);  % 4 inner freq points

% Compute kappa for each frequency point using absorption formula
% (requires WaterVapor_Loss.m helper or manual computation)
% kfVec = [ka(fcVec(1),T,p), ka(fcVec(2),T,p), ...];
% Placeholder: replace with actual values from ka() function
kfVec = zeros(1, length(fcVec));  % REPLACE WITH ACTUAL VALUES
AEE_kappa = zeros(4, length(fcVec));

for ii = 1:length(fcVec)
    clearvars -except T epsilon N AEE_kappa fcVec kfVec ii resultsDir figuresDir repoRoot
    addpath(genpath(fullfile(repoRoot,'src'))); addpath(fullfile(repoRoot,'simulations'));
    fc = fcVec(ii); kf = kfVec(ii);
    SystemParams_fc; Proposed; AEE_kappa(1,ii) = AEE(end); close all;

    clearvars -except T epsilon N AEE_kappa fcVec kfVec ii resultsDir figuresDir repoRoot fc kf
    addpath(genpath(fullfile(repoRoot,'src'))); addpath(fullfile(repoRoot,'simulations'));
    fc = fcVec(ii); kf = kfVec(ii);
    SystemParams_fc; Benchmark1; AEE_kappa(2,ii) = AEE(end); close all;

    clearvars -except T epsilon N AEE_kappa fcVec kfVec ii resultsDir figuresDir repoRoot fc kf
    addpath(genpath(fullfile(repoRoot,'src'))); addpath(fullfile(repoRoot,'simulations'));
    fc = fcVec(ii); kf = kfVec(ii);
    SystemParams_fc; Benchmark2; AEE_kappa(3,ii) = AEE(end); close all;

    clearvars -except T epsilon N AEE_kappa fcVec kfVec ii resultsDir figuresDir repoRoot fc kf
    addpath(genpath(fullfile(repoRoot,'src'))); addpath(fullfile(repoRoot,'simulations'));
    fc = fcVec(ii); kf = kfVec(ii);
    SystemParams_fc; ite_index = 1; Feasible_Init; AEE_kappa(4,ii) = AEE(end); close all;
end

Flightconstants;
scale = 1e8 / ((Po+Pi)*2);
AEE_kappa_scaled = scale * AEE_kappa / 1e3;
save(fullfile(resultsDir,'mAEE_freq.mat'), 'AEE_kappa_scaled', 'fcVec', 'kfVec');

figV = figure(11);
lgd_name4 = {'Proposed JTCD','Benchmark I - CD','Benchmark II - TD','Benchmark III - IFTR'};
colours = {'#0072BD','#77AC30','#A2142F','#7E2F8E'};
linestyle = {'-','--',':','-.'}; marker = {'s','h','o','+'};
for i = 1:4
    qq(i) = plot(fcVec/1e9, AEE_kappa_scaled(i,:));
    qq(i).LineStyle = linestyle{i}; qq(i).Marker = marker{i};
    qq(i).MarkerSize = 5; qq(i).Color = colours{i}; qq(i).LineWidth = 2; hold on;
end
xlabel('Carrier frequency [GHz]','Interpreter','latex');
ylabel('minimum Average Energy Efficiency [kbits/Joule]','Interpreter','latex');
legend(lgd_name4,'Location','northeast'); grid on;
saveas(figV, fullfile(figuresDir, sprintf('mAEEfreq_%deps%d', T, 100*epsilon)), 'fig');
saveas(figV, fullfile(figuresDir, sprintf('mAEEfreq_%deps%d', T, 100*epsilon)), 'epsc');

fprintf('\n All experiments complete. Figures saved to: %s\n', figuresDir);
