
clc 
clear all
close all

filename = 'myResults.mat';  %N = 30

 SystemParams   % loading system parameters
 FeasiblePoint  % feasible initialization
main_SEEoptim_MIBCD
Itr_MIBCD = Itr;
SEE_MIBCD = SEE;
save(filename, 'Itr_MIBCD')
save(filename, 'SEE_MIBCD','-append')

 SystemParams   % loading system parameters
 FeasiblePoint  % feasible initialization
main_SEEoptim_SeqBCD
Itr_SeqBCD = Itr;
SEE_SeqBCD = SEE;
save(filename, 'Itr_SeqBCD', '-append')
save(filename, 'SEE_SeqBCD', '-append')

 SystemParams   % loading system parameters
 FeasiblePoint  % feasible initialization

main_SEEoptim_FixedTrj
Itr_FixedTrj = Itr;
SEE_FixedTrj = SEE;
save(filename, 'Itr_FixedTrj', '-append')
save(filename, 'SEE_FixedTrj', '-append')

 SystemParams   % loading system parameters
 FeasiblePoint  % feasible initialization

main_SEEoptim_FixedPow
Itr_FixedPow = Itr;
SEE_FixedPow = SEE;
save(filename, 'Itr_FixedPow', '-append')
save(filename, 'SEE_FixedPow', '-append')

 SystemParams   % loading system parameters
 FeasiblePoint  % feasible initialization

main_ASRoptim_NAPCL
Itr_ASRoptim = Itr;
SEE_ASRoptim = SEE;
save(filename, 'Itr_ASRoptim', '-append')
save(filename, 'SEE_ASRoptim', '-append')

%%
clc
close all

load(filename);

Itr_arr = {Itr_MIBCD, Itr_SeqBCD, Itr_FixedTrj, Itr_FixedPow, Itr_ASRoptim};
SEE_arr = {SEE_MIBCD, SEE_SeqBCD, SEE_FixedTrj, SEE_FixedPow, SEE_ASRoptim};



Itr_init = Itr_arr{1};
Qu = Itr_init(1).Trj';
m=0; 
Visulization

for i=1:length(Itr_arr)
    Itr = Itr_arr{i};
    l = length(SEE_arr{i});
    Qu = Itr(end).Trj';
    m=m+1;
    uavTrj(Qu(1,:),Qu(2,:), fig1_ax,lgd, m, l)
end

%% --- SEE comparison

figure
colours = {'k', '#0072BD', '#D95319', '#7E2F8E','#77AC30', '#A2142F'};
marker = {'none','^', '*', '+', 'x','p'}; 
lgd_name = {'SEE-MI',...
      'SEE-Seq',...
        'SEE-FTrj Trj. ',...
        'SEE-FPow Trj. ',...
        'ASR-Seq Trj.'}; 
linestyle = {'-','--',':','-.'};

for i=1:length(SEE_arr)-1
    SEE = SEE_arr{i};
    l_tem(i) = numel(SEE);    
end
l_offset = max(l_tem);

for i=1:length(SEE_arr)-1
    SEE = SEE_arr{i};
    plt(i) = plot(0:l_offset-1, BW*[SEE, repmat(SEE(end),1,l_offset-length(SEE))] /P_lim/1e6, '-or','LineWidth',2);
    plt(i).LineStyle = linestyle{mod(i,4)+1};
    plt(i).Marker= marker{i};
    plt(i).MarkerIndices = 1:3:l_offset;
    plt(i).Color = colours{i};
    plt(i).LineWidth = 2;
    hold on
end
xlabel('Iteration index');
ylabel('Minimum secrecy energy efficiency (Mbits/Joule)');
lgd = legend(lgd_name{1:length(SEE_arr)-1});
lgd.Location='northeast';
lgd.FontSize=10;
lgd.FontName='Times New Roman';

xlim([0,l_offset-1])

%% ---  Velocity and flight power comparisons


fig_Vel = figure(4);
colours = {'#0072BD', '#D95319', '#7E2F8E', '#A2142F', '#77AC30','#A2142F'};
marker = {'s', 'h', '^', 'p', 's', 'o'};
linestyle = {'-','--',':','-.'};
lgd_name = {'SEE-MI, Vel.', 'SEE-Seq, Vel.',  'SEE-FTrj, Vel.',...
        'SEE-FPow, Vel.', 'SEE-MI, IFPC','SEE-Seq, IFPC', ...
        'SEE-FTrj, IFPC', 'SEE-FPow, IFPC'}; 

tt = linspace(0,T,N);

for i=1:length(Itr_arr)-1
    Itr = Itr_arr{i};
    hold on
    yyaxis left
    pp(i)=plot(tt,vecnorm(Itr(end).Vel',2));
    pp(i).LineStyle = linestyle{i};
    pp(i).Marker= marker{i};
    pp(i).MarkerIndices = 1:10:N;
    pp(i).LineWidth = 2;  
    hold on
    yyaxis right
    qq(i) = plot(tt,Itr(end).prop);    
    qq(i).LineStyle = linestyle{i};
    qq(i).Marker= marker{i};
    qq(i).MarkerIndices = 1:10:N;
    qq(i).LineWidth = 2;
end


hold on
yyaxis left
ylabel("UAV's velocity [m/s]")

hold on
yyaxis right
ylabel("UAV's instantaneous flight power consumption [W]")
    
xlim([0 T]);
legend(lgd_name)
xlabel('Time [s]');

grid on


%% Tranmission powers  comparison for different scenarios
close all
colours = {'#0072BD', '#D95319', '#7E2F8E', '#A2142F', '#77AC30','#A2142F'};
marker = {'s', '^', 'p', 'd', 'o'};
linestyle = {'-','--',':','-.'};
title_name = {'SEE-MI', 'SEE-Seq', 'SEE-FTrj', 'SEE-FPow', 'ASR-Seq'}; 

tt = linspace(0,T,N);

for i=1:length(Itr_arr)
    fig_TxPow(i)= figure;
    Itr = Itr_arr{i};
    zeta = Itr(end).usrSch;
    zeta = zeta';
    Pu = Itr(end).uavPow;
    Pu = Pu';
    Pa = Itr(end).usrPow;
    Pa = Pa';
    Pb = Itr(end).bsPow;
    Pb = Pb';
    numUsr = size(zeta,1);
    zeta_n = (zeta>0.5);
    for j=1:numUsr
    hold on
        plt(j) = plot(tt, zeta_n(j,:).*(zeta(j,:).*Pa(j,:)));
        plt(j).LineStyle = linestyle{mod(j,4)+1};
        plt(j).Marker= marker{j};
        plt(j).MarkerIndices = 1:10:N;
        plt(j).LineWidth = 2;  
    end
    hold on
    plt(numUsr+1) = plot(tt,Pu, '-or');
    plt(numUsr+1).MarkerIndices = floor(linspace(1,N,20));
    plt(numUsr+1).LineWidth = 2;
    hold on
    plt(numUsr+2) = plot(tt,Pb, '--pb');
    plt(numUsr+2).MarkerIndices = floor(linspace(1,N,20));
    plt(numUsr+2).LineWidth = 2;
    
    Legend=cell(numUsr+2,1);
    for usr=1:numUsr
   Legend{usr}=strcat('P_', num2str(usr));
    end
    Legend{numUsr+1} = 'P_u';
    Legend{numUsr+2} = 'P_b';

    
    lgd =  legend(Legend);
    lgd.Location='northeast';
    lgd.FontSize=10;
    lgd.FontName='Times New Roman';

    ylabel('Transmit power [W]');
    xlabel('Time [s]')
    title (title_name{i});
    
    grid on
    hold off
    
    clear Pa Pb Pu zeta

end  

%% --- ASR vs APC for different schemes over the iteration index

close all
colours = {'#0072BD', '#D95319', '#7E2F8E', '#A2142F', '#77AC30','#A2142F'};
marker = {'s', '^', 'p', 'd', 'o'};
linestyle = {'-',':','-.'};
Legend_name = {'SEE-MI, mASR', 'SEE-Seq, mASR', 'SEE-FTrj, mASR', ...
               'SEE-FPow, mASR', 'ASR-Seq, mASR', ...
               'SEE-MI, AFPC', 'SEE-Seq, AFPC', 'SEE-FTrj, AFPC', ...
               'SEE-FPow, AFPC', 'ASR-Seq, AFPC'}; 

for i=1:length(SEE_arr)
    SEE = SEE_arr{i};
    l_tem(i) = numel(SEE);    
end
l_offset = max(l_tem);


for i=1:length(SEE_arr)
    APC=[];
    ASR = [];
    SEE = SEE_arr{i};
    Itr = Itr_arr{i};
    for j=1:l_tem(i)
        APC(j) = mean(Itr(j).prop);
    end
    ASR = BW*SEE.*APC/P_lim;
    yyaxis left
    hold on
    lplt(i) = plot(0:l_offset-1, [ASR, repmat(ASR(end),1,l_offset-length(ASR))]/1e9);
    lplt(i).LineStyle = linestyle{1};
    lplt(i).Marker= marker{i};
    lplt(i).MarkerIndices = i:3:l_offset;
    %lplt(i).Color = colours{i};
    lplt(i).LineWidth = 2;
    ylabel('Minimum average secrecy rate (Gbits/s)');

    yyaxis right
    hold on
    rplt(i) = plot(0:l_offset-1, [APC, repmat(APC(end),1,l_offset-length(APC))]);
    rplt(i).LineStyle = linestyle{2};
    rplt(i).Marker= marker{i};
    rplt(i).MarkerIndices = i:3:l_offset;
    %rplt(i).Color = colours{i};
    rplt(i).LineWidth = 2;
    ylabel("UAV's average flight power consumption (W)");
    
end
xlabel('Iteration index (l)');
lgd = legend(Legend_name);
lgd.Location='northeast';
lgd.FontSize=10;
lgd.FontName='Times New Roman';

xlim([0,l_offset-1])
    

%%
pName = "G:\My Drive\MyPhDResearch\ICC2021_NEW\Results";
dName = sprintf("SEEoptim_MIBCD%d", N);
mkdir (fullfile(pName, dName));
Fig = {fig_Trj, fig_SEE, fig_Pow, fig_Vel, fig_AP };
Names={'Traj', 'Conv', 'Pow', 'Vel', 'AP'};
for i=1:5
path = fullfile(pName,dName, sprintf('%s',Names{i}));
saveas(Fig{i},path, 'fig')
end

