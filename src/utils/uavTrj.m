function uavTrj(qr,qj, fig1_ax, lgd, mm)
    
    colours = {'#A2142F','#0072BD', '#D95319', '#7E2F8E','#77AC30', '#A2142F','#000000'};
    marker = {'<','v','o','h','d','^'};
    lgd_name = {sprintf('UIRS-Trj, ite = %d',mm-1), sprintf('UCJ-Trj, ite = %d',mm-1)};    
    hold on
    xu = qr(1,:);
    yu = qr(2,:);
    R = plot(fig1_ax, xu,yu);
    R.LineStyle = '-';
    R.Marker= marker{mod(mm,length(marker))+1};
    R.MarkerSize = 2.5;
    R.MarkerIndices = 2:length(xu)-1;
    R.Color = colours{mod(mm,length(colours))+1};
    R.LineWidth = 1.5;
    
    hold on
    xu = qj(1,:);
    yu = qj(2,:);
    R = plot(fig1_ax, xu,yu);
    R.LineStyle = ':';
    R.Marker= marker{mod(mm,length(marker))+1};
    R.MarkerSize = 2.5;
    R.MarkerIndices = 2:length(xu)-1;
    R.Color = colours{end-mod(mm,length(colours))};
    R.LineWidth = 1.5;
    
    hold off
    hFigure = findall(0,'type','figure','name','Traj');
    TrjcurvesCnt = findall(hFigure,'type','line');  
    
    mmm = 2*mm+1;
    lgd.String(mmm) = {lgd_name{1}};
    lgd.String(mmm+1) = {lgd_name{2}};
    legend
    
end

%  Legend=cell(N,1)
%  for iter=1:N
%    Legend{iter}=strcat('Your_Data number', num2str(iter));
%  end
%  legend(Legend)