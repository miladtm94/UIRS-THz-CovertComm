%------ 2D Visulaization of network locations

theta = linspace(0,2*pi, N);

xm = R1*cos(theta);
ym = R1*sin(theta);

xM = R2*cos(theta);
yM = R2*sin(theta);

fig1 = figure(1);
fig1.Name = 'Traj';
fig1_ax = axes;
plot(fig1_ax, xm, ym, ':k',  'LineWidth', 1);

hold on;
plot(fig1_ax, xM, yM, '--k', 'LineWidth', 1 );

%---- Linear Traj
% xa = [1 1    2  2.5 2    3     4   3.5 1.5 3.5   ];
% ya = [1 -0.5  -1 -2  1   2.5  -1.5 1 0 0];
% xu = linspace(0,x_final, N);
% yu = zeros(1,N);

%fprintf("displacement = %.3fm\n",max(vecnorm([diff(xu); diff(yu)])))

hold on
%%%%%%  Users  Location  
U=scatter(fig1_ax, xk,yk);
hold on
for i=1:length(xk)
UsrLabel{i} = sprintf('\n#%d',i);
end
text(xk,yk,UsrLabel,'HorizontalAlignment','center')
U.LineWidth = 2;
U.MarkerEdgeColor = 'b';
U.MarkerFaceColor = 'm';

hold on
%%%%%%  Source Location  
B=scatter(fig1_ax, qa(1),qa(2));
B.LineWidth = 2;
B.MarkerEdgeColor = 'g';
B.MarkerFaceColor = 'y';
A = sprintf('\n AP');
text(qa(1),qa(2),A,'HorizontalAlignment','center')


hold on
%%%%%%  Initial and Final UAV's Location  
F=scatter(fig1_ax, qr_I(1),qr_I(2),'H');
F.LineWidth = 2;
F.MarkerEdgeColor = 'r';
F.MarkerFaceColor = 'c';
%  txt = '\Uparrow UAV-IRS';
%  text(qr_I(1)+5,qr_I(2)+5,txt)
 
FF=scatter(fig1_ax,qj_I(1),qj_I(2),'H');
FF.LineWidth = 2;
FF.MarkerEdgeColor = 'b';
FF.MarkerFaceColor = 'g';
% txt = 'UCJ-Init';
% text(qj_I(1),qj_I(2),txt,'HorizontalAlignment','left')


lgd = legend(fig1_ax, [F, FF],{'UIRS-Init', 'UCJ-Init' });
lgd.Location='northeast';
lgd.FontSize=10;
lgd.FontName='Times New Roman';
axis([-(1.2*R2) (1.2*R2) -(1.2*R2) (1.2*R2) ]);


hold on
xlabel ('x [m]')
ylabel ('y [m]')
