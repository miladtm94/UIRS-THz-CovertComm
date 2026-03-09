%% --- Flight power constants for both UAVs
omega = 300; % rad/s
r = 0.4; % m
rho = 1.225;  % kg/m3
s = 0.05;   %
A = 0.503;  %m^2
v0 = 4.03;  % m/s
d0 = 0.6; 
delta_coeff = 0.012;
k_corr = 0.1;
Wgt = 20;
Po = delta_coeff*rho*s*A*omega^3*r^3/8; %79.86
Pi = (1+k_corr)*Wgt^(3/2)/(sqrt(2*rho*A)); %88.63
c0 = 3/(omega^2)/(r^2);
c1 = 0.5*d0*rho*s*A;
c2 = 1/2/(v0^2);
