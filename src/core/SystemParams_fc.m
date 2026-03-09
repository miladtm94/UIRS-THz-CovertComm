%% Setting Parameters %%%
epsilon = 1e-2;  % covertness constraint
flag = 'false';
eps_algo = 1e-1;  % Algorithm's termination fractional error criteria
eps_frac_algo = 1e-2;  % dinklebach's termination criteria
%% --- communication constants
% fc = 300e9;  % 0.3THz
C = 3e8;
lambda_c = C/fc;
PSD = db2pow(-250);  % -190 dBm/Hz
BW = 10e9; 
delta2 = PSD*BW; % noise power at users
BW = 1e2;
hj0 = (lambda_c/4/pi)/sqrt(delta2);  % -71dB channel gain at 1 meter
hr0 = (lambda_c/(8*(sqrt(pi^3))))/sqrt(delta2);
% kf = 3.2094e-04; % molecular channel absorption
Pa_max = 1;
P_tot = 40;
Pj_max = 1;
pl = 2.3;
%% --- IRS-related parameters
delta_x = 1e-3;
delta_y = 1e-3;
Lx = 6;
Ly = 5;
L= Lx*Ly;
%% -- Geometric parameters and network topology
H = 50;
Ds = 0.1*H;
Hr  = H;
Hj = H;  % Hj = Hr + Ds
R  = 2*H; %100
R1 = 2*H-1; R2 = 4*H;

Rr = R; Rj = R - 4*Ds;
qa = zeros(3,1);
qr_I = [Rr; 0; Hr];
qj_I = [Rj; 0; Hj];


Amax = 6; %m/s^2
Vmax = 25; % m/s
dt = 1; % T/N = 1

% j= 3;
% while (true)
%     theta = linspace(0,2*pi, j);
%     x = R*cos(theta);
%     y = R*sin(theta);
%     y(end) = 0;
%     q = [x;y];
%     tmp= diff(q,[],2)/dt;
%     v = [tmp, tmp(:,end)];
%    j = j + 1;
%     if(max(norms(v))<= Vmax)
%         break;
%     end
% end
% Nmin = j


numUsr = 5;

rng default
[xk, yk]= UsrRandDist(0,0,R1,R2,numUsr);
qk = [xk;yk;zeros(1,numUsr)];

