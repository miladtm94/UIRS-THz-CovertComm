%% -- Feasible Initialization --------------%%
Flightconstants

%%
% rep = floor(N/Nmin);
theta = linspace(0,2*pi, N);
xr = qr_I(1)*cos(theta);
yr = qr_I(1)*sin(theta);
yr(end)=0;
qr = [sign(xr).*ceil(abs(xr));sign(yr).*ceil(abs(yr));repmat(Hr,1,N)];

tmp = diff(qr,[],2)/dt;
vr = [tmp(1:2,:), tmp(1:2,end)];

xj = qj_I(1)*cos(theta);
yj = qj_I(1)*sin(theta);
yj(end)=0;
qj = [sign(xj).*floor(abs(xj));sign(yj).*floor(abs(yj));repmat(Hj,1,N)];

tmp = diff(qj,[],2)/dt;
vj = [tmp(1:2,:), tmp(1:2,end)];

% % Inputs - change as desired
% pad_size = rep;
% % Process the inputs
% initial_size = size(v,2); 
% final_size = initial_size + (pad_size * (initial_size));
% % Assign a variable with zeros for the final size of the output
% vr = zeros(2,final_size); % modified to be a row instead of a column
% vr(:,1:pad_size+1:end) = v;
% vj=vr; % display the result

%%

% flag1 = true;
% i = 0;
% while(flag1)
% 
% %strInput = sprintf('Enter the value of N (integer multiple #%d):\t', numUsr);
% %N = input(strInput);
%     N = 80;
% if (N >= Nmin)
%     theta = linspace(0,2*pi, N);
%     xr = qr_I(1)*cos(theta);
%     yr = qr_I(1)*sin(theta);
%     yr(end) = qr_I(2);
% %     vxr = diff(xr)/dt;
% %     vxr = [vxr, vxr(end)];
% %     vyr = diff(yr)/dt;
% %     vyr = [vyr, vyr(end)];
%     qr = [xr;yr;repmat(Hr,1,N)];
%     
%     tmpv=(qr(1:2,2:N) - (qr(1:2,1:N-1)))/dt;
%     vr = [tmpv, tmpv(end)];
% 
%     xj = (qj_I(1))*cos(theta);
%     yj = (qj_I(1))*sin(theta);
%     yj(end) = qj_I(2);
%     qj = [xj;yj;repmat(Hj,1,N)];
%     tmpv=(qj(1:2,2:N) - (qj(1:2,1:N-1)))/dt;
%     vj = [tmpv, tmpv(end)];
%     
%     flag1 = false;
%  
% else
%     disp('Increase N to make circular flight plausible :)');
%     
% end
% end

%%
% ==== 2D visulazation
fprintf("Maximum UAV-IRS's velocity per time slot is %2.3f m/s \n", max(norms(diff(qr,1,2)/dt)));
fprintf("Maximum UAV-IRS's acceleration per time slot is %2.3f m/s^2 \n", max(norms(diff(vr,1,2))));
fprintf("Maximum UAV-CJ's velocity per time slot is %2.3f m/s \n", max(norms(diff(qj,1,2)/dt)));
fprintf("Maximum UAV-CJ's acceleration per time slot is %2.3f m/s^2 \n", max(norms(diff(vj,1,2))));
% Visulization

% uavTrj(xu,yu, fig1_ax, lgd, 1, 1)

%%
% =====  User scheduling initializaiton
% alpha = zeros(numUsr, N);
% for i=1:N
%     [~, idx(i)] = min(vecnorm(qr(:,i)-qk));
%     alpha(idx(i),i)=1; 
% 
% end

%% Random user scheduling
idx = randi([1,numUsr],1,N);
alpha = zeros(numUsr,N);
for n=1:N
    alpha(idx(n),n)=1;
end
% sum(alpha,2)
% c = floor(N/numUsr);
% base = [ones(1,c) zeros(1,N-c)];
% alpha = zeros(numUsr, N);
% for i=1:numUsr
%     alpha(i,:) = base;
%     base = circshift(base,[1, c]);
% end
% 
% for n=1:N
%    [~, idx(n)] = max(alpha(:,n));
% end

%%
Phi = exp(1i*2*pi*rand(L,numUsr,N));

theta = zeros(L,N);
for lx = 1:Lx
    for ly = 1:Ly
        l=(ly-1)*Lx+lx;
        delta_lxly = [(lx-1)*delta_x; (ly-1)*delta_y; 0];
        temp = (2*pi*(qr-qa)'*delta_lxly)./lambda_c./reshape(norms(qr-qa), N, 1);
         theta(l, :) = wrapTo2Pi(temp');
    end
end
Ea = exp(-1i*theta);

beta = zeros(numUsr, L, N);
for usr=1:numUsr
    for lx = 1:Lx
        for ly = 1:Ly
            l=(ly-1)*Lx+lx;
            delta_lxly = [(lx-1)*delta_x; (ly-1)*delta_y; 0];
            temp = (2*pi*(qk(:,usr)-qr)'*delta_lxly)./lambda_c./reshape(norms(qk(:,usr)-qr), N, 1);
            beta(usr, l, :) = wrapTo2Pi(temp');
        end
    end
        
end

for k=1:numUsr
    Phi(:,k,n) = exp(-1i*(theta(:,n) - beta(k,:,n)'));
end

%%

% =====  Transmit power initializaiton

pa = (P_tot/N/4)*ones(1,N); 
% pa = (P_tot/N/8)*ones(1,N); 

pj = (P_tot/N/2)*ones(1,N);
% r = a + (b-a).*rand(N,1)

AEE_Calc

