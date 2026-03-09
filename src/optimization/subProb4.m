Channels

% G_jw = zeros(numUsr-1,N);
% g_jb = zeros(1,N);
% 
% for n=1:N
%     temp1 = G_jk(:,n);
%     g_jb(n) =  G_jk(idx(n),n);  
%     temp1(idx(n)) = [];
%     G_jw(:,n) = temp1; 
% end


% A ln(1+B z * )
A = alpha*BW/log(2);
% B = (hr0^2)*repmat(pa,numUsr,1).*beam./(0.5*repmat(pj,numUsr,1).*G_jk+1);
% D = zeros(numUsr-1,numUsr,N);
% for n=1:N
%     for k=1:numUsr
%      D(:,k,n) = hr0^2.*(repmat(pa(n),numUsr - 1,1))./repmat(pj(n),numUsr-1,1)./G_jw(:,k,n);
%     end
% end

B = L^2*(hr0^2)*repmat(pa,numUsr,1)./(0.5*repmat(pj,numUsr,1).*G_jk+1);
D = zeros(numUsr-1,numUsr,N);
for n=1:N
    for k=1:numUsr
     D(:,k,n) = L^2*hr0^2.*(repmat(pa(n),numUsr - 1,1))./repmat(pj(n),numUsr-1,1)./G_jw(:,k,n);
    end
end
% D = hr0^2.*repmat(pa,numUsr-1,1)./repmat(pj,numUsr-1,1)./G_jw;

Consts_Trj_UIRS = {Po, Pi, c0, c1,c2, epsilon , Ds, kf,  qr_I,  pl, Hr ,Vmax, Amax, dt, qk,qa,lambda_c, Lx, Ly, delta_x, delta_y,R2};

[~, PF_j] = flight_pow(vr,vj);
% psi_frac = 0;
lambda_r = zeros(numUsr,N,20);
% Num_rTrj = A*log(1+pa.*g_aub./(0.5*pj.*g_jb + 1));
x_lo = repmat(norms(qr-qa),numUsr,1);
y_lo = zeros(numUsr,N);
for k=1:numUsr
    y_lo(k,:) = norms(qr - qk(:,k));
end
% z_lo = zeros(1,N);
% for n= 1:N
%     trm1_lo = (qk(:,idx(n))-qr(:,n))./norms(qk(:,idx(n))-qr(:,n));
%     trm2_lo = (qr(:,n) - qa)/norms(qr(:,n) - qa);
%     w_lo = zeros(1,L);
%     for lx = 1:Lx
%         for ly = 1:Ly
%             l=(ly-1)*Lx+lx;
%             delta_lxly = [(lx-1)*delta_x; (ly-1)*delta_y; 0];
%             w_lo(l) = (2*pi/lambda_c)*(trm1_lo - trm2_lo)'*delta_lxly;
% 
%         end
%     end
%     z_lo(n) = (sum(cos(w_lo)))^2 + (sum(sin(w_lo)))^2;  
% end
Num_rTrj  = A.*log(1+B.*exp(-kf*(x_lo+y_lo))./(x_lo.^pl)./(y_lo.^pl));

% Num_rTrj = zeros(numUsr,N);
% for k=1:numUsr
% Num_rTrj(k,:) = A(k,:).*log(1+B(k,:).*(exp(-(kf/pl)*(x_lo + y_lo(k,:)))./x_lo./y_lo(k,:)).^pl);
% end

t_lo = sqrt(sqrt(1+c2^2*norms(vr).^4) - c2*norms(vr).^2);
[PF_r, ~] = flight_pow_cvxapprox(vr,vj,t_lo);

Den_rTrj = repmat((PF_r + PF_j)/(Po+Pi)/2,numUsr,1);
m=1;
lambda_r(:,:,m) =  sqrt(Num_rTrj)./Den_rTrj;

optVal_transformed(m) =  min(sum(Num_rTrj./Den_rTrj,2)/N);
optVal(m) = min(sum(2*lambda_r(:,:,m).*sqrt(Num_rTrj) - (lambda_r(:,:,m).^2).*Den_rTrj,2)/N);

%%
while (true)
[qr_new,vr_new,x_lo,y_lo,t_lo] = Trj_uavIRS_optim(x_lo,y_lo,t_lo,qr, vr, alpha, A, B, D, PF_j, qj,lambda_r(:,:,m), Consts_Trj_UIRS);
m = m+1;
qr = qr_new;
vr = vr_new;   

% x_lo = repmat(norms(qr-qa),numUsr,1);
% y_lo = zeros(numUsr,N);
% for k=1:numUsr
%     y_lo(k,:) = norms(qr - qk(:,k));
% end

Num_rTrj  = A.*log(1+B.*exp(-kf*(x_lo+y_lo))./(x_lo.^pl)./(y_lo.^pl));

[PF_r, ~] = flight_pow_cvxapprox(vr,vj,t_lo);

Den_rTrj = repmat((PF_r + PF_j)/(Po+Pi)/2,numUsr,1);

lambda_r(:,:,m) =  sqrt(Num_rTrj)./Den_rTrj;

optVal_transformed(m) =  min(sum(Num_rTrj./Den_rTrj,2)/N);
optVal(m) = min(sum(2*lambda_r(:,:,m).*sqrt(Num_rTrj) - (lambda_r(:,:,m).^2).*Den_rTrj,2)/N);

    ite_index = ite_index + 1;
    AEE_Calc

err_uIRCTrj = (optVal(m)-optVal(m-1))/optVal(m-1);
fprintf('Concave-Convex FP of sum-ratio maximization error at iteration #%d is %2.3e!\n', m-1, err_uIRCTrj);

if (abs(err_uIRCTrj) <= eps_frac_algo)
    break;
end
end

