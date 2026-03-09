Channels

% G_auw = zeros(numUsr-1,N);
% G_jw = zeros(numUsr-1,N);
% g_aub = zeros(1,N);
% for n=1:N
%     temp1 = G_auk(:,n);
%     temp2 = G_jk(:,n);
%     
%     g_aub(n) = G_auk(idx(n),n);
% 
%     temp1(idx(n)) = [];
%     G_auw(:,n) = temp1;
%     
%     temp2(idx(n)) = [];
%     G_jw(:,n) = temp2; 
% end

A = alpha.*BW/log(2);
B = repmat(pa,numUsr,1).*G_auk;
C = repmat(0.5*pj,numUsr,1);
D = zeros(numUsr-1,numUsr,N);
for n=1:N
    for k=1:numUsr
     D(:,k,n) = (repmat(pa(n),numUsr - 1,1).*G_auw(:,k,n)) ./(repmat(pj(n),numUsr - 1, 1));
    end
end
m=1;

[PF_r, ~] = flight_pow(vr,vj);

lambda_j = zeros(numUsr,N,20);

% A.*log(1+B.*L.*(exp(-(kf/pl)*(x_lo + y_lo))./x_lo./y_lo).^pl)
w_lo = zeros(numUsr,N);
for k=1:numUsr
    w_lo(k,:) = (hj0^2)*exp(-kf*norms(qj-qk(:,k)))./((norms(qj-qk(:,k))).^(pl));
end
Num_jTrj = A.*log(1+B./(C.*w_lo+1));

t_lo = sqrt(sqrt(1+c2^2*norms(vj).^4) - c2*norms(vj).^2);
[~, PF_j] = flight_pow_cvxapprox(vr,vj,t_lo);
Den_jTrj = repmat((PF_r + PF_j)/(Po+Pi)/2,numUsr,1);

lambda_j(:,:,m) =  sqrt(Num_jTrj)./Den_jTrj;

optVal_transformed_j(m) = min(sum(Num_jTrj./Den_jTrj,2)/N);
optVal_j(m) = min(sum(2*lambda_j(:,:,m).*sqrt(Num_jTrj) - (lambda_j(:,:,m).^2).*Den_jTrj,2)/N);

Consts_Trj_CJU = {Po, Pi, c0, c1,c2,  epsilon , Ds, kf,  qj_I,  pl, Hj ,Vmax, Amax, dt, qk, numUsr, N,qa, R2,pa,pj,hj0};
%%
while (true)
    [qj,vj,t_lo,w_lo] = Trj_CJU_optim(t_lo,w_lo,alpha,qj, vj, A, B, C, D, PF_r, qr, reshape(lambda_j(:,:,m),numUsr,N), Consts_Trj_CJU);
    m = m+1;  
    
%     w_lo = zeros(numUsr,N);
%     for k=1:numUsr
%         w_lo(k,:) = (hj0^2)*exp(-kf*norms(qj-qk(:,k)))./((norms(qj-qk(:,k))).^(pl));
%     end
    
    Num_jTrj =  A.*log(1+B./(C.*w_lo+1));

%     [PF_r, PF_j] = flight_pow(vr,vj);
    
    [~, PF_j] = flight_pow_cvxapprox(vr,vj,t_lo);
    Den_jTrj = repmat((PF_r + PF_j)/(Po+Pi)/2,numUsr,1);
    
    lambda_j(:,:,m) =  sqrt(Num_jTrj)./Den_jTrj;

    optVal_transformed_j(m) = min(sum(Num_jTrj./Den_jTrj,2)/N);
    optVal_j(m) = min(sum(2*lambda_j(:,:,m).*sqrt(Num_jTrj) - (lambda_j(:,:,m).^2).*Den_jTrj,2)/N);
    
    ite_index = ite_index + 1;
    AEE_Calc

    err_ucjTrj = (optVal_j(m)-optVal_j(m-1))/optVal_j(m-1);
    fprintf('Concave-Convex FP of sum-ratio maximization error at iteration #%d is %2.3e!\n', m-1, err_ucjTrj);
    
    if (abs(err_ucjTrj) <= eps_frac_algo)
        break;
    end

end



