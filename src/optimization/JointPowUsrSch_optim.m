function [alpha_new, pa_new, pj_new] = JointPowUsrSch_optim(alpha, pa,pj, A, G_auk,G_jk, epsilon,Pj_max,Pa_max,P_tot)

[numUsr, N] = size(alpha);

w_lo = log(1+G_auk.*repmat(pa,numUsr,1)./(1+0.5*G_jk.*repmat(pj,numUsr,1)));

s_lo = zeros(numUsr,N);
for k=1:numUsr

    G_auw = zeros(numUsr-1,N);
    G_jw = zeros(numUsr-1,N);
    g_aub = zeros(1,N);
    g_jb = zeros(1,N);

    for n=1:N
    temp1 = G_auk(:,n);
    temp2 = G_jk(:,n);

    g_aub(n) = G_auk(k,n);
    g_jb(n) =  G_jk(k,n);

    temp1(k) = [];
    G_auw(:,n) = temp1;

    temp2(k) = [];
    G_jw(:,n) = temp2; 
    end
    
    
    s_lo(k,:) = min(1-(repmat(pa,numUsr-1,1).*G_auw)./(repmat(pj,numUsr-1,1).*G_jw),[],1);
end
%%
mu_MAX = 2e10;
mu = 2e-5; itr_algo = 0; cte = 2;
while (true)   
cvx_begin  
cvx_solver mosek

    variable alpha_opt(numUsr,N) nonnegative
    variable pa_opt(1,N) nonnegative
    variable pj_opt(1,N) nonnegative 
    variable w(numUsr,N) nonnegative
    variable s(numUsr,N) nonnegative
    variable r(1,N) nonnegative

    variable eta nonnegative
    expressions g_lb(numUsr,N) objFunc

    % min(A,[],2) column vector and min(A,[],1) row vector
    objFunc = min(sum(repmat(0.25*A,numUsr,1).*(-(alpha+w_lo).^2 + 2*(alpha+w_lo).*(alpha_opt+w) - square_pos(alpha_opt-w)),2)) ;
    %max sum (min_k alpha_k[n] Aln(1+B_k,n pa[n]/(C_k,n pj[n]+1)))
    maximize (objFunc - eta * mu)

    subject to 

%% --- constraints         
    
    sum(pa_opt+pj_opt) <= P_tot; % 1
    0 <= pj_opt <= Pj_max; % N
    0 <= pa_opt <= Pa_max; % N 


for k=1:numUsr

    B = G_auk(k,:);
    C = G_jk(k,:)/2;
    
    g_lb(k,:) = log(1+C.*pj) + (C./(1+pj.*C)).*(pj_opt - pj);
    log(1+B.*pa_opt+C.*pj_opt) - g_lb(k,:) >= w(k,:);  

    for n=1:N
    temp1 = G_auk(:,n);
    temp2 = G_jk(:,n);

    temp1(k) = [];
    G_auw = temp1;

    temp2(k) = [];
    G_jw = temp2; 
    
    D = G_auw./G_jw;
    for m=1:numUsr-1
        1 -(D(m)*quad_over_lin(pa_opt(n),r(n))) >= s(k,n);
    end
    pa_opt(n).*G_auw <= pj_opt(n).*G_jw;
    
    end 
    
end

    4*r <= (-(pa+pj).^2 + 2*(pa+pj).*(pa_opt+pj_opt) - square_pos(pa_opt-pj_opt));
    sum(( 2*(alpha+s_lo).*(alpha_opt-alpha+s-s_lo)+(alpha+s_lo).^2 - square_pos(alpha_opt-s)),1) >= 4*(1-epsilon);
    
%      sum((-(alpha+s_lo).^2 + 2*(alpha+s_lo).*(alpha+s_lo) - (alpha-s_lo).^2),1)
%      sum((-(alpha+s_lo).^2 + 2*(alpha+s_lo).*(alpha_opt+s) - (alpha_opt-s).^2),1) 
%      sum(s.*alpha_opt,1)
%      
     
    0 <= alpha_opt <= 1;
    sum(alpha_opt,1) <= 1;
    
    % (a - a^2) <= 0  OR  a>=1 && a<=0
    sum(sum((alpha_opt + alpha.^2 - 2*alpha.*alpha_opt)))  <= eta;
    
cvx_end

%%=========================================================================
alpha = full(alpha_opt);
pa = pa_opt;
pj = pj_opt;
% w_lo = w;
% s_lo = s;
fprintf('itr_algo = #%d, eta = %2.3e, mu = %2.3e, obj = %2.3e\n', itr_algo, eta, mu, objFunc);
 if (eta <= 1e-6 || mu == mu_MAX)
    break;
 else
    itr_algo = itr_algo +1;
    mu = min(cte*mu,mu_MAX);
end

end
alpha_new = full(alpha_opt)';
pa_new = pa_opt';
pj_new = pj_opt';
end