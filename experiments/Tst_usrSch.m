clc
clear all
close all

T = 30; epsilon = 0.01; N = 30; 
filename = sprintf('myUsrSch_T%deps%d.mat',T,100*epsilon);  
pName = "C:\Users\mtat0004\Google Drive\MyPhDResearch\IEEEIoT_IRSCovert_NEW";
dName = "Simulations";

% Run proposed algorithm for joint design
SystemParams
ite_index = 1;
Feasible_Init  % initilizing a feasible point

alpha = rand(numUsr,N);
print_results(ATR,APC,AEE)
% Visulization


Channels
[PF_r, PF_j] = flight_pow(vr,vj);
 
Rkn = BW*log2(1+repmat(pa,numUsr,1).*G_auk./ (0.5*repmat(pj,numUsr,1).*G_jk + 1));
Akn = Rkn./(repmat((PF_r+PF_j)/(Po+Pi)/2,numUsr,1));

Bkn = zeros(numUsr,N);
for n=1:N
    for k=1:numUsr

           Bwn = 1-repmat(pa(n),numUsr-1,1).*G_auw(:,k,n) ...
                 ./repmat(pj(n),numUsr-1,1)./G_jw(:,k,n); 

           Bkn(k,n) = min(Bwn);
    end
end

mu_MAX = 2e15;
mu = power(2,-10); itr_algo = 1; cte = 2;
alpha_old = alpha;

%%
while (true)

cvx_begin  quiet
cvx_solver mosek

    variable alpha_opt(numUsr,N) nonnegative
    variable eta nonnegative
    variable s_opt 

    maximize (s_opt - eta * mu)

    subject to 
    (sum(alpha_opt.*Akn,2)) >= s_opt;
    sum(alpha_opt.*Bkn,1) >= 1 - epsilon;  
    0 <= alpha_opt <= 1;
    sum(alpha_opt,1) <= 1;
    sum(sum( (1 - 2*alpha).*alpha_opt+ alpha.^2))  <= eta;
    
cvx_end

%%=========================================================================
alpha = full(alpha_opt);
% alpha(alpha <=0) = 0;
% alpha_new = Sch_optim_bin(alpha, Akn,zeta_min, epsilon);
fprintf('itr_algo = #%d, mu = %2.3e, eta = %2.3e\n',itr_algo, mu, eta);
 if (eta ==0)
    Channels
    itr_algo = itr_algo +1;
    ITR = alpha*BW.*log2(1+repmat(pa,numUsr,1).*G_auk./(0.5*repmat(pj,numUsr,1).*G_jk + 1));
    [PF_r, PF_j] = flight_pow(vr,vj);
    IFP= (PF_r+ PF_j)/(Po+Pi)/2; 
    AEE(itr_algo) = min(mean(ITR./repmat(IFP,numUsr,1),2));    
    break;
 % optimal alpha cannot be achieved, so we keep previous  value of alpha
 elseif (mu == mu_MAX)
    alpha = alpha_old;
    break;
 else
    Channels
    itr_algo = itr_algo +1;
    ITR = alpha*BW.*log2(1+repmat(pa,numUsr,1).*G_auk./(0.5*repmat(pj,numUsr,1).*G_jk + 1));
    [PF_r, PF_j] = flight_pow(vr,vj);
    IFP= (PF_r+ PF_j)/(Po+Pi)/2; 
    AEE(itr_algo) = min(mean(ITR./repmat(IFP,numUsr,1),2));    
    mu = min(cte*mu,mu_MAX);
 end

end

   

save(filename, 'AEE')

%%
figure
plot(1:itr_algo,AEE)