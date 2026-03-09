function [pa_new, alpha_new] =  UsrSchPow_optim(pa, pj, alpha, Pa_tot, ...
                               Pa_max, Akn, G_auk,G_jk, epsilon)

[numUsr, N] = size(alpha);

alpha_lo = alpha;


G_auw = zeros(numUsr-1,N);
G_jw = zeros(numUsr-1,N);
s_lo = zeros(numUsr,N);
t_lo = log(1+repmat(pa,numUsr,1).*Akn) ;

for k = 1: numUsr
   for n=1:N
       temp1 = G_auk(:,n);
       temp2 = G_jk(:,n);
       temp1(k) = [];
       G_auw(:,n) = temp1;
       temp2(k) = [];
       G_jw(:,n) = temp2; 
       Bwn = G_auw(:,n)./pj(n)./G_jw(:,n); 
       s_lo(k, n) = min(1 - pa(n) .* Bwn);
   end         
end
     
tst = zeros(2,20);

%% 
% itr_algo_MAX = 1.5e1;
% mu = 0.01; itr_algo = 0; cte = 5;
% while (itr_algo < itr_algo_MAX)   
cvx_begin  quiet
cvx_solver mosek

    variable pa_opt(1,N) nonnegative
    variable alpha_opt(numUsr,N) nonnegative
 
    variable t(numUsr,N) nonnegative
    variable s(numUsr,N) nonnegative
        
    variable myPhi  
     variable eta

    maximize myPhi

    subject to 

%% --- constraints
   
     sum(-(alpha_lo + t_lo).^2 + ...
          2*(alpha_lo.* alpha_opt +  ...
               t_lo.*t) - square_pos(alpha_opt - t), 2) >= 4*myPhi;
          
    
        
     log(1+repmat(pa_opt, numUsr, 1).*Akn) >= t;
     
      
     
     for k = 1: numUsr
           for n=1:N
               temp1 = G_auk(:,n);
               temp2 = G_jk(:,n);
%                g_aub = G_auk(k,n);
%                g_jb =  G_jk(k,n);
               temp1(k) = [];
               G_auw(:,n) = temp1;
               temp2(k) = [];
               G_jw(:,n) = temp2; 
               Bwn = G_auw(:,n)./pj(n)./G_jw(:,n); 
               for m=1:numUsr-1
                    1 - pa_opt(n) .* Bwn >= s(k,n);    
               end
           end         
     end

     
      sum(-(alpha_lo + s_lo).^2 + ...
            2*(alpha_lo.* alpha_opt +  ...
               s_lo.*s) - ...
                  square_pos(alpha_opt - s), 1) >= 4*(1 - epsilon);
     
    sum(pa_opt) <= Pa_tot;
    0 <= pa_opt <= Pa_max;
    
    0 <= alpha_opt <= 1;
    sum(alpha_opt,1) <= 1;
    
    sum(sum(alpha_opt + alpha_lo.^2 - 2*alpha_lo.*alpha_opt))  <= 0;
   

  
    
    
cvx_end


% alpha_lo = full(alpha_opt);
% t_lo = t;
% s_lo = s;
% fprintf('itr_algo = #%d, eta = %2.3e\n',itr_algo, eta);
% itr_algo = itr_algo +1;
% mu = cte*mu;
%   tst(1,itr_algo) = myPhi;
%    tst(2,itr_algo) = eta;
% 
% 
% end
% 
% figure
% plot(1:itr_algo,tst(1,1:itr_algo))
% hold on
% plot(1:itr_algo,tst(2,1:itr_algo))
% legend('penalty parameter')


alpha_new = full(alpha_opt)';
pa_new = pa_opt';
end



%% feasiblity check
% log(1+repmat(pa, numUsr, 1).*Akn) >= t_lo
% 
% 
%  for k = 1: numUsr
%        for n=1:N
%            temp1 = G_auk(:,n);
%            temp2 = G_jk(:,n);
%            temp1(k) = [];
%            G_auw(:,n) = temp1;
%            temp2(k) = [];
%            G_jw(:,n) = temp2; 
%            Bwn = G_auw./pj./G_jw; 
%            for m=1:numUsr-1
%                 1 - pa(n) .* Bwn(m,n) >= s_lo(k,n) 
%            end
%        end         
%  end
% 
%   sum(-(alpha_lo + s_lo).^2 + ...
%         2*(alpha_lo.* alpha_lo +  ...
%            s_lo.*s_lo) - ...
%               (alpha_lo - s_lo).^2, 1) >= 4*(1 - epsilon);
% sum(pa) <= Pa_tot
% 0 <= pa <= Pa_max
% 
% 0 <= alpha <= 1
% sum(alpha_lo,1) <= 1
% 
% sum(sum(alpha_lo + alpha_lo.^2 - 2*alpha_lo.*alpha_lo))  <= 0
