function [pa_new, pj_new] = JointPowerWoRecast_optim(pa,pj,A,B,C,D, epsilon,Pj_max,Pa_max,P_tot, alpha)

[numUsr, N] = size(alpha);


cvx_begin quiet
cvx_solver mosek
    variable pa_opt(1,N) nonnegative
    variable pj_opt(1,N) nonnegative
    variable s(numUsr,N) nonnegative
    variable eta
%     expressions g_lb(1,N)
    
    % cvx - ccv
     g_lb = log(1+C.*repmat(pj,numUsr,1)) + (C./(1+repmat(pj,numUsr,1).*C)).*(repmat(pj_opt - pj,numUsr,1));

     maximize eta
     
    
    subject to
    
        sum(A.*(log(1+B.*repmat(pa_opt,numUsr,1)+C.*repmat(pj_opt,numUsr,1)) - g_lb),2) >= eta;
    
        0 <= pa_opt <= Pa_max;
        sum(pa_opt+pj_opt) <= P_tot; % 1
        0 <= pj_opt <= Pj_max; % N

        sum(alpha.*s,1) >= 1 - epsilon;
        s >= 0;
        
        
        % sum (alpha * min_m(1- D pa/pj)) >= 1-epsilon
   for n=1:N               
        for k=1:numUsr
%             s_lo = min(1 - D(:,k,n).*pa(n)./pj(n)); 
            for m=1:numUsr-1
                log(pj_opt(n)) + log(1-s(k,n)) >= log(D(m,k,n)*pa(n)) ...
                    + (pa_opt(n)-pa(n))/(pa(n));
%                (pj(n) - s_lo + 1) ...
%                *(2*pj_opt(n) - pj(n) - 2*s(k,n) + s_lo + 1) ...
%                -square_pos(1-s(k,n)-pj_opt(n)) >= 4*D(m,k,n)*pa_opt(n); 
            end
        end
   end
        
cvx_end
pa_new = pa_opt;
pj_new = pj_opt;

end
