function pj_new = JammingPow_optim(pj,A,B,C,D, epsilon,Pj_max, alpha)

[numUsr, N] = size(alpha);
pj_lo = pj;
cvx_begin quiet
cvx_solver mosek
    variable pj_opt(1,N) nonnegative
    variable rho(1,N)
    variable x(1,N) nonnegative
    
    maximize sum(rho)
    
    subject to
    
        A.*log(1+B./(pj_lo.*C+1)) ...
         - A.*B.*C.*(pj_opt-pj_lo)./(C.*pj_lo+1)./(1+B+C.*pj_lo) >= rho; 
        D.*repmat(inv_pos(pj_opt), numUsr-1, 1) <= repmat(x, numUsr-1, 1);
        x <= 1 - epsilon;
        pj_opt <= Pj_max;

cvx_end


pj_new = pj_opt';

end

