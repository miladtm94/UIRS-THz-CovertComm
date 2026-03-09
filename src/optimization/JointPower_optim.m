function [pa_new, pj_new] = JointPower_optim(pa,pj,A,B,C,D, epsilon,Pj_max,Pa_max,P_tot, alpha)

[numUsr, N] = size(alpha);
rho_lo = pa./pj;
p_lo = pj;
% t_lo = 1./(pj.*C+pa.*B);
t_lo = sqrt(B.*pa);
% s_lo = sqrt(B.*pa);
% total constraints: 5N+1
% total variables: 2N
cvx_begin quiet
cvx_solver mosek
    variable rho_opt(1,N) nonnegative
    variable p_opt(1,N) nonnegative
    variable t(1,N) nonnegative
    variable s(1,N) nonnegative

    expressions g_lb(1,N)
    
    % cvx - ccv
%     f_lb = log(1+1./t_lo) - (t-t_lo)./(t_lo.*(1+t_lo));
     g_lb = log(1+p_lo.*C) + (C./(1+p_lo.*C)).*(p_opt - p_lo);


%     maximize (sum(A.*(log(1+t+C.*p_opt) - g_lb)))
    
    maximize (sum(A.*(s - g_lb)))
    subject to
    
%     4*t./B <= (-(p_lo+rho_lo).^2+2*(p_lo + rho_lo).*(p_opt+rho_opt) ...
%                   -square_pos(p_opt-rho_opt));
             
    for n= 1:N
        B(n)*p_opt(n) >= quad_over_lin(t(n),rho_opt(n));   
    end
    
        1-t_lo.^2+2*t_lo.*t + C.*p_opt >= exp(s);
%     1+C.*p_opt <= (s_lo.^2./t_lo) + (2*s_lo./t_lo).*(s-s_lo) + ...
%                                            (-s_lo.^2./(t_lo.^2)).*(t-t_lo);

%         log(p_opt) + log(t) + log(C+p_opt.*B) >= 0;
        
        D.*repmat(rho_opt,numUsr-1,1) <= epsilon; % N
        
        sum(p_opt+0.25*(square_pos(p_opt+rho_opt) -  ...
         (-(p_lo-rho_lo).^2 + 2.*(p_lo-rho_lo).*(p_opt-rho_opt)))) <= P_tot; % 1
        
        p_opt <= Pj_max; % N
        
        % 
        square_pos(p_opt+rho_opt) - ...
        (-(p_lo-rho_lo).^2 + 2.*(p_lo-rho_lo).*(p_opt-rho_opt)) <= 4*min(Pa_max, min(epsilon*Pj_max./D,[],1)); % N
        

cvx_end
pa_new = (rho_opt.*p_opt)';
pj_new = p_opt';

end
