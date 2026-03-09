function [alpha_new, eta_new] = Sch_optim(alpha, Akn,Bkn,epsilon, mu)


[numUsr, N] = size(alpha);

% tmp = reshape(alpha, [1, numUsr*N]);
% tmp(tmp==1) = 1-epsilon;
% s_lo = reshape(tmp, [numUsr, N]);

%%
   
cvx_begin  quiet
cvx_solver mosek

    variable alpha_opt(numUsr,N) nonnegative
%   variable t(1,N) nonnegative
    variable eta nonnegative
    %variable s_opt(numUsr,N) nonnegative

    maximize (min((sum(alpha_opt.*Akn,2))) - eta * mu)

    subject to 

%% --- constraints        
%     alpha_opt.*Akn >= repmat(t,numUsr,1); 
     
    sum(alpha_opt.*Bkn,1) >= 1 - epsilon;
%     sum (alpha_opt.*s_lo,1) >= 1 - epsilon ;   
    % alpha*s =(1/4)*((alpha+s)^2-(alpha-s)^2))
%     sum(-(alpha + s_lo).^2 + ...
%             2*(alpha_lo.* alpha_opt +  ...
%                s_lo.*s_opt) - square_pos(alpha_opt - s_opt), 1) >= 4*(1 - epsilon);
    
    0 <= alpha_opt <= 1;
    sum(alpha_opt,1) <= 1;
    
    % (a - a^2) <= 0  OR  a>=1 && a<=0
    sum(sum( (1 - 2*alpha).*alpha_opt+ alpha.^2))  <= eta;
    
cvx_end

%%=========================================================================
alpha_new = full(alpha_opt);
eta_new = eta;
end