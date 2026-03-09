function alpha_new = Sch_optim_bin(alpha, Akn,zeta_min, epsilon)


[numUsr, N] = size(alpha);

 
cvx_begin  
cvx_solver mosek

    variable alpha_opt(numUsr,N) binary
%     variable t(1,N) nonnegative
%     variable eta nonnegative
    %variable s_opt(numUsr,N) nonnegative

    maximize (min(sum(alpha_opt.*Akn,2)))

    subject to 

%% --- constraints        
          
%     alpha_opt.*Akn >= repmat(t,numUsr,1); 
     
    sum(alpha_opt.*zeta_min,1) >= 1 - epsilon;  
    sum(alpha_opt,1) <= 1;
      
cvx_end

alpha_new = full(alpha_opt)';
end