function PHI_new = IRSbeamformingNEW_optim(alpha,Akn, Bkn, Cmkn, Ek, Ea, epsilon,idx)

[numUsr, L, N] = size(Ek);
PHI_new = zeros(L,N);
itr_algo = 0;

cvx_begin  sdp quiet
cvx_solver  mosek
    variable W(L,L,N) hermitian semidefinite
    expression R(L-1,L-1)

    expressions  s_opt(numUsr,N) t_opt(numUsr-1,numUsr,N)
    
 maximize (min(sum(Akn.*log(1+Bkn.*s_opt),2)))
%  maximize (min(min(s_opt,2)))


subject to


    for n=1:N
         

        for k = 1:numUsr
            ek = reshape(Ek(k,:,n),L,1);
            ea = reshape(Ea(:,n),L,1);
            Ak = diag(ek')*ea; % L-by-1
            AAk = Ak*Ak'; % L-by-L 
            s_opt(k,n) = trace(AAk*W(:,:,n)); 
            j=1;
            for m=1:numUsr
                if (m~=k)
                    em = reshape(Ek(m,:,n),L,1);
                    Bm = diag(em')*ea; 
                    BBm = Bm*Bm';  
                    t_opt(j,k,n) = 1 - Cmkn(j,k,n)*trace(BBm*W(:,:,n));   
                    j = j + 1;
                else
                    continue;
                end
                
                
            end
        end
        
%         sum(alpha(:,n).*t_opt(:,n)) >= 1 - epsilon;
        W(:,:,n) == hermitian_semidefinite(L) ;
        diag(W(:,:,n)) <= 1;
        
    end
    
    
    for n = 1:N
      for m = 1: numUsr-1   
        t_opt(m,:,n)*alpha(:,n) >= 1 - epsilon;

      end       
    end
    
    
    
cvx_end

%%
% W_lo = W;
% V_lo = zeros(L,L-1,N);
% for n=1:N
%     [V_lo(:,:,n), ~] =  eigs(W_lo(:,:,n),L-1,'smallestabs') ;
% end
% itr_algo = itr_algo + 1;
% 
% eps_IRM = 1e-2;
% cte = 2;
% eta = 10;
% 
% while(true)
% cvx_begin  sdp quiet
% cvx_solver  mosek
%     variable W(L,L,N) hermitian semidefinite
%     variable r nonnegative
% 
%     expression R(L-1,L-1,N)
%     expressions  s_opt(numUsr,N) t_opt(numUsr-1,numUsr,N)
%     
%  maximize (min(sum(Akn.*log(1+Bkn.*s_opt),2)) - r*eta)
% %  maximize (min(min(s_opt,2)))
% 
% 
% subject to
% 
% for n=1:N
% 
%         for k = 1:numUsr
%             ek = reshape(Ek(k,:,n),L,1);
%             ea = reshape(Ea(:,n),L,1);
%             Ak = diag(ek')*ea; % L-by-1
%             AAk = Ak*Ak'; % L-by-L 
%             s_opt(k,n) = trace(AAk*W(:,:,n)); 
%             j=1;
%             for m=1:numUsr
%                 if (m~=k)
%                     em = reshape(Ek(m,:,n),L,1);
%                     Bm = diag(em')*ea; 
%                     BBm = Bm*Bm';  
%                     t_opt(j,k,n) = 1 - Cmkn(j,k,n)*trace(BBm*W(:,:,n));   
%                     j = j + 1;
%                 else
%                     continue;
%                 end
%                 
%                 
%             end
%         end
%         
% %         sum(alpha(:,n).*t_opt(:,n)) >= 1 - epsilon;
%         W(:,:,n) == hermitian_semidefinite(L) ;
%         diag(W(:,:,n)) <= 1;
%         
%         R(:,:,n) = r*eye(L-1) - V_lo(:,:,n)'*W_lo(:,:,n)*V_lo(:,:,n);
%         0.5*(R(:,:,n)+R(:,:,n)') >= 0;
% end
%     
% 
% for n = 1:N
%   for m = 1: numUsr-1   
%     t_opt(m,:,n)*alpha(:,n) >= 1 - epsilon;
% 
%   end       
% end
%     
%     
%     
% cvx_end
%         
% %%      
% if(r <= eps_IRM || itr_algo >=5)
%     break;
% else
%     fprintf('r = %f\n', r);
%     W_lo = W;
%     V_lo = zeros(L,L-1,N);
%     for n=1:N
%         [V_lo(:,:,n), ~] =  eigs(W_lo(:,:,n),L-1,'smallestabs') ;
%     end
%     itr_algo = itr_algo + 1;
%     eta = eta*cte;
% 
% end
% 
% end
%%
for n = 1: N
    [U,S,V] = svd(W(:,:,n));
    u = U(:,1);
    s = S(1,1);
    v = V(:,1);
    w = u*sqrt(s);
    PHI_new(:,idx(n),n) = w;
end

end
