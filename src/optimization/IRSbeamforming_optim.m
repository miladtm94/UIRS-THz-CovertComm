function PHI_new = IRSbeamforming_optim(C, Ek, Ea, idx, epsilon,n)

[numUsr, L, ~] = size(Ek);
eta_MAX = power(2,20);
eps_IRM = 1e-2;
itr_algo = 0;
cte = 2;
eta = power(2,-10);
% eta_MAX = power(2,20);
% cte = 2;
% 
% eta = power(2,-10); itr_algo = 0;
%   vl = cos(Phi(:,n))-1i*sin(Phi(:,n));
%     v_lo = rand(L,1) + 1i*rand(L,1);

%%%%=====  Initialize a positive semidefinite complex matrix 
%     matrixSize = L;
%     while true
%       AA = rand(matrixSize, matrixSize)+1i*rand(matrixSize, matrixSize);
%       if rank(AA) == matrixSize; break; end    %will be true nearly all the time
%     end
%     BB = AA' * AA;
%     W_lo = BB./max(max(abs(BB)));
    
%     rank(W_lo);
%     try  chol(W_lo);
%         disp('Matrix is symmetric positive definite.')
%     catch ME
%         disp('Matrix is not symmetric positive definite')
%     end


eb = reshape(Ek(idx(n),:,n),L,1);
ea = reshape(Ea(:,n),L,1);
Ab = diag(eb')*ea; % L-by-1
AAb = Ab*Ab'; % L-by-L

ek = reshape(Ek(:,:,n),numUsr,L);
em = zeros(numUsr-1, L);
for l=1:L
  tmp1 = ek(:,l);
  tmp1(idx(n))=[];
  em(:,l) = tmp1;
end

AAm = zeros(L,L,numUsr-1);
for m=1:numUsr-1
   Am = diag(em(m,:)')*ea;
   AAm(:,:,m) = Am*Am';
end
     
cvx_begin  sdp quiet
cvx_solver  mosek
    %cvx_precision('low')
    variable W(L,L) complex semidefinite

    maximize trace(AAb*W)

    subject to

        % 1 - Cm Tr(Bm * W) >= w, forall m 
        for m=1:numUsr-1
            C(m,n)*trace(AAm(:,:,m)*W) <= epsilon;
        end

        % Wl,l <= 1 for all l=1 .. L
        diag(W) <= 1;

        % rank(W)==1    Tr(W) - eigval_max(W) <= 0  OR [wmax_lo, ~] = eigs(W,1);
        %trace(W) - trace(wmax_lo*wmax_lo'*W) <= mYmu;
        W >= 0;
        %W-1e-12*eye(L) >= 0; % LMI  Tr(W) - eigval_max(W) >= 0
        % min(eig(W)) >= 0
        % W == hermitian_semidefinite(L);
cvx_end

% try  chol(W);
%         fprintf('Y ')
%     catch ME
%         fprintf('N ')
% end
%     
%     [V_lo, ~] =  eigs(W,L-1,'smallestabs') ;
%     itr_algo = itr_algo + 1;
%     W_lo = W;
%     [wmax_lo, ~] = eigs(W_lo,1);
%    trace(W_lo) - trace((wmax_lo*wmax_lo')*W_lo);

%     fprintf("Initial: rank(W) = %d\n", rank(W_lo));
    

%      for m=1:numUsr-1
%            Am = diag(em(m,:)')*ea;
%            AAm = Am*Am';
%            1-C(m,n)*trace(AAm*W_lo) >= 1-epsilon;
%      end
        
% while (eta  < eta_MAX)
%%
% while(true)
% cvx_begin sdp quiet
% cvx_solver sedumi
%     %cvx_precision('low')
%     variable omega 
%     variable r nonnegative
%     variable W(L,L) hermitian semidefinite
%     expression R(L-1,L-1)
%     maximize (trace(AAb*W)- eta*r)
%     
%     subject to
%         
%         omega >= 1- epsilon;
%         
%         % 1 - Cm Tr(Bm * W) >= w, forall m 
%         for m=1:numUsr-1
%            1-C(m,n)*trace(AAm(:,:,m)*W) >= omega;
%         end
%               
%         % Wl,l <= 1 for all l=1 .. L
%         diag(W) <= 1;
%           
%         % rank(W)==1    Tr(W) - eigval_max(W) <= 0  OR [wmax_lo, ~] = eigs(W,1);
% %         trace(W) - trace(wmax_lo*wmax_lo'*W) <= mYmu;
%         
%         R = r*eye(L-1) - V_lo'*W_lo*V_lo;
%         0.5*(R+R') >= 0;
%         
%         W-1e-15*eye(L) >= 0; % LMI  Tr(W) - eigval_max(W) >= 0
%         % min(eig(W)) >= 0
%         W == hermitian_semidefinite(L);
%         
%        
% cvx_end
% 
% [V_lo, ~] =  eigs(W,L-1,'smallestabs') ;
% W_lo = W;
% % [wmax_lo, ~] = eigs(W_lo,1);
% itr_algo = itr_algo +1;
% eta = min(eta*cte,eta_MAX);
% [U,S,V] = svd(W);
% u = U(:,1);
% s = S(1,1);
% v = V(:,1);
% w = u*sqrt(s);
% W_app = w*w';   
% fprintf('Iteration: #%d, eta = %2.3e,  r = %2.3e, err = %2.2e, Rank(W)=%d\n',itr_algo, eta, r, norm(W-W_app), rank(W));
% 
% if(r <= eps_IRM)
%     break;
% end
% % if(rank(W)==1)
% %     fprintf('Mission accomplished for time slot #%d :) Yaay!\n', n);
% %     break;
% % elseif(rank(W)~=1 && r <= eps_IRM)
% %     printf('Termination criteria is reached but searching rank-one matrix W failed :(\n');
% % elseif (eta == eta_MAX)
% %     printf('Termination criteria is reached but searching rank-one matrix W failed :(\n');
% %     break;
% % end
% 
% end
%%


    [U,S,V] = svd(W);
    u = U(:,1);
    s = S(1,1);
    v = V(:,1);
    w = u*sqrt(s);
    PHI_new = w;
end
% try  chol(W);
%         fprintf('Y ')
%     catch ME
%         fprintf('N ')
% end
   %     if(rank(W)==1)
%         fprintf('Time slot %d: W is PSD and Rank-one. Yaay!\n', n);
%         break;       
%     end  

    
    % = v* v^H


% assert(norm(W_app-W)<1e-10)
%%
   
%     w = u*sqrt(s);
%     W = w*w';
%     W = W_lo;
%     [U,S,V] = svd(W);
%     u = U(:,1);
%     s = S(1,1);
%     v = V(:,1);

%     [wmax_lo, ~] = eigs(W,1);
%     break;

%     rank(W);
 
%     fprintf('Time slot: #%d, rho = %2.3e,  mYeps = %2.3e\n',n,rho, mYeps);
%     PHI_new(:,n) = w;

% figure
% semilogy(1:itr_algo,tst(1,1:itr_algo))
% hold on
% semilogy(1:itr_algo,tst(2,1:itr_algo))
% hold on
% semilogy(1:itr_algo,tst(3,1:itr_algo))
% legend('function value with penalty', 'optimizaiton value','penalty parameter')



