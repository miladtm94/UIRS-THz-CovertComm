% theta = zeros(L,N);
% for lx = 1:Lx
%     for ly = 1:Ly
%         l=(ly-1)*Lx+lx;
%         delta_lxly = [(lx-1)*delta_x; (ly-1)*delta_y; 0];
%         temp = (2*pi*(qr-qa)'*delta_lxly)./lambda_c./reshape(norms(qr-qa), N, 1);
%          theta(l, :) = temp';
%     end
% end
% 
% beta = zeros(L, N);
% for n=1:N
%     for lx = 1:Lx
%         for ly = 1:Ly
%             l=(ly-1)*Lx+lx;
%             delta_lxly = [(lx-1)*delta_x; (ly-1)*delta_y; 0];
%             temp = (2*pi*(qk(:,idx(n))-qr(:,n))'*delta_lxly)./lambda_c./norms(qk(:,idx(n))-qr(:,n));
%             beta(l, n) = temp;
%         end
%     end
%         
% end
% 
% Phi = exp(-1i*(theta - beta));


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

% theta = zeros(L,N);
% for lx = 1:Lx
%     for ly = 1:Ly
%         l=(ly-1)*Lx+lx;
%         delta_lxly = [(lx-1)*delta_x; (ly-1)*delta_y; 0];
%         temp = (2*pi*(qr-qa)'*delta_lxly)./lambda_c./reshape(norms(qr-qa), N, 1);
%          theta(l, :) = temp';
%     end
% end
% 
% Beta = zeros(numUsr,N,L);
% for k = 1: numUsr
%     for n=1:N
%         for lx = 1:Lx
%             for ly = 1:Ly
%                 l=(ly-1)*Lx+lx;
%                 delta_lxly = [(lx-1)*delta_x; (ly-1)*delta_y; 0];
%                 temp = (2*pi*(qk(:,k)-qr(:,n))'*delta_lxly)./lambda_c./norms(qk(:,k)-qr(:,n));
%                 Beta(k, n, l) = temp;
%             end
%         end
% 
%     end
% end
% 
% for k=1:numUsr
%      tmp = reshape(Beta(k, :,:), N,L);
%      beta = tmp';
%      Phi = exp(-1i*(theta - beta));
%      Channels
%      
%      Rkn = BW*log2(1+repmat(pa,numUsr,1).*G_auk./ (0.5*repmat(pj,numUsr,1).*G_jk + 1));
%      Akn = Rkn./(repmat((PF_r+PF_j)/(Po+Pi)/2,numUsr,1));
%      
%      for n=1:N
%        G_auw = G_auk(:,n);
%        G_jw = G_jk(:,n);
%        G_auw(k) = [];
%        G_jw(k) = [];
%        Bwn = 1-repmat(pa(n),numUsr-1,1).*G_auw ...
%              ./repmat(pj(n),numUsr-1,1)./G_jw; 
% 
%        Bkn(k,n) = min(Bwn);
%     end         
% end


%%
mu_MAX = 2e15;
mu = power(2,-10); itr_algo = 0; cte = 2;
alpha_old = alpha;
while (true)

[alpha, eta] = Sch_optim(alpha, Akn,Bkn, epsilon, mu);
alpha(alpha <=0) = 0;

% alpha_new = Sch_optim_bin(alpha, Akn,zeta_min, epsilon);
fprintf('itr_algo = #%d, mu = %2.3e, eta = %2.3e\n',itr_algo, mu, eta);
 if (eta ==0)
    break;
 % optimal alpha cannot be achieved, so we keep previous  value of alpha
 elseif (mu == mu_MAX)
    alpha = alpha_old;
    break;
 else
    itr_algo = itr_algo +1;
    mu = min(cte*mu,mu_MAX);
 end

end
ite_index = ite_index + 1;
    
for n=1:N
   [~, idx(n)] = max(alpha(:,n));
end
