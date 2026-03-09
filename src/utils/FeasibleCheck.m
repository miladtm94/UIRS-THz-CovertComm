%% Feasibly check
zeta = zeros(numUsr,N);
for n=1:N
    for k=1:numUsr

           Bwn = 1-repmat(pa(n),numUsr-1,1).*G_auw(:,k,n) ...
                 ./repmat(pj(n),numUsr-1,1)./G_jw(:,k,n); 

           zeta(k,n) = min(Bwn);
    end
end
C1 = sum(alpha.*zeta,1) >= 1-epsilon
sum(alpha.*zeta,1) - (1-epsilon)
C2 = [qj(:,1)'== qj(:,N)', qr(:,1)'== qr(:,N)', norms(qr-qj) >= Ds]
C3 =[norms(vr) <= Vmax, norms(vj) <= Vmax, norms((vr(:,2:N) - vr(:, 1:N-1))) <= Amax , ...
norms((vj(:,2:N) - vj(:, 1:N-1))) <= Amax]

C5 = [qr(1:2,2:N) == (qr(1:2,1:N-1) + vr(:,1:N-1)*dt), ...
qj(1:2,2:N) == qj(1:2,1:N-1) + vj(:,1:N-1)*dt]

C7 = [pa<=Pa_max; pj<=Pj_max; repmat(sum(pa+pj) <= P_tot,1,N)]
%% subproblem 1 feasiblity check - User schedilng

alpha_lo = alpha;
zeta_min = zeros(numUsr,N);
for k = 1: numUsr
   for n=1:N
       temp1 = G_auk(:,n);
       temp2 = G_jk(:,n);
       temp1(k) = [];
       G_auw(:,n) = temp1;
       temp2(k) = [];
       G_jw(:,n) = temp2; 

       Bwn = 1-repmat(pa(n),numUsr-1,1).*G_auw(:,n) ...
             ./repmat(pj(n),numUsr-1,1)./G_jw(:,n); 

       zeta_min(k,n) = min(Bwn);
    end         
end
sum(alpha.*zeta_min,1) >= 1 - epsilon
0 <= alpha <= 1
sum(alpha,1) <= 1
sum(sum((alpha + alpha.^2 - 2*alpha_lo.*alpha))) <= 0

%%
% Subproblem 3

Channels
G_auw = zeros(numUsr-1,N);
G_jw = zeros(numUsr-1,N);
g_aub = zeros(1,N);
g_jb = zeros(1,N);
for n=1:N
    temp1 = G_auk(:,n);
    temp2 = G_jk(:,n);
   
    g_aub(n) = G_auk(idx(n),n);
    g_jb(n) =  G_jk(idx(n),n);

    temp1(idx(n)) = [];
    G_auw(:,n) = temp1;
    
    temp2(idx(n)) = [];
    G_jw(:,n) = temp2; 
end

[PF_r, PF_j] = flight_pow(vr,vj);

A = BW./log(2)./((PF_r + PF_j)/(Po+Pi));
B = g_aub;
C = g_jb/2;
D= G_auw./G_jw;

[numUsr, N] = size(alpha);
rho_lo = pa./pj;
p_lo = pj;
t_lo = B.*pa;

C1 =4*t_lo/B <= (-(p_lo+rho_lo).^2+2*(p_lo + rho_lo).*(p_lo + rho_lo) ...
    -(p_lo-rho_lo).^2)

4*t_lo/B - (-(p_lo+rho_lo).^2+2*(p_lo + rho_lo).*(p_lo + rho_lo) ...
    -(p_lo-rho_lo).^2)

C2 = D.*repmat(rho_lo,numUsr-1,1) <= epsilon % N

C4 = sum(p_lo+0.25*((p_lo+rho_lo).^2 -  ...
(-(p_lo-rho_lo).^2 + 2.*(p_lo-rho_lo).^2))) <= P_tot % 1

C5 = p_lo <= Pj_max % N
C6 = (p_lo+rho_lo).^2 - ...
(-(p_lo-rho_lo).^2 + 2.*(p_lo-rho_lo).^2) <= 4*Pa_max % N



        
        % sum (alpha * min_m(1- D pa/pj)) >= 1-epsilon
   s_lo = zeros(N,numUsr);     
   for n=1:N               
        for k=1:numUsr
            for m=1:numUsr-1
                s_lo(n,k) = min(1-D(:,k,n)*pa(n)/pj(n));
            end
        end
   end

    for n=1:N               
        for k=1:numUsr
            for m=1:numUsr-1
                log(pj(n)) + log(1-s(k,n)) >= log(D(m,k,n)*pa(n)) ...
                    + (pa(n)-pa(n))/(pa(n));
            end
        end
   end


        
        sum(alpha.*s_lo',1) >= 1 - epsilon;
        0 <= s_lo <= 1;





%%
% Subproblem 4    UIRS trajectory
C1 =  -(exp(-(kf*(x_lo + y_lo))/pl).*(pl*x_lo.*y_lo - kf*x_lo.^2.*y_lo - ...
  kf*x_lo.*y_lo.^2 + pl*x_lo.*y_lo - 3*pl.*x_lo.*y_lo + kf*x_lo.*x_lo.*y_lo + ...
      kf*x_lo.*y_lo.*y_lo))./(pl*x_lo.^2.*y_lo.^2) >= (exp(log(v_lo)/pl)  ...
                                    + exp(log(v_lo)/pl).*(v_lo-v_lo)/pl./v_lo)

-(exp(-(kf*(x_lo + y_lo))/pl).*(pl*x_lo.*y_lo - kf*x_lo.^2.*y_lo - ...
  kf*x_lo.*y_lo.^2 + pl*x_lo.*y_lo - 3*pl.*x_lo.*y_lo + kf*x_lo.*x_lo.*y_lo + ...
   kf*x_lo.*y_lo.*y_lo))./(pl*x_lo.^2.*y_lo.^2) - (exp(log(v_lo)/pl)  ...
                                    + exp(log(v_lo)/pl).*(v_lo-v_lo)/pl./v_lo)


C2 =  -norms(qj-qr).^2 + 2*sum((qj-qr).*(qj - qr),1) >= Ds^2
%% 
% 1./(x_lo.*y_lo) >= exp((kf/pl).*(x_lo+y_lo) -s_lo/pl)./((B.*z_lo).^(1/pl))
% 
% 1./(x_lo.*y_lo) -( exp((kf/pl).*(x_lo+y_lo) -s_lo/pl)./((B.*z_lo).^(1/pl)))


t1_lo = norms(qa-qr);
for n=1:N
    % || qa-qw || >= t2
    qw = qk;
    qw(:,idx(n))=[];
    for m=1:numUsr-1
        t2_lo(m,n) = norms(qw(:,m)-qr(:,n));
    end
end

% zw_lo = zeros(numUsr-1,N);
% for n=1:N
%     qw = qk;
%     qw(:,idx(n))=[];
%     for m = 1:numUsr-1
%         
%         for n= 1:N
%             trm1_lo = (qw(:,m)-qr(:,n))./norms(qw(:,m)-qr(:,n));
%             trm2_lo = (qr(:,n) - qa)/norms(qr(:,n) - qa);
%             w_lo = zeros(1,L);
%             for lx = 1:Lx
%                 for ly = 1:Ly
%                     l=(ly-1)*Lx+lx;
%                     delta_lxly = [(lx-1)*delta_x; (ly-1)*delta_y; 0];
%                     w_lo(l) = (2*pi/lambda_c)*(trm1_lo - trm2_lo)'*delta_lxly;
%                 end
%             end
%             zw_lo(m,n) = (sum(cos(w_lo)))^2 + (sum(sin(w_lo)))^2;  
%         end
%     end
% end 
for m = 1:numUsr-1
   C3 =log(t1_lo) + log(t2_lo(m,:)) + log((epsilon./D(m,:))).^(1/pl) +(kf/pl)*(t1_lo+t2_lo(m,:)) >= 0
    log(t1_lo) + log(t2_lo(m,:)) + log((epsilon./D(m,:)).^(1/pl)) +(kf/pl)*(t1_lo+t2_lo(m,:))
    log(D(m,:)./epsilon)-kf.*(t1_lo+t2_lo(m,:)) <= pl*(log(t1_lo)+log(t2_lo(m,:)))

end

%% 

C4 = t_lo.^2 + c2.*norms(vr).^2 >= t_lo.^(-2)
t_lo.^2 + 2*c2.*norms(vr).^2 - t_lo.^(-2)
% 
%% 
C5 = qr(:,1) == qr_I(:)
C6 = qr(:,N) == qr_I(:)
C7 = qr(1:2, 2:N) - qr(1:2,1:N-1) == vr(1:2,1:N-1)*dt
C8 = norms(vr) <= Vmax;
C9 = norms(vr(:,2:N) - vr(:,1:N-1)) <= Amax 
C10 = norms(qr(1:2,:)-repmat([qa(1);qa(2)],1,N)) <= R2


%%
% Subproblem 5 - UCJ traj feasiblity check
for k= 1:numUsr
    for n=1:N
        v_lo(k,n) = norms(qj(:,n)-qk(:,k)); 
    end
end
    
C1 = log(w_lo)+ pl*log(v_lo) +kf*v_lo  >= 0
C2 = t_lo.^2 + 2*c2.*norms(vj).^2 >= 1./(t_lo.^2)
t_lo.^2 + 2*c2.*norms(vj).^2 - 1./(t_lo.^2)

C3 = qj(:,1) == qj_I(:)
C4 = qj(:,N) == qj_I(:)
C5 = qj(1:2, 2:N) - qj(1:2,1:N-1) == vj(:,1:N-1)*dt
C6 = norms(vj) <= Vmax
C7 = norms(vj(:,2:N) - vj(:,1:N-1)) <= Amax
C9 = norms(qj(1:2,:)-repmat([qa(1);qa(2)],1,N)) <= R2
C10 = norms(qj-qr)>= Ds

%%
% for  n=1:N
%    qw = qk;
%    qw(:,idx(n))=[];
%   for m=1:numUsr-1 
%         s_lo(m,n) = norms(qj(:,n)-qw(:,m)); 
%   end
% end
% C8 = D.*s_lo.^pl.*exp(kf*s_lo) <= epsilon

s= zeros(N, numUsr,numUsr-1);
u= zeros(N, numUsr,numUsr-1);
z = zeros(N, numUsr);

for n=1:N
    for k=1:numUsr
       
       qw = qk;
       qw(:,k)=[];
       for m=1:numUsr-1 
             s(n,k,m)=  norms(qj(:,n)-qw(:,m)); 
             u(n,k,m)= (s(n,k,m)*exp((kf/pl)*s(n,k,m)))^pl;    
       end 
       
       z(n,k) = min(1 - D(:,k,n).*reshape(u(n,k,:),numUsr-1,1));
    end  
    
end
 
sum(alpha.*z',1) >= 1- epsilon
z <= 1

%log(s) + log(z) >= exp(-kf*s/pl);
% repmat(u,numUsr-1,1) + log(D./(1-epsilon)) >=0;
% ||qj - qw|| >= s









%%
%Joint user sh=cheduling and power allocation 
min(sum(repmat(0.25*A,numUsr,1).*(alpha.*w_lo),2)) 


    sum(pa+pj) <= P_tot % 1
    mean(0 <= pj <= Pj_max) % N
    mean(0 <= pa <= Pa_max) % N 
    
    g_lb = zeros(numUsr,N);
    
    cond_1 = zeros(numUsr,N);
    cond_2 = zeros(numUsr,N);
    r = pa.*pj;

for k=1:numUsr

    B = G_auk(k,:);
    C = G_jk(k,:)/2;
    
    g_lb(k,:) = log(1+C.*pj) + (C./(1+pj.*C)).*(pj - pj);
    cond_0(k) = mean(log(1+B.*pa+C.*pj) - g_lb(k,:) >= w_lo(k,:));
    log(1+B.*pa+C.*pj) - g_lb(k,:) - w_lo(k,:);
    
    G_auw = zeros(numUsr-1,N);
    G_jw = zeros(numUsr-1,N);
    g_aub = zeros(1,N);
    g_jb = zeros(1,N);

    for n=1:N
    temp1 = G_auk(:,n);
    temp2 = G_jk(:,n);

    g_aub(n) = G_auk(k,n);
    g_jb(n) =  G_jk(k,n);

    temp1(k) = [];
    G_auw(:,n) = temp1;

    temp2(k) = [];
    G_jw(:,n) = temp2; 
    
    D = G_auw(:,n)./G_jw(:,n);
    cond_1(k,n) = min(1-(D.*pa(n)^2/r(n)))>= s_lo(k,n);
    cond_2(k,n) = min(pa(n).*G_auw(:,n) <= pj(n).*G_jw(:,n));
    min(pa(n).*G_auw(:,n) - pj(n).*G_jw(:,n))
    end

end
    cond_0
    cond_1
    min(1-(D.*pa.^2./r(n)))- s_lo

    cond_2
    sum(alpha.*s_lo,1) >= 1-epsilon
    
    0 <= alpha <= 1
    sum(alpha,1) <= 1
    
    % (a - a^2) <= 0  OR  a>=1 && a<=0
    sum(sum((alpha - alpha.^2)))  
