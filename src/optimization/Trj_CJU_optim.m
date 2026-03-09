function [qj_new,vj_new,t_lo,w_lo] = Trj_CJU_optim(t_lo,w_lo,alpha,qj, vj, A, B, C,D,  PF_r, qr, lambda_jm, Consts_Trj_CJU)

[Po, Pi, c0, c1,c2,  epsilon , Ds, kf,  qj_I,  pl, Hj ,Vmax, Amax, dt, qk, numUsr, N, qa, R2,pa,pj,hj0] = Consts_Trj_CJU{:};
lambda_jm = reshape(lambda_jm,numUsr,N);
% lambda_jm = reshape(lambda_j(:,:,m),numUsr,N);
% lambda_jm = lambda_j(:,:,m);
% w_lo = zeros(numUsr,N);
% for k=1:numUsr
%     w_lo(k,:) = (hj0^2)*exp(-kf*norms(qj-qk(:,k)))./((norms(qj-qk(:,k))).^(pl));
% end
% t_lo = sqrt(sqrt(1+c2^2*norms(vj).^4) - c2*norms(vj).^2);

% s_lo = zeros(numUsr-1,numUsr,N);

% f =  log(1+a exp(-bx) / (x^c))
% df = - a(c+bx)/x/(a+x^c exp(bx))

% f0 = log(1+B+C.*exp(-kf.*xx_lo./(xx_lo.^pl)));
% bb = C./(1+B);
% f1 = -bb.*(pl+kf.*xx_lo)./xx_lo./(bb+xx_lo.^pl.*exp(kf.*xx_lo));
 
f0 = A.*log(1+B./(C.*w_lo+1));
f1 = -(A.*B.*C)./((C.*w_lo+1).*(B+C.*w_lo+1)); 
    
cvx_begin quiet
cvx_solver mosek
    variable qj_opt(2,N) 
    variable vj_opt(2,N)
    variable v(numUsr,N) nonnegative
    variable w(numUsr,N) nonnegative
%     variables u1(numUsr,N) u2(numUsr,N)
%     variable xx(numUsr,N) nonnegative
%     variable yy(numUsr,N) nonnegative
    variable t(1,N) nonnegative
    variable minACR
    variable s(numUsr-1,numUsr,N) nonnegative
    variable r(numUsr-1,numUsr,N) nonnegative
    variable z(numUsr,N) nonnegative
%     variable u(N,numUsr,numUsr-1) nonnegative
    variable y(numUsr-1,numUsr,N) nonnegative

    expressions Num(numUsr,N) Den(numUsr,N)  obj(numUsr,N)

%==== Objective function  ==========================================%%%%%
   
%     norms([qj_opt;repmat(Hj,1,N)]-qk(:,idx)) <= v;
%     {rr,v,w_lo.^(-pl) - pl*w_lo.^(-(pl+1)).*(w-w_lo)} == exponential(1);
%     (kf/pl)*square_pos(v) <= rr;
        
    Num = f0 + f1.*(w - w_lo) ;
    
    Den = repmat((Po*(1+c0*square_pos(norms(vj_opt))) + c1*pow_pos(norms(vj_opt), 3) + Pi*t + PF_r)/(Po+Pi)/2,numUsr,1);

    % 2*lambda_m sqrt(f_num) - 
    
    maximize minACR
%     maximize 0

    subject to
    for n=1:N
        for k=1:numUsr
            obj(k,n) = 2.*lambda_jm(k,n).*sqrt(Num(k,n)) - (lambda_jm(k,n).^2).*Den(k,n);
        end
    end
    sum(obj,2)/N >= minACR;
    
%     u1 + u2 <= 1;
%     for k=1:numUsr
%         for n=1:N
%            {v(k,n)-w(k,n),1,u1(k,n)} == exponential(1);
%            {-w(k,n),1,u2(k,n)} == exponential(1);
%         end
%     end
%     pl.*log(yy) + kf*yy + v -log(C) >=0;
       
    %  v^pl = exp(pl log(v))    exp(-kf v)/(v^pl) <= w  v>=0  w>=0
    log(w)+ pl*log(v) +kf*v - 2*log(hj0)  >= 0;

    % v <= ||qj - qk||
    for k= 1:numUsr
        for n=1:N
            square_pos(v(k,n)) <= -norms(qj(:,n)-qk(:,k)).^2 + ...
                   2*sum((qj(:,n)-qk(:,k)).*([qj_opt(:,n);Hj]-qk(:,k)));
%             xx(k,n) >= norms([qj_opt(:,n);Hj]-qk(:,k));
        end
    end
%%  t^ + c2||vj||^2 >= 1/t^2        ||qj - qr|| >= Ds

%% =================== flight power non-convex term constraint ============
% for n=1:N
% -t_lo(n).^2 +2*t(n).*t_lo(n) + 2*c2.*(-norms(vj(:,n)).^2 + ...
%                          2*vj(:,n)'*vj_opt(:,n)) >= pow_p(t(n),-2);
% end

-t_lo.^2 +2*t.*t_lo + 2*c2.*(-norms(vj).^2 + ...
                                  2*sum(vj.*vj_opt,1)) >= pow_p(t,-2);
% for n=1:N
%        -norms(qj(:,n)-qr(:,n)).^2 + ...
%            2*((qj(:,n)-qr(:,n))'*([qj_opt(:,n);Hj] - qr(:,n))) >= Ds^2;
% end

%% ==== Mobility constraints ==========================================%%%%
    qj_opt(:,1) == qj_I(1:2);
    qj_opt(:,N) == qj_I(1:2);
    qj_opt(:, 2:N) ==  vj_opt(:,1:N-1)*dt + qj_opt(:,1:N-1);
    norms(vj_opt)' <= Vmax;
    norms(vj_opt(:,2:N) - vj_opt(:,1:N-1))' <= Amax; 
    

    
    
%%======= Covert constraints =====================================%%%%%%%%
      
% D ||qj - qw||^pl *exp(kf * ||qj - qw||) <= epsilon
% 
% ||qj - qw||^pl <= s  -->  ||qj - qw|| <= s^(1/pl)
% D s^pl exp(kf s) <= epsilon
% s exp((kf*s^(1+1/pl))/s) <= z


for n=1:N
    for k=1:numUsr
       
       qw = qk;
       qw(:,k)=[];
%        for m=1:numUsr-1 
%             norms([qj_opt(:,n);Hj]-qw(:,m)) <= s(n,k,m); 
%             {r(n,k,m),s(n,k,m),y(n,k,m)} == exponential(1);
%             (kf/pl)*square_pos(s(n,k,m)) <= r(n,k,m);
%             1- pow_p((D(m,k,n)^(1/pl))*y(n,k,m),pl) >= z(k,n);      
%        end 

       for m=1:numUsr-1 
%            s_lo(m,k,n) = norms(qj(:,n)-qw(:,m));
            norms([qj_opt(:,n);Hj]-qw(:,m)) <= s(m,k,n); 
        
%             log(1-z(k,n)) - log(D(m,k,n)/(hj0^2)) +kf*s(m,k,n) ...
%              - pl*(log(s_lo(m,k,n))+ (s(m,k,n) - s_lo(m,k,n))/(s_lo(m,k,n))) >= 0;
%          
            sqrt(y(m,k,n)) >= sqrt((kf/pl))*(s(m,k,n));
%             y > 0 , y*exp(x/y) <= z
%                 as
%              {x,y,z} == exponential(1)
            {y(m,k,n),s(m,k,n),r(m,k,n)} == exponential(1);
            1- pow_p(((D(m,k,n)/(hj0^2))^(1/pl))*r(m,k,n),pl) >= z(k,n);  
            s(m,k,n) >= 0; y(m,k,n)>=0;
       end 
       
    end  
    
end
sum(alpha.*z,1) >= 1 - epsilon;
z<=1;

%% ================= Flying zone ============================================
norms(qj_opt-repmat([qa(1);qa(2)],1,N)) <= R2;

%% ==================== Safety flying Distance ============================
% || qr - qj|| >= Ds;
% for n=1:N
%     -norms(qj(:,n)-qr(:,n)).^2 + 2*(qj(:,n)-qr(:,n))'*([qj_opt(:,n);Hj] - qr(:,n)) >= Ds^2;
% end
 -norms(qr-qj).^2 + 2*sum((qr-qj).*(qr - [qj_opt;repmat(Hj,1,N)]),1) >= Ds^2;

%%=========================================================================
cvx_end
%%
    w_lo = full(w);
    t_lo = full(t);
    qj_new = [qj_opt; repmat(Hj,1,N)];
    vj_new = vj_opt;
end

