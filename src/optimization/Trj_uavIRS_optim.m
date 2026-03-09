function [qr_new,vr_new,x_lo,y_lo,t_lo] = Trj_uavIRS_optim(x_lo,y_lo,t_lo,qr, vr, alpha, A, B, D, PF_j, qj,lambda_rm, Consts_Trj_UIRS)
       
[Po, Pi, c0, c1,c2,  epsilon , Ds, kf,  qr_I,  pl, Hr ,Vmax, Amax, dt, qk,qa,lambda_c, Lx, Ly, delta_x, delta_y,R2] = Consts_Trj_UIRS{:};
[numUsr, N] = size(alpha);

lambda_rm = reshape(lambda_rm,numUsr,N);
% lambda_rm = reshape(lambda_r(:,:,m),numUsr,N);

% x_lo = repmat(norms(qr-qa),numUsr,1);

% y_lo = zeros(numUsr,N);
% % v_lo = zeros(numUsr,N);
% for k=1:numUsr
%     y_lo(k,:) = norms(qr - qk(:,k));
% end

% v_lo= (exp(-(kf/pl)*(x_lo + y_lo))./x_lo./y_lo).^pl-1e-10;

% t_lo = sqrt(sqrt(1+c2^2*norms(vr).^4) - c2*norms(vr).^2);


  
% beta_sel = zeros(L,N);
% for n=1:N
%     beta_sel(:,n) = reshape(beta(idx(n),:,n),L,1);
% end
 

% exp(-kf(x+y))/(x^pl * y^pl) >= 1/v
% f0 = exp(-kf*((x_lo) + (y_lo)))./((x_lo.^(pl)).*(y_lo.^(pl)));
% f1x = f0.*(pl+kf*(x_lo))./x_lo;
% f1y = f0.*(pl+kf*(y_lo))./y_lo;

% ww_lo = f0;

% bim_lo = zeros(1,N);
% for n= 1:N
%     for k = 1:numUsr
%     trm1_lo = (qk(:,k)-qr(:,n))./norms(qk(:,k)-qr(:,n));
%     trm2_lo = (qr(:,n) - qa)/norms(qr(:,n) - qa);
%     w_lo = zeros(1,L);
%     for lx = 1:Lx
%         for ly = 1:Ly
%             l=(ly-1)*Lx+lx;
%             delta_lxly = [(lx-1)*delta_x; (ly-1)*delta_y; 0];
%             w_lo(l) = wrapTo2Pi((2*pi/lambda_c)*(trm1_lo - trm2_lo)'*delta_lxly);
%         end
%     end
%     
% %     tmp = 0;
% %     for i=1:L
% %         for j=1:L
% %             if (i~=j)
% %                 tmp = tmp + real(exp(1i*(w_lo(i) - w_lo(j))));
% %             end
% %         
% %         end
% %     end
%     z_lo(n) = (sum(cos(w_lo)))^2 + (sum(sin(w_lo)))^2;  
%     end
% end

% s_lo = -log(B.*z_lo.*f0);
% r_lo = sqrt(B.*z_lo.*exp(-kf*((x_lo) + (y_lo)))./((x_lo.^(pl)).*(y_lo.^(pl))));


% v_lo = (exp(-(kf/pl)*(x_lo + y_lo))./x_lo./y_lo).^pl;
% tt_lo = sqrt(min(f0,[],1).*z_lo)-1e-10;
% 
% zzzz_lo = 1./z_lo;
% yyyy_lo = log(exp(-kf*(x_lo+y_lo))./((x_lo.*y_lo).^pl));
%%
cvx_begin quiet
cvx_solver mosek
% cvx_precision best
variable t(1,N) nonnegative
% variable u(1,N) nonnegative
% variable tt(1,N) nonnegative
% variable ww(1,N) nonnegative
variable x(1,N) nonnegative
variable y(numUsr,N) nonnegative
% variable r(1,N) nonnegative
variable z(numUsr,N) 
variable t2(numUsr-1,numUsr,N) nonnegative
variable t1(1,N)  nonnegative
variables qr_opt(2,N) vr_opt(2,N)
% variable u(numUsr-1,N) 
% variable r(1,N) nonnegative
variable v(numUsr,N) nonnegative
% variable u(1,N) nonnegative
variable minACR
% variable zzzz_opt nonnegative
% variable yyyy_opt 

% A.*log(1+B.*L.*(exp(-(kf/pl)*(x_lo + y_lo))./x_lo./y_lo).^pl)

% Num_jTrj = A*log(1+B./(C.*(exp(-kf*norms(qj-qk(:,idx)))./((norms(qj-qk(:,idx))).^(pl)))+1));

expressions Num(numUsr,N) Den(numUsr,N)  obj(numUsr,N) 


% f(z) is an approx linear finction of "qr" corresponding to the IRS phase
% for n= 1:N
%     trm1 = (qk(:,idx(n))-[qr_opt(:,n);Hr])./norms(qk(:,idx(n))-qr(:,n));
%     trm2 = ([qr_opt(:,n);Hr] - qa)./norms(qr(:,n) - qa);
%     
%     trm1_lo = (qk(:,idx(n))-qr(:,n))./norms(qk(:,idx(n))-qr(:,n));
%     trm2_lo = (qr(:,n) - qa)/norms(qr(:,n) - qa);
% 
%     w_lo = zeros(1,L);
%     for lx = 1:Lx
%         for ly = 1:Ly
%             l=(ly-1)*Lx+lx;
%             delta_lxly = [(lx-1)*delta_x; (ly-1)*delta_y; 0];
%             w(n,l) = (2*pi/lambda_c).*(trm1 - trm2)'*delta_lxly;
%             w_lo(l) = (2*pi/lambda_c)*(trm1_lo - trm2_lo)'*delta_lxly;
%         end
%     end
% %     Phi_opt(:,n) + theta(:,n) - beta_sel(:,n) == w';
% 
%     z(n) = (sum(cos(w_lo)))^2 + (sum(sin(w_lo)))^2 ...
%     -2*sum(sum(cos(w_lo)).*sin(w_lo).*(w(n,:)-w_lo)) ...
%     +2*sum(sum(sin(w_lo)).*cos(w_lo).*(w(n,:)-w_lo)); 
% end

Den = repmat((Po*(1+c0*square_pos(norms(vr_opt))) + c1*pow_pos(norms(vr_opt), 3) + Pi*t + PF_j)/(Po+Pi)/2,numUsr,1);

% objFunc = sum(2*lambda_rm.*sqrt(Num) - (lambda_rm.^2).*Den);

% Num = -A.*rel_entr(1,1+B.*s_slack);

f0  = A.*log(1+B.*exp(-kf*(x_lo+y_lo))./(x_lo.^pl)./(y_lo.^pl));
f1x = -(A.*((B.*pl.*exp(-kf.*(x_lo + y_lo)))./(x_lo.^(pl + 1).*y_lo.^pl)  ...
+ (B.*kf.*exp(-kf.*(x_lo + y_lo)))./(x_lo.^pl.*y_lo.^pl)))./((B.*exp(-kf.* ...
 (x_lo + y_lo)))./(x_lo.^pl.*y_lo.^pl) + 1);
f1y = -(A.*((B.*pl.*exp(-kf.*(x_lo + y_lo)))./(x_lo.^pl.*y_lo.^(pl + 1)) ...
+ (B.*kf.*exp(-kf.*(x_lo + y_lo)))./(x_lo.^pl.*y_lo.^pl)))...
./((B.*exp(-kf.*(x_lo + y_lo)))./(x_lo.^pl.*y_lo.^pl) + 1);

Num = f0  + f1x.*(repmat(x,numUsr,1)-x_lo)+f1y.*(y-y_lo);


% Num =  A.*log(B.*v+1);
% Num = u;
% g0 = log(1+exp(-s_lo)) + s_lo.*exp(-s_lo)./(exp(-s_lo)+1)

% Num = A.*( - (exp(-s_lo).*(s-s_lo))./(exp(-s_lo)+1));
% Num = A.*(log(v+B.*z) + log(u));

maximize minACR

subject to

for n=1:N
    for k=1:numUsr
        obj(k,n) = 2.*lambda_rm(k,n).*sqrt(Num(k,n)) - (lambda_rm(k,n).^2).*Den(k,n);
    end
end
    sum(obj,2)/N >= minACR;

%% =================== Objective function ===========================
%  sum(A ln(1+B*f(qr))./(g(vr) +pfj));
% log(1+B.*z_lo.*tt)>= s;

%  A.*log(1+B.*z_lo.*tt) >= s;

% u >= inv_pos(v);
% u+s <= log((B.*z_lo));
% t^2 >= s
% 1+ (-r_lo.^2 + 2*r_lo.*r) >= exp(u);
% Num = u;
%  1./(x_lo.*y_lo) - (x - x_lo)./(x_lo.^2.*y_lo) - (y - y_lo)./(x_lo.*y_lo.^2) ...
%      >= exp((kf/pl).*(x+y) -s);
% 
%  B.*z_lo.*(1./u_lo -(u-u_lo)./(u_lo.^2)) >= exp(pl*s);

% exp(-kf/pl (x+y))/x/y >= v^(1/pl) ===> cvx >= ccv
% for k=1:numUsr
% -(exp(-(kf*(x_lo + y_lo(k,:)))/pl).*(pl*x.*y_lo(k,:) - kf*x_lo.^2.*y_lo(k,:) - ...
%   kf*x_lo.*y_lo(k,:).^2 + pl*x_lo.*y(k,:) - 3*pl.*x_lo.*y_lo(k,:) + kf*x.*x_lo.*y_lo(k,:) + ...
%       kf*x_lo.*y(k,:).*y_lo(k,:)))./(pl*x_lo.^2.*y_lo(k,:).^2) >= (exp(log(v_lo(k,:))/pl)  ...
%                                     + exp(log(v_lo(k,:))/pl).*(v(k,:)-v_lo(k,:))/pl./v_lo(k,:));
% end


% -(-rel_entr(1,x) -rel_entr(1,y) - (1/pl)*rel_entr(1,tt)) >= (kf/pl)*(x+y);

% ww >= 0;



% A*ln(1+B*z*ww) >= u   ---->  ww*z >= tt^2   tt^2 >= (exp(u/A)-1)/B
%    -A.*rel_entr(1,B.*(-tt_lo.^2+2*tt_lo.*tt)+1) >= u;
% A.*log(B.*s_slack+1) >= u;

% for n=1:N
%     quad_over_lin(r(n),z(n)) <= B(n).*v(n);
% end


% f0 - f1x.*(x-x_lo) -  f1y.*(y-y_lo) >= inv_pos(v);
% 
% 2*log(r) + s >=0 ;
%  for n=1:N
%     quad_over_lin(tt(n),ww(n)) <= z(n);
%     ww(n)>=0 ;
%  end



% Num = A*log(1 + B.*exp(yyyy_lo)./zzzz_lo) +  ...
%     A.*B.*((exp(yyyy_lo).*(yyyy_opt - yyyy_lo))./(zzzz_lo + exp(yyyy_lo)) ...
%        - (exp(yyyy_lo).*(zzzz_opt - zzzz_lo))./(zzzz_lo.*(zzzz_lo + exp(yyyy_lo))));
% 
% z >= inv_pos(zzzz_opt);

% ||qr - qa|| <= x
norms(repmat(qa,1,N) - [qr_opt;repmat(Hr,1,N)]) <= x;
% -norms(qr - qa).^2 + 2*sum((qr-qa).*([qr_opt;repmat(Hr,1,N)]-repmat(qa,1,N)),1) >= square_pos(x);

% ||qr - qk|| <= y
for k=1:numUsr
    norms(repmat(qk(:,k),1,N) - [qr_opt;repmat(Hr,1,N)]) <= y(k,:);
end
% -norms(qr - qk(:,idx)).^2 + 2*sum((qr-qk(:,idx)).*([qr_opt;repmat(Hr,1,N)]-qk(:,idx)),1) >= square_pos(y);

% (x_lo.*y_lo - pl.*x.*y_lo - pl.*x_lo.*y + 2*pl.*x_lo.*y_lo) ...
%                ./(x_lo.*y_lo.*(x_lo.*y_lo).^pl) >= exp(yyyy_opt + kf* (x+y));

%% ==================== Safety flying Distance ============================
% || qr - qj|| >= Ds;
 -norms(qj-qr).^2 + 2*sum((qj-qr).*(qj - [qr_opt;repmat(Hr,1,N)]),1) >= Ds^2;

%% ===================== Covertness constraint ============================
 % min_w (1-  D * exp(-kf(||qr - qa|| + ||qr - qw||))/(||qr - qa|| * ||qr -
 % qw||)^(pl)) >= 1- epsilon; for all w <IN> K/{selected user k}
%  t1 >= 1;
%  t2 >= 1;
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
%     
    % t1 * t2 >= exp(yy)  exp(-kf(t1+t2))/(t1 * t2)^(pl) <= r
    % exp(-kf(t1+t2))/(t1 * t2)^(pl) <= exp(-kf(t1+t2))/exp(yy*pl) <= r
% for m=1:numUsr-1 
% %     log(t1) + log(t2(m,:)) + log((epsilon./D(m,:)).^(1/pl)) +(kf/pl)*(t1+t2(m,:)) >= 0;
%     
%     log(z(m,:))+pl*(log(t1)+log(t2(m,:)))+kf.*(t1+t2(m,:)) >= 0;
%     D(m,:).*z(m,:) <= epsilon;
%     z(m,:) >= 0;
% %     log(D(m,:)./epsilon)-kf.*(t1+t2(m,:)) <= pl*(log(t1)+log(t2(m,:)))
%      %     -rel_entr(1,t1)-rel_entr(1,t2(m,:))>= yy(m,:);
%     
% %      {-(kf*(t1+t2(m,:)) + (pl)*yy(m,:)),1,zzz(m,:)} == exponential(1);
% 
% %      D(m,:).*(L^2).* zzz(m,:) <= r;
% 
% %       D(m,:).*(L^2).* exp(-(kf*(t1+t2(m,:)) + (pl)*yy(m,:))) <= r;
% %   log(D(m,:).*(L^2))+ -(kf*(t1+t2(m,:)) + (pl)*yy(m,:)) <= log(r);
% end

% for m = 1:numUsr-1
%     % t1 * t2 >= exp(yy)  exp(-kf(t1+t2))/(t1 * t2)^(pl/2) <= r
%     % exp(-kf(t1+t2))/(t1 * t2)^(pl/2) <= exp(-kf(t1+t2))/exp(yy*pl/2) <= r
%    t1 .* t2(m,:) >= exp(yy(m,:));
%    D(m,:).*(L^2).* exp(-(kf*(t1+t2(m,:)) + (pl/2)*yy(m,:))) <= r;
% end

 
%     r <= epsilon;

%  t1_lo = norms(qa-qr);

    % ||qa-qr|| >= t1
    -norms(repmat(qa,1,N)-qr).^2 + 2*sum((repmat(qa,1,N)-qr).*(repmat(qa,1,N) - [qr_opt;repmat(Hr,1,N)]),1) >= square_pos(t1);
for n=1:N
    % || qr-qw || >= t2
    for k=1:numUsr
        qw = qk;
        qw(:,k)=[];
        for m=1:numUsr-1
            -norms(qw(:,m)-qr(:,n)).^2 + 2*(qw(:,m)-qr(:,n))'*(qw(:,m) - [qr_opt(:,n);Hr]) >= square_pos(t2(m,k, n));       
            log((1-z(k,n))/D(m,k,n))+pl*(log(t1(n))+log(t2(m,k, n)))+kf*(t1(n)+t2(m,k, n)) >= 0;      
        end
    end
end

sum(alpha.*z,1) >= 1 - epsilon;
z<=1;

%% =================== flight power non-convex term constraint ============
-t_lo.^2 +2*t.*t_lo + 2*c2.*(-norms(vr).^2 + ...
                                  2*sum(vr.*vr_opt,1)) >= pow_p(t,-2);

%% =================== Mobility constraints ===============================
    qr_opt(:,1) == qr_I(1:2);
    qr_opt(:,N) == qr_I(1:2);
    qr_opt(:, 2:N)== vr_opt(:,1:N-1)*dt + qr_opt(:,1:N-1);
    norms(vr_opt) <= Vmax;
    norms(vr_opt(:,2:N) - vr_opt(:,1:N-1)) <= Amax; 
    
    
%% ================= Flying zone ============================================
norms(qr_opt-repmat([qa(1);qa(2)],1,N)) <= R2;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
cvx_end
%%
t_lo = full(t);
y_lo = full(y);
x_lo = full(repmat(x,numUsr,1));
qr_new = [qr_opt; repmat(Hr,1,N)];
vr_new = vr_opt;
% Phi_new = exp(1i*Phi_opt)';
end



%%
% 
% % f0 - f1x.*(x-x_lo) -  f1y.*(y-y_lo) >= ww;
% 
% for n= 1:N
%     trm1_lo = (qk(:,idx(n))-qr(:,n))./norms(qk(:,idx(n))-qr(:,n));
%     trm2_lo = (qr(:,n) - qa)./norms(qr(:,n) - qa);
% 
% w_lo = zeros(1,L);
%     for lx = 1:Lx
%         for ly = 1:Ly
%             l=(ly-1)*Lx+lx;
%             delta_lxly = [(lx-1)*delta_x; (ly-1)*delta_y; 0];
%             w_lo(l) = (2*pi/lambda_c)*(trm1_lo - trm2_lo)'*delta_lxly;
% 
%         end
%     end
% 
% z_lo(n) = (sum(cos(w_lo)))^2 + (sum(sin(w_lo)))^2;  
% end
% 
% 
% for n=1:N
%     -norms(qj(:,n)-qr(:,n)).^2 + 2*(qj(:,n)-qr(:,n))'*(qj(:,n) - qr(:,n)) - Ds^2
% end
% 
% r >= 1 - epsilon;
% 
% for m = 1:numUsr-1
% % t1 * t2 >= exp(y)  exp(-kf(t1+t2))/(t1 * t2)^pl/2 <=t
% log(t1) + log(t2(m,:)) >= yy(m,:);
% 1-D(m,:).*(L^2).*exp(-(kf*(t1+t2(m,:)) + (pl/2)*yy(m,:))) >= r;  
% end
% 
% % ||qa-qr|| >= t1
% % || qa-qw || >= t2
% 
% for n=1:N
% -norms(qa-qr(:,n)).^2 + 2*(qa-qr(:,n))'*(qa - qr(:,n)) >= square_pos(t1(n));
% qw = qk;
% qw(:,idx(n))=[];
% for m=1:numUsr-1
% -norms(qw(:,m)-qr(:,n)).^2 + 2*(qw(:,m)-qr(:,n))'*(qw(:,m) - qr_opt(:,n)) >= square_pos(t2(m, n));
% end
% end
% 
% 
% for n=1:N
% -t_lo(n).^2 +2*t(n).*t_lo(n) + c2.*(-norms(vr(:,n)).^2 + ...
%                       2*vr(:,n)'*vr_opt(:,n)) >= pow_p(t(n),-2);
% end
% 
% 
% qr(:,1) == qr_I
% qr(:,N) == qr_I
% qr(3,:) == H
% vr(3,:) == 0
% qr(:, 2:N) - qr(:,1:N-1) == vr(:,1:N-1)*dt
% norms(vr) <= Vmax
% norms(vr(:,2:N) - vr(:,1:N-1)) <= Amax 

