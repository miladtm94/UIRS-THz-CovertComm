x_lo = norms(qr-qa);
y_lo = norms(qr - qk(:,idx));

beta_sel = zeros(L,N);
for n=1:N
    beta_sel(:,n) = reshape(beta(idx(n),:,n),L,1);
end
 
% exp(-kf(x+y))/(x^pl/2 * y^pl/2) >= w
f0 = exp(-kf*((x_lo) + (y_lo)))./(x_lo.^(pl).*y_lo.^(pl));
f1x = f0.*(pl+kf*(x_lo))./x_lo;
f1y = f0.*(pl+kf*(y_lo))./y_lo;

ww_lo = f0;

z_lo = zeros(1,N);
for n= 1:N
    trm1_lo = (qk(:,idx(n))-qr(:,n))./norms(qk(:,idx(n))-qr(:,n));
    trm2_lo = (qr(:,n) - qa)/norms(qr(:,n) - qa);
    w_lo = zeros(1,L);
    for lx = 1:Lx
        for ly = 1:Ly
            l=(ly-1)*Lx+lx;
            delta_lxly = [(lx-1)*delta_x; (ly-1)*delta_y; 0];
            w_lo(l) = (2*pi/lambda_c)*(trm1_lo - trm2_lo)'*delta_lxly;

        end
    end
    z_lo(n) = (sum(cos(w_lo)))^2 + (sum(sin(w_lo)))^2;  
end


tt_lo = sqrt(min(f0,[],1).*z_lo)-1e-10;
%% =================== Objective function ===========================
%  sum(min_k BW*log2(1+B*f(qr)))./(g(vr) +Pfj);


% f0 - f1x.*(x-x_lo) -  f1y.*(y-y_lo) >= ww_lo;

% f(z) is an approx linear finction of "qr" corresponding to the IRS phase
for n= 1:N   
    trm1_lo = (qk(:,idx(n))-qr(:,n))./norms(qk(:,idx(n))-qr(:,n));
    trm2_lo = (qr(:,n) - qa)/norms(qr(:,n) - qa);

    w_lo = zeros(1,L);
    for lx = 1:Lx
        for ly = 1:Ly
            l=(ly-1)*Lx+lx;
            delta_lxly = [(lx-1)*delta_x; (ly-1)*delta_y; 0];
            w_lo(l) = (2*pi/lambda_c)*(trm1_lo - trm2_lo)'*delta_lxly;
        end
    end
%     Phi_opt(:,n) + theta(:,n) - beta_sel(:,n) == w';

    z_lo(n) = (sum(cos(w_lo)))^2 + (sum(sin(w_lo)))^2;
end

% A*ln(1+B*z*ww) >= u   ---->  ww*z >= tt^2   tt^2 >= (exp(u/A)-1)/B
%   tt_lo.^2 >= (exp(u/A)-1)./B;
 (tt_lo.^2./z_lo) <= ww_lo
   (tt_lo.^2./z_lo) - ww_lo


% % ||qr - qa|| <= x
% norms(repmat(qa,1,N) - qr) <= x_lo
% 
% % ||qr - qk|| <= y
% norms(qk(:,idx) - qr) <= y_lo

%% ==================== Safety flying Distance ============================
% || qr - qj|| >= Ds;
norms(qj-qr).^2  >= Ds^2

%% ===================== Covertness constraint ============================
 % 1- min_w D * exp(-kf(||qr - qa|| + ||qr - qw||))/(||qr - qa|| * ||qr -
 % qw||)^(pl) < 1- epsilon; for all w <IN> K/{selected user k}
%  t1 >= 1;
%  t2 >= 1;
t1_lo = norms(qa-qr);
t2_lo = zeros(numUsr-1, N);
yy_lo = zeros(numUsr-1, N);
for n=1:N
    qw = qk;
    qw(:,idx(n))=[];
    for m=1:numUsr-1
        t2_lo(m, n) = norms(qw(:,m)-qr(:,n));
        yy_lo(m,:) = log(t1_lo) + log(t2_lo(m,:)) ;

    end
end

for m = 1:numUsr-1
    % t1 * t2 >= exp(yy)  exp(-kf(t1+t2))/(t1 * t2)^(pl/2) <= r
    % exp(-kf(t1+t2))/(t1 * t2)^(pl) <= exp(-kf(t1+t2))/exp(yy*pl) <= r
    yy_lo(m,:) = log(t1_lo) + log(t2_lo(m,:)) ;
    m
    D(m,:).*(L^2).* exp(-(kf*(t1_lo+t2_lo(m,:)) + (pl)*yy_lo(m,:))) <= epsilon
%   log(D(m,:).*(L^2))+ -(kf*(t1+t2(m,:)) + (pl)*yy(m,:)) <= log(r);
end

%% =================== flight power non-convex term constraint ============
t_lo.^2 + c2.*(-norms(vr).^2 + ...
                                  2*sum(vr.*vr,1)) >= 1./(t_lo.^(-2))

%% =================== Mobility constraints ===============================
    qr(:,1) == qr_I
    qr(:,N) == qr_I
    qr(1:2, 2:N) - qr(1:2,1:N-1) == vr(:,1:N-1)*dt
    norms(vr) <= Vmax
    norms(vr(:,2:N) - vr(:,1:N-1)) <= Amax
      
%% ================= Flying zone ============================================
norms(qr-qa) <= R2
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    