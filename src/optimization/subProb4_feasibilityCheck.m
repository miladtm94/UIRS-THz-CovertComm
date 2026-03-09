SPr = sum(Po*(1+c0*(norms(vr)).^2) + c1*power(norms(vr), 3) + Pi*t_lo)/N/100 + Pfj0;

objFunc = eta_lo - psi_frac * SPr;

%% =================== Objective function ===========================
% min_k mean(BW*log2(1+B*f(qr))))/(g(vr) +Pfj0);

x_lo = norms(repmat(qa,1,N) - qr);
 
t_lo = (200 - (Po*(1+c0*(norms(vr).^2))...
      + c1*(norms(vr)).^3))/ Pi;

SPr = sum(Po*(1+c0*(norms(vr)).^2) + c1*(norms(vr).^3) + Pi*t_lo)/N/100 + Pfj0;


y_lo = norms(qk(:,idx) - qr);
% exp(-kf(x+y))/(x^pl/2 * y^pl/2) >= w
f0 = exp(-kf*((x_lo) + (y_lo)))./(x_lo.^(pl/2).*y_lo.^(pl/2));

A*ones(1,N)*u_lo' >= eta;

ww_lo = f0;

% f(z) is an approx linear finction of "qr" corresponding to the IRS phase
for n= 1:N
    trm1_lo = (qk(:,idx(n))-qr(:,n))./norms(qk(:,idx(n))-qr(:,n));
    trm2_lo = (qr(:,n) - qa)./norms(qr(:,n) - qa);

    w_lo = zeros(1,L);
    for lx = 1:Lx
        for ly = 1:Ly
            l=(ly-1)*Lx+lx;
            delta_lxly = [(lx-1)*delta_x; (ly-1)*delta_y; 0];
            w_lo(l) = (2*pi/lambda_c).*(trm1_lo - trm2_lo)'*delta_lxly;
        end
    end

    z_lo(n) = (sum(cos(w_lo)))^2 + (sum(sin(w_lo)))^2;
end

% A*ln(1+B*z*ww) >= u   ---->  ww*z >= tt^2   tt^2 >= (exp(u/A)-1)/B
  tt_lo.^2 >= (exp(u/A)-1)./B;
  for n=1:N
    quad_over_lin(tt(n),z_lo(n)) <= ww_lo(n);
  end


%% ==================== Safety flying Distance ============================
% || qr - qj|| >= Ds;
norms(qj-qr).^2  >= Ds^2


%% ===================== Covertness constraint ============================
 % 1- min_w D * exp(-kf(||qr - qa|| + ||qr - qw||))/(||qr - qa|| * ||qr -
 % qw||)^(pl/2) < 1- epsilon; for all w <IN> K/{selected user k}

for n=1:N
    % ||qa-qr|| >= t1
    t1_lo(n) = sqrt(-norms(qa-qr(:,n)).^2 + 2*(qa-qr(:,n))'*(qa - qr(:,n)));

    % || qa-qw || >= t2
    qw = qk;
    qw(:,idx(n))=[];
    for m=1:numUsr-1
        t2_lo(m, n) = -norms(qw(:,m)-qr(:,n)).^2 + 2*(qw(:,m)-qr(:,n))'*(qw(:,m) - qr(:,n));
    end
end

 t1_lo >= 1
 t2_lo >= 1
 
 for m = 1:numUsr-1
    % t1 * t2 >= exp(yy)  exp(-kf(t1+t2))/(t1 * t2)^(pl/2) <= r
    % exp(-kf(t1+t2))/(t1 * t2)^(pl/2) <= exp(-kf(t1+t2))/exp(yy*pl/2) <= r
   yy_lo(m,:) = log(t1_lo) + log(t2_lo(m,:));
   r1(m,:) = D(m,:).*(L^2).* exp(-(kf*(t1_lo+t2_lo(m,:)) + (pl/2)*yy_lo(m,:)));
 end

  for m = 1:numUsr-1
    % t1 * t2 >= exp(yy)  exp(-kf(t1+t2))/(t1 * t2)^(pl/2) <= r
    % exp(-kf(t1+t2))/(t1 * t2)^(pl/2) <= exp(-kf(t1+t2))/exp(yy*pl/2) <= r
   yy_lo(m,:) = log(t1_lo) + log(t2_lo(m,:));
   r1(m,:) = D(m,:).*(L^2).* exp(-(kf*(t1_lo+t2_lo(m,:))));
end
r_lo = min(r1,[],1);

r_lo <= epsilon
%% =================== flight power non-convex term constraint ============
% t_lo = (200 - (Po*(1+c0*(norms(vr).^2))...
%       + c1*(norms(vr)).^3))/ Pi;
for n=1:N
ccc (n) = -t_lo(n).^2 +2*t_lo(n).*t_lo(n) + c2.*(-norms(vr(:,n)).^2 + ...
                                  2*vr(:,n)'*vr(:,n)) - power(t_lo(n),-2)
end
mean(ccc)
%% =================== Mobility constraints ===============================
    qr(:,1) == qr_I
    qr(:,N) == qr_I
    qr(:, 2:N) - qr(:,1:N-1) == [vr(:,1:N-1)*dt; zeros(1,N-1)]
    norms(vr) <= Vmax
    norms(vr(:,2:N) - vr(:,1:N-1)) <= Amax
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 