theta = zeros(L,N);
for lx = 1:Lx
    for ly = 1:Ly
        l=(ly-1)*Lx+lx;
        delta_lxly = [(lx-1)*delta_x; (ly-1)*delta_y; 0];
        temp = (2*pi*(qr-qa)'*delta_lxly)./lambda_c./reshape(norms(qr-qa), N, 1);
         theta(l, :) = temp';
    end
end
Ea = exp(-1i*theta);

beta = zeros(numUsr, L, N);
for usr=1:numUsr
    for lx = 1:Lx
        for ly = 1:Ly
            l=(ly-1)*Lx+lx;
            delta_lxly = [(lx-1)*delta_x; (ly-1)*delta_y; 0];
            temp = (2*pi*(qk(:,usr)-qr)'*delta_lxly)./lambda_c./reshape(norms(qk(:,usr)-qr), N, 1);
            beta(usr, l, :) = temp';
        end
    end
        
end

Ek = exp(-1i*beta); 

H_ark = zeros(numUsr,N);
for usr = 1: numUsr
%--- BS-IRS-user channel gains
    d_ar = norms(qa-qr);
    d_rk =  norms(qk(usr)-qr);
    d_ark = d_ar + d_rk;
    H_ark(usr,:) = hr0.*exp(-1i*2*pi*d_ark/lambda_c).*exp(-kf.*d_ark/2).*(d_ar.*d_rk).^(-pl/2);
%--- Jamming user channel gain
    d_jk =  norms(qk(usr)-qj);
    H_jk(usr,:) = hj0.*exp(-1i*2*pi*d_jk/lambda_c).*exp(-kf.*d_jk/2).*(d_jk).^(-pl/2);  

end

[PF_r, PF_j] = flight_pow(vr,vj);

A = BW./log(2)./(PF_r+ PF_j)/2; % 1-by-N
B= zeros(1,N);
for n = 1:N
    B(n) = abs(H_ark(idx(n),n)).^2*pa(n)./(0.5*pj(n)*abs(H_jk(idx(n),n)).^2 + 1); % 1-by-N
end
Ch = abs(H_ark).^2.*repmat(pa, numUsr, 1)./repmat(pj, numUsr,1)./(abs(H_jk).^2);
C =  zeros(numUsr-1, N); % (K-1)-by-N

for n=1:N
    tmp = Ch(:,n);
    tmp(idx(n)) = [];
    C(:,n) = tmp;      
end

% min(C*L^2,[],1) <= epsilon
%%

parfor n=1:N
    Phi(:,n) = IRSbeamforming_optim(C, Ek, Ea, idx, epsilon,n);
%     Phi(:,n) = abs(V).*exp(1i*wrapTo2Pi(-angle(V)));
%    Phi(:,n) = abs(V).*exp(1i*wrapTo2Pi(+angle(V)));
%     Phi(:,n) = V;
end

%%
AEE_Calc
