function AEE = myfunc_NLP(x)

pa_opt = x(:,1);
pj_opt = x(:,2);
b1=3+numUsr-1;
alpha_opt = x(:,3:b1);
b2=3+numUsr-1+2*L;
phi_opt = x(:,b1+1:b2);
qr_opt = x(:,b2+1:b2+2);
qj_opt = x(:,b2+3:b2+4);
vr_opt = x(:,b2+5:b2+6);
vj_opt = x(:,b2+7:b2+8);



theta = zeros(L,N);
for lx = 1:Lx
    for ly = 1:Ly
        l=(ly-1)*Lx+lx;
        delta_lxly = [(lx-1)*delta_x; (ly-1)*delta_y; 0];
        temp = (2*pi*(qr_opt-qa)'*delta_lxly)./lambda_c./reshape(norms(qr_opt-qa), N, 1);
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
            temp = (2*pi*(qk(:,usr)-qr_opt)'*delta_lxly)./lambda_c./reshape(norms(qk(:,usr)-qr_opt), N, 1);
            beta(usr, l, :) = temp';
        end
    end
        
end

Ek = exp(-1i*beta);

H_ark = zeros(numUsr,N);
for usr = 1: numUsr
    d_ar = norms(qa-qr_opt);
    d_rk =  norms(qk(usr)-qr_opt);
    d_ark = d_ar + d_rk;
    H_ark(usr,:) = hr0.*exp(-1i*2*pi*d_ark/lambda_c).*exp(-kf.*d_ark/2).*((d_ar.*d_rk).^(-pl/2));

    h_ark = H_ark(usr,:);
    for n=1:N
        ek = reshape(Ek(usr,:,n),L,1);
        ea = reshape(Ea(:,n),L,1);
        a = diag(ek')*ea;
        A = a*a';
        u = reshape(phi_opt(n,1:L),L,1)+1i*reshape(phi_opt(n,L+1:2*L),L,1);
        G_auk(usr,n) =  abs(h_ark(n))^2 * abs(u'*A*u)^2;
    end
    
%--- Jamming user channel gain
    d_jk =  norms(qk(usr)-qj_opt);
    H_jk(usr,:) = hj0.*exp(-1i*2*pi*d_jk/lambda_c).*exp(-kf.*d_jk/2).*(d_jk).^(-pl/2);  

end
G_jk = H_jk.*conj(H_jk);


ITR = min(BW*alpha_opt.*log2(1+repmat(pa_opt,numUsr,1).*G_auk./(0.5*repmat(pj_opt,numUsr,1).*G_jk + 1)),[],1);
[PF_r, PF_j] = flight_pow(vr_opt,vj_opt);
IFP= (PF_r+ PF_j)/(Po+Pi); 
AEE = mean(ITR./IFP);
end