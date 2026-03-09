theta = zeros(L,N);
for lx = 1:Lx
    for ly = 1:Ly
        l=(ly-1)*Lx+lx;
        delta_lxly = [(lx-1)*delta_x; (ly-1)*delta_y; 0];
        temp = (2*pi*(qr-qa)'*delta_lxly)./lambda_c./reshape(norms(qr-qa), N, 1);
         theta(l, :) = temp';
    end
end

beta = zeros(L, N);
for n=1:N
    for lx = 1:Lx
        for ly = 1:Ly
            l=(ly-1)*Lx+lx;
            delta_lxly = [(lx-1)*delta_x; (ly-1)*delta_y; 0];
            temp = (2*pi*(qk(:,idx(n))-qr(:,n))'*delta_lxly)./lambda_c./norms(qk(:,idx(n))-qr(:,n));
            beta(l, n) = temp;
        end
    end
        
end

Phi = exp(-1i*(theta - beta));

Channels
 
ITR = alpha*BW.*log2(1+repmat(pa,numUsr,1).*G_auk./(0.5*repmat(pj,numUsr,1).*G_jk + 1));
ATR(ite_index) = mean(sum(ITR,1));
[PF_r, PF_j] = flight_pow(vr,vj);
IFP= (PF_r+ PF_j)/(Po+Pi)/2; 
APC(ite_index) = mean(IFP);
AEE(ite_index) = mean(sum(ITR./repmat(IFP,numUsr,1),1));