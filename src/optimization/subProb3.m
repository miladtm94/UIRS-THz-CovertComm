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

A = BW./log(2)./(PF_r + PF_j);
B = pa.*g_aub;
C = g_jb/2;
D = pa.*G_auw./G_jw;

pj_new = JammingPow_optim(pj,A,B,C,D, epsilon,Pj_max, alpha);
pj = pj_new';

ITR = BW*log2(1+pa.*g_aub./(pj.*g_jb + 1));
ATR(ite_index) = mean(ITR);
[PF_r, PF_j] = flight_pow(vr,vj);
IFP= PF_r+ PF_j; 
APC(ite_index) = mean(IFP);
AEE(ite_index) = mean(ITR./IFP);