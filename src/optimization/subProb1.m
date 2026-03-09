Channels
Akn = G_auk./ (0.5*pj.*G_jk + 1);
[pa_new, alpha_new] = UsrSchPow_optim(pa, pj, alpha, Pa_tot, ...
                               Pa_max, Akn, G_auk,G_jk, epsilon);

alpha = alpha_new';
for n=1:N
   [~, idx(n)] = max(alpha(:,n));
end
pa = pa_new';

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

ATR(ite_index) = mean(BW*log2(1+pa.*g_aub./(pj.*g_jb + 1)));

[PF_r, PF_j] = flight_pow(vr,vj);
AFP(ite_index) = mean([PF_r, PF_j])/100; % normalized by some constant value

AEE(ite_index) = ATR(ite_index)/AFP(ite_index);
