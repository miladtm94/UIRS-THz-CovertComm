Channels
[PF_r, PF_j] = flight_pow(vr,vj);
A = BW./log(2)./((PF_r + PF_j)/(Po+Pi)/2)/N;
%%
% c1 = sum(pa+pj) <= P_tot;
% for m=1:numUsr-1
%     c2(m)= mean(pa <= pj.*epsilon./D(m,:) );% N
% end
% c3 = mean(0 <= pj <= Pj_max); % N
% c4 = mean(0 <= pa <= Pa_max); % N
% 
% [c1 c2 c3 c4]
%%
% [pa_new, pj_new] = JointPower_optim(pa,pj,A,B,C,D, epsilon,Pj_max,Pa_max,P_tot, alpha);
[alpha_new, pa_new, pj_new] = JointPowUsrSch_optim(alpha, pa,pj, A, G_auk,G_jk, epsilon,Pj_max,Pa_max,P_tot);
pj = pj_new';
pa = pa_new';
alpha = alpha_new';

for n=1:N
   [~, idx(n)] = max(alpha(:,n));
end
AEE_Calc