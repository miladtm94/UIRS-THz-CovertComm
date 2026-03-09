Channels

[PF_r, PF_j] = flight_pow(vr,vj);

A = alpha*BW./log(2)./((PF_r + PF_j)/(Po+Pi)/2)/N;
B = G_auk;
C = G_jk/2;
D = G_auw./G_jw;

%%
% [pa_new, pj_new] = JointPower_optim(pa,pj,A,B,C,D, epsilon,Pj_max,Pa_max,P_tot, alpha);
[pa, pj] = JointPowerWoRecast_optim(pa,pj,A,B,C,D, epsilon,Pj_max,Pa_max,P_tot, alpha);

AEE_Calc