SystemParams   % loading system parameters
ite_index = 1;
Feasible_Init  % initilizing a feasible point

pa_opt = optimvar('pa_opt',N,'LowerBound',0,'UpperBound',Pa_max);
pj_opt = optimvar('pa_opt',N,'LowerBound',0,'UpperBound',Pj_max);

alpha_opt = optimvar('alpha_opt',N,numUsr,'Type','integer','LowerBound',0,'UpperBound',1);

phi_opt = optimvar('phi_opt',N,2*L);

qr_opt = optimvar('qr_opt',N,2);
qj_opt = optimvar('qj_opt',N,2);

vr_opt = optimvar('vr_opt',N,2);
vj_opt = optimvar('vj_opt',N,2);

vars = {'P1','P2','I1','I2','C','LE1','LE2','HE1','HE2',...
    'HPS','MPS','LPS','BF1','BF2','EP','PP'};

x = {pa_opt,pj_opt,alpha_opt,phi_opt,qr_opt,qj_opt,vr_opt,vj_opt};


rng default % For reproducibility
gs = GlobalSearch;
x0 = [zeros(N,2*L+numUsr+2),qr(1:2,:)',qj(1:2,:)',vr',vj'];
problem = createOptimProblem('fmincon','x0',x0,...
    'objective',myfunc_NLP(x));
x = run(gs,problem)

