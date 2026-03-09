syms alpha Pa A

%-P1
f1 = @(alpha, Pa, A) alpha*log(1+Pa*A);
f1d= gradient(f1,[alpha, Pa]);
f1d2 = hessian(f1,[alpha,Pa]);

pretty(simplify(f1d))
pretty(simplify(f1d2))

%%
syms  x y k b a z
f2 = @(x,y,k, b, a) log(1+z*exp(-k*(x+y))/x^a/y^a);
f2d= gradient(f2,[x, y, z]);
f2d2 = hessian(f2,[x,y, z]);

pretty(simplify(f2d))
pretty(simplify(f2d2))
simplify(det(f2d2))

%%

syms  z a t
f2 = @(z,t, a) log(1+z*exp(a*t));
f2d= gradient(f2,[z,t]);
f2d2 = hessian(f2,[z, t]);

pretty(simplify(f2d))
pretty(simplify(f2d2))
simplify(det(f2d2))

  
%%
syms  u z
f3 = @(u,z) (exp(u)-1)/z;
f3d= gradient(f3,[x, y, z]);
f3d3 = hessian(f3,[x,y, z]);

pretty(simplify(f3d))
pretty(simplify(f3d3))
simplify(det(f3d3))

%%
syms  x y k a b c
f4 = @(x,y) log(1+a*exp(-b*(x+y))/(x^(c))/(y^(c)));
f4d= gradient(f4,[x, y]);
f4d2 = hessian(f4,[x,y]);
assume(a,{'real', 'positive'})
assume(b,{'real', 'positive'})
assume(c,{'real', 'positive'})
pretty(simplify(f4d2,'IgnoreAnalyticConstraints', true,'Steps',50))
pretty(simplify(det(f4d2)))
%%
syms  x k a
f5 = @(x) x^(a/2)*exp(k*(x));
f5d= gradient(f5,[x]);
f5d5 = hessian(f5,[x]);

pretty(simplify(f5d))
pretty(simplify(det(f5d5)))
%%

syms x s alpha kf
syms x_lo s_lo
f6 = @(x,s) (x*s)^(-alpha/2)*exp(-kf*(x+s));
f6d= gradient(f6,[x,s]);
f6d2 = hessian(f6,[x,s]);
pretty(simplify(f6d))
pretty(simplify(f6d2))
pretty(simplify(det(f6d2)))

simplify(taylor((x*s)^(alpha/2)*exp(kf*(x+s)),[x, s],[x_lo s_lo],'Order',2))


simplify(taylor((x)^(-alpha/2)*exp(-kf*x)),[x],[x_lo],'Order',2)

simplify(taylor((x*s)^(alpha/2)*exp(kf*(x+s)),[x, s],[x_lo s_lo],'Order',2))

%%
syms  x y k a
f4 = @(x,y,k, a) exp(-k*(x+y))/(x^(a))/(y^(a));
f4d= gradient(f4,[x, y]);
f4d4 = hessian(f4,[x,y]);

pretty(simplify(f4d))
pretty(simplify(det(f4d4)))


%%


syms a x y
f = @(w)  log(1+a*x/y^2);
fd= gradient(f,[x,y]);
fd2 = hessian(f,[x,y]);
pretty(simplify(fd))
pretty(simplify(fd2))
% pretty(simplify(taylor(f(w),[w],[w_lo],'Order',2)))



%%  jammer Traj
syms A B C w w_lo W
f = @(w) A*log(1+B/(C*w+1));
fd= gradient(f,[w]);
fd2 = hessian(f,[w]);
pretty(simplify(fd))
pretty(simplify(fd2))
pretty(simplify(taylor(f(w),[w],[w_lo],'Order',2)))



syms s a kf 
f = @(s) s^a*exp(kf * s);
fd= gradient(f,[s]);
fd2 = hessian(f,[s]);
pretty(simplify(fd))
pretty(simplify(fd2))



%%
syms A B t 
f = @(t) A*log(1+B*t^2);
fd= gradient(f,[t]);
fd2 = hessian(f,[t]);
pretty(simplify(fd))
pretty(simplify(fd2))

%%
%%
syms x y a x_lo y_lo
f = @(x,y) 1./((x*y)^a);
fd= gradient(f,[x,y]);
fd2 = hessian(f,[x,y]);
pretty(simplify(fd))
pretty(simplify(fd2))
pretty(simplify(det(fd2)))

pretty(simplify(taylor(f(x,y),[x,y],[x_lo, y_lo],'Order',2)))

%%

syms x y 
f = @(x,y) log(1+(1/x)*exp(y));
fd= gradient(f,[x,y]);
fd2 = hessian(f,[x,y]);
pretty(simplify(fd))
pretty(simplify(fd2))
pretty(simplify(det(fd2)))
pretty(simplify(taylor(f(x,y),[x,y],[x_lo, y_lo],'Order',2)))



%%

syms x y pl kf
syms x_lo y_lo
f = @(x,y) (x*y)^(-1)*exp(-(kf/pl)*(x+y));
fd= gradient(f,[x,y]);
fd2 = hessian(f,[x,y]);
pretty(simplify(fd))
pretty(simplify(fd2))
pretty(simplify(det(fd2)))

taylor(f(x,y),[x, y],[x_lo y_lo],'Order',2)

pretty(simplify(taylor(f(x,y),[x, y],[x_lo y_lo],'Order',2)))


%% 

syms s a kf 
f = @(s) s^a/exp(-kf * s);
fd= gradient(f,[s]);
fd2 = hessian(f,[s]);
pretty(simplify(fd))
pretty(simplify(fd2))


%%
syms   t t_lo
f = @(t) log(1+1/t);
fd= gradient(f,[t]);
fd2 = hessian(f,[t]);
pretty(simplify(fd))
pretty(simplify(fd2))
pretty(simplify(taylor(f,[t],[t_lo],'Order',2)))

%%
syms  a b t y
f = @(t,y) log(1+(a*t^2)/(1+b*y));
fd= gradient(f,[t,y]);
fd2 = hessian(f,[t,y]);
pretty(simplify(fd))
pretty(simplify(fd2))
pretty(simplify(det(fd2)))
%%


syms x y x_lo y_lo
f = @(x,y) 1/x/y

pretty(hessian(f,[x,y]))
taylor(f,[x,y],[x_lo,y_lo],'Order',2)
pretty(taylor(f,[x,y],[x_lo,y_lo],'Order',2))

syms s_lo b s
g = @(s) log(1+b*exp(-s));
pretty(taylor(g,s,s_lo,'Order',2))
taylor(g,s,s_lo,'Order',2)

%%

syms  x x_lo k a b c
f = @(x) log(1+a*exp(-b*x)/(x^c));
fd= gradient(f,[x]);
fd2 = hessian(f,[x]);
pretty(simplify(fd2))

pretty(simplify(taylor(f,x,x_lo,'Order',2)))

%%
syms  x x_lo a b c pl kf
f = @(x) a*log(1+b/(c*(exp(-kf*x)/(x^pl))+1));
fd= gradient(f,[x]);
fd2 = hessian(f,[x]);
pretty(simplify(fd2))

pretty(simplify(taylor(f,x,x_lo,'Order',2)))

%%

syms s s_lo pj pj_lo

f = @(s,pj) (1-s+pj)^2
fd= gradient(f,[s,pj])
pretty(simplify(taylor(f,[s pj],[s_lo, pj_lo],'Order',2)))

