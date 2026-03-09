function [PF_r, PF_j] = flight_pow_cvxapprox(vr,vj,t)

Flightconstants

PF_r = Po*(1+c0*norms(vr).^2) + c1*norms(vr).^3 + Pi*t;

PF_j = Po*(1+c0*norms(vj).^2) + c1*norms(vj).^3 + Pi*t;
end
