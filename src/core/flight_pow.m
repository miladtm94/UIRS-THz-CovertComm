function [PF_r, PF_j] = flight_pow(vr,vj)

Flightconstants

PF_r = Po*(1+c0*norms(vr).^2) + c1*norms(vr).^3 ...
    + Pi*sqrt(sqrt(1+c2^2*norms(vr).^4) - c2*norms(vr).^2);

PF_j = Po*(1+c0*norms(vj).^2) + c1*norms(vj).^3 ...
    + Pi*sqrt(sqrt(1+c2^2*norms(vj).^4) - c2*norms(vj).^2);
end
