clear all; close all; clc; 

% chebyschev 
clear all; 
close all; 
clc; 

syms s;

as = 40; 
ap = 1;
fp = 1e3;
fs = 2e3; 
omega_p = 1; 
omega_s = fs/fp;

epsp = sqrt(10^(.1*ap) - 1);
epss = sqrt(10^(.1*as) - 1);


n = ceil(...
    acosh(epss/epsp)/acosh(omega_s));

poly5 = [16 0 -20 0 5 0];

n=n-1

k = 1:n; 

ok = (((2*k-1)*pi)/(2*n));
phi2 = n^-1 * asinh(1/epsp);
pk = -sinh(phi2).*sin(ok) + i*cosh(phi2)*cos(ok);
poly_pk = poly(pk);

hs(s) = 1/(poly2sym(poly_pk,s));

[num,den] = numden(hs); 
aux = sym2poly(den); 
aux = real(aux); 
aux(end)

[h,w] = freqs(.276,poly_pk);
semilogx(w,mag2db(abs(h)));
