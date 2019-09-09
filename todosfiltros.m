%MEU CODIGO - filtros


clear all;
close all;
clc;

syms s; %simbolico s
N = 0; % N =  ordem do filtro

as = 35;
ap = 2;
fp = 3.4e3;
fs = 4e3;
omega_p = 1;
Op = omega_p;
omega_s = fs/fp;
Os = omega_s;

%Subs de variaveis
Wp = Op; % Omega P
Ws = Os; % Omega S
Ap = ap; % Atenuaçao de descida ou subida
As = as; % Atenuação de stop




%% Butter
[N,Wn] = buttord(Wp, Ws, Ap, As,'s');
[b,a] = butter(N,Wn, 's');

%% Chebyshev I
N = cheb1ord(Wp, Ws, Ap, As,'s');
[b,a] = cheby1(N, Ap, Wp, 's');

%% Chebyshev II
N = cheb2ord(Wp, Ws, Ap, As,'s');
[b,a] = cheby2(N,As, Ws, 's');

%% Elliptic - Cauer
[N, Wn] = ellipord(Wp, Ws, Ap, As,'s');
[b,a] = ellip(N,Ap,As, Wn, 's');

%% Transformação de frequência
syms p;
Hp = (poly2sym(b,p))/(poly2sym(a,p));
Hs = subs(Hp,(s/(2*pi*fp))); % Alterar fp para mudar a freq.

[bs,as] = numden(Hs);
[h,w] = freqs(real(sym2poly(bs)),real(sym2poly(as)));

subplot(211)
semilogx(w/(2*pi),mag2db(abs(h))); 
hold on; 
grid;

plot([fp,fs], -[Ap, As], 'xk');%marcar X preto em (fp,fs) e -(ap,as) 
hold off;
subplot(212)
semilogx(w/(2*pi),unwrap(angle(h))/pi); 
grid;
