clear all;
close all;
clc;

fa = 8e3;%freq amostragem
as = 10;%atenuaçao stop
ap = 1;%atenuaçao passagem
fp = 941;
fs = 1209;
wp = 2*pi*fp;
ws = 2*pi*fs;

%transformação bilinear
%z = exp(1j*ws);
%s = 2*fa*( (z-1) / (z+1) );
N = 0; % N =  ordem do filtro



teta_p = (fp*2)/fa; %teta = fp/(fa/2)
teta_s = (fs*2)/fa;
lambda_p = 2*tan((teta_p*pi)/2);
lambda_s = 2*tan((teta_s*pi)/2);

omega_p = 1;
Op = omega_p;
omega_s = lambda_s/lambda_p;
Os = omega_s;

%Subs de variaveis
Wp = Op; % Omega P
Ws = Os; % Omega S
Ap = ap; % Atenuaçao de descida ou subida
As = as; % Atenuação de stop




%% Butter
[N,Wn] = buttord(Wp, Ws, Ap, As,'s');
[b,a] = butter(N,Wn, 's');


%% Transformação de frequência
syms s;
syms p;
Hp = (poly2sym(b,p))/(poly2sym(a,p));
Hs = subs(Hp,(s/(2*pi*fp))); % Alterar fp para mudar a freq.
pretty(Hp)


%%
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