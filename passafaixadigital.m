%FILTRO DIGITAL


clear all;
close all;
clc;

fa = 8e3;
as = 30;
ap = 1;

fp1 = 811;
fp2 = 895.5;
fs1 = 770;
fs2 = 1209;

wp1 = 2*pi*fp1;
ws1 = 2*pi*fs1;
wp2 = 2*pi*fp2;
ws2 = 2*pi*fs2;


%transformação bilinear
%z = exp(1j*ws);
%s = 2*fa*( (z-1) / (z+1) );

teta_p1 = (fp1*2)/fa; %teta = fp/(fa/2)
teta_s1 = (fs1*2)/fa;
teta_p2 = (fp2*2)/fa; %teta = fp/(fa/2)
teta_s2 = (fs2*2)/fa;


lambda_p1 = 2*tan((teta_p1*pi)/2);
lambda_s1 = 2*tan((teta_s1*pi)/2);
lambda_p2 = 2*tan((teta_p2*pi)/2);
lambda_s2 = 2*tan((teta_s2*pi)/2);
Wp = 1;

lambda_s = min(lambda_s1,lambda_s2);

B = lambda_p2 - lambda_p1;

lambda_0 = sqrt(lambda_p2*lambda_p1);

Ws = abs((-(lambda_s^2) + lambda_0^2)/(B*lambda_s));




%Subs de variaveis
Ap = ap; % Atenuaçao de descida ou subida
As = as; % Atenuação de stop


%Butter
[N,Wn] = buttord(Wp, Ws, Ap, As,'s');
[b,a] = butter(N,Wn, 's');


syms s;
syms p;

P = ((s^2) +(lambda_0^2))/(B*s);

Hp = (poly2sym(b,p))/(poly2sym(a,p));
pretty(vpa(Hp,3))
[h,w] = freqs(b,a,10e3);
semilogx(w,mag2db(abs(h)));
title('Hp')
grid on
ylim([-20 0])
hold on
plot([Wp,Ws], -[Ap, As], 'xk');
hold off
figure(2)
zplane(b,a)
title('Zplane - Hp')

%% Transformacao de frequencia

Hs(s) = subs(Hp,P); % Alterar fp para mudar a freq.
pretty(vpa(Hs,3))

[bs_5,as_5] = numden(Hs);
bs = sym2poly(bs_5);
as = sym2poly(as_5);
[h,w] = freqs(bs,as,10e3);
figure(3)
semilogx(w,mag2db(abs(h)));
title('Hs')
grid on
ylim([-20 0])
hold on
plot([lambda_p1, lambda_p2,lambda_s1, lambda_s2], -[Ap, As], 'xk');
hold off
figure(4)
zplane(bs,as)
title('Zplane - Hs')

%% Transformaçao bilinear 
syms z
Hz(z) = subs(Hs,(2*1*(z-1)/(z+1))); 
pretty(vpa(collect(Hz),3))
[bz_5,az_5] = numden(collect(Hz));
bz = sym2poly(bz_5);
az = sym2poly(az_5);

%normalização
an = az(1);
bzn = bz/an;
azn = az/an;


[h,w] = freqz(bzn,azn,1000);
figure(5)
plot((w/(2*pi))*fa,mag2db(abs(h)));
title('Hz')
ylim([-20 0])
grid on
hold on
plot([fp1, fp2,fs1, fs2], -[Ap, As], 'xk');
hold off
figure(6)
zplane(bzn,azn)
title('Zplane - Hz')
