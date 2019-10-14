%% Filtro passa-faixa (Frequência de passagem 12 - 19 Hz) Cheby1
% Atenuação na banda de passagem = 1
% Atenuação na banda de rejeição = 50

%%
syms s; %simbolico s
N = 0; % N =  ordem do filtro

As = 20;
Ap = 1;
fp1 = 10e3;
fp2 = 19e3;
fs1 = 10e3;
fs2 = 25e3;

wp1 = fp1*2*pi;
ws1 = fs1*2*pi;
wp2 = fp2*2*pi;
ws2 = fs2*2*pi;

B = (wp2 - wp1);


w0 = sqrt(wp2*wp1);
omega_p = 1;
Op = omega_p;
omega_s1 = abs((((-ws1^2) + (w0^2)))/(B*ws1));
omega_s2 = abs((((-ws2^2) + (w0^2)))/(B*ws2));
Os = min(omega_s1,omega_s2)

%Subs de variaveis
Wp = Op; % Omega P
Ws = Os; % Omega S





%% Butter
%[N,Wn] = buttord(Wp, Ws, Ap, As,'s');
%[b,a] = butter(N,Wn, 's');

%% Chebyshev I
figure (1)
N = cheb1ord(Wp, Ws, Ap, As,'s')
[b,a] = cheby1(N, Ap, Wp, 's');
[h,w] = freqs(b,a);
semilogx(w/(2*pi),mag2db(abs(h))); 


%% Chebyshev II
%N = cheb2ord(Wp, Ws, Ap, As,'s');
%[b,a] = cheby2(N,As, Ws, 's');

%% Elliptic - Cauer
%[N, Wn] = ellipord(Wp, Ws, Ap, As,'s');
%[b,a] = ellip(N,Ap,As, Wn, 's');

%% Transformação de frequência
syms p;

P = (((s^2) + (w0^2)))/(B*s);



Hp = (poly2sym(b,p))/(poly2sym(a,p));

Hs = subs(Hp,P); % Alterar fp para mudar a freq.

[bs,as] = numden(Hs);
[h,w] = freqs(real(sym2poly(bs)),real(sym2poly(as)), 10000);

figure(2)
subplot(211)
semilogx(w/(2*pi),mag2db(abs(h))); 
ylim([-80 5])
hold on; 
grid;

plot([fp1,fs1,fp2,fs2], -[Ap, As, Ap, As], 'xk');%marcar X preto em (fp,fs) e -(ap,as) 
hold off;
subplot(212)
semilogx(w/(2*pi),unwrap(angle(h))/pi); 
grid;