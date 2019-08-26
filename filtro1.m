fp = 1e3;
fs = 2e3;
ws = 2*pi*fs;
wp = 2*pi*fp;

Ws = ws/wp;
Wp = wp/wp;

As = 40; % dB


Ap = 3; % dB

eps = sqrt((10^(0.1*Ap))-1);

n = (log((10^(0.1*As))-1))/(2*log(Ws));
n = ceil (n);

for k = 1:n
    p(k) = exp((i*pi*((2*k + n -1)/(2*n))));
end

D = poly(p);
 
[h,w] = freqs(1,D,1000);
subplot(211);
semilogx(w,20*log10 (abs(h)));
grid on;

subplot(212);
semilogx(w,unwrap(angle(h))/pi);
grid on;

