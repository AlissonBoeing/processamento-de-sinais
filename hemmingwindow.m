%% passa baixa
m = 50;
n= -m:m;

wc = 0.2*pi;

w = 0.5 + 0.5.*cos((2*pi.*n)/(2*m + 1))

h = (sin(wc.*n)./(pi.*n)).*w



h(m+1) = (wc/pi)*w(m+1)

stem(n,h)
%freqz(h,1)
%zplane(h,1)