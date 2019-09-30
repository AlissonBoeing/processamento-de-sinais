%% passa baixa
m = 5;
n= -m:m;

fp = 10e3;
fs = 15e3;
fa = 50e3;

wp = 2*pi*fp;
ws = 2*pi*fs;
wa = 2*pi*fa;

wc = sqrt(ws*wp);

teta_c = wc/(wa/2);

%w = 0.5 + 0.5.*cos((2*pi.*n)/(2*m + 1)) %hemming

w = 0.5 + 0.5.*cos((2*pi.*n)/(2*m +1)); %hann

h = (sin(teta_c*pi.*n)./(pi.*n)).*w;

h(m+1) = (teta_c).*w(m+1);

wsm = 0.766; %medidos
wpm = 0.33;

wse = 0.6; %especificados
wpe = 0.4;

dwm = (wsm - wpm);

n2 = ((dwm)*(m*2))/(wse-wpe);


stem(n,h)
freqz(h,1)
%zplane(h,1)