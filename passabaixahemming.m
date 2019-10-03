%% passa baixa
m = 15;
n= -m:m;


wp = 0.2*pi;
ws = 0.2*pi;
Ap = 0.2;
As = 50;
wc = sqrt(ws*wp);

%teta_c = wc/(wa/2);

w = 0.5 + 0.5.*cos((2*pi.*n)/(2*m + 1)); %hemming

%w = 0.5 + 0.5.*cos((2*pi.*n)/(2*m +1)); %hann

h = (sin(wc.*n)./(pi.*n)).*w;

h(m+1) = (wc/pi).*w(m+1);

h = h*10^((-Ap/2)/20); %correcao 

[H, w] = freqz(h,1);

wsm = sum(mag2db(abs(H))>-0.2)



%wsm = 0.766; %medidos
%wpm = 0.33;

wse = 0.3*pi; %especificados
wpe = 0.2*pi;

dwm = (wsm - wpm);

n2 = ((dwm)*(m*2))/(wse-wpe);


stem(n,h)
figure(1)
freqz(h,1)
hold on
plot([0 wp wp], -[Ap Ap As+20], ':r');
hold on;
plot([0 ws ws 1],-[0 0 As As], ':m');
hold off;

%zplane(h,1)
