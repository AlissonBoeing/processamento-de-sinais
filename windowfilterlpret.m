%% passa baixa
m = 50;
n= -m:m;

wc = 0.2*pi;

h = sin(wc.*n)./(pi.*n)

h(m+1) = (wc/pi)

stem(n,h)
%freqz(h,1)
%zplane(h,1)

%% passa alta
m = 50;
n= -m:m;

wc = 0.2*pi;

h = -(sin(wc.*n)./(pi.*n))

h(m+1) = 1- (wc/pi)

stem(n,h)


%% passa faixa

m = 50;
n= -m:m;

wc1 = 0.5*pi;
wc2 = 0.55*pi;

h = (1/pi.*n).*(sin(wc2.*n) - sin(wc1.*n))

h(m+1) = (wc2-wc1)/(pi)

stem(n,h)