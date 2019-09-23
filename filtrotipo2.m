%% type 2
b = [-1 -1 2 1 1 2 -1 -1]
b = [b flip(b)]

figure(1)
freqz(b,1)

figure(2)
subplot(211)
grpdelay(b,1)
subplot(212)
zplane(b,1)

