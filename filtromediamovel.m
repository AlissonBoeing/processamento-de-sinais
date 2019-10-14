r = 2*(rand(100000, 1)-0.5);

%N = 400;
%b = ones(N,1)/N;
N = length(b)-1;

i = 1;

x = zeros(N,1);
s = zeros(length(r),1);

while i< length(r)
    
    k = N;
    while k > 1
       x(k) = x(k-1);        
       k = k - 1;
    end
    
    x(1) = r(i);
    
    k = 1;
    while k <= N
        s(i) = s(i) + b(k)*x(k);
        k = k + 1;
    end
    
    i = i + 1;    
end

t = 1:1:size(r);

plot(t, r, ':b', t, s, '+r');

S = fftshift(fft(s)/length(s));
R = fftshift(fft(r)/length(r));

f = linspace(-1,1,length(r));

figure(2);
subplot(2,1,1);
plot(f, 20*log10(abs(R)));

subplot(2,1,2);
plot(f, 20*log10(abs(S)));

figure(3);
freqz(b,1);