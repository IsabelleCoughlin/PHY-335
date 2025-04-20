clear
N = 2^11;
d = (1-(-1))/(N-1);
t = -1:d:1;
n = -N/2: d: N/2;
f = n/(N*d);
w = 25;
f0 = 25;

h = 2*w*sqrt(pi)*exp(-1*(pi^2)*(w^2)*(t.^2)).*cos(2*pi*f0*t);

figure(1)
plot(t, h)
%axis([-0.05 0.05 -20 100])


