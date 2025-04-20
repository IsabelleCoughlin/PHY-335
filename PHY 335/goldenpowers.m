phi(1) = single(1);
phi(2) = (sqrt(5)-1)/2;
phim(1) = single(1);
phim(2) = (sqrt(5)-1)/2;
error(1) = abs(phi(1)-phim(1))/phim(1);
error(2) = abs(phi(2)-phim(2))/phim(2);
for i = 3:40
    phi(i) = phi(i-2) - phi(i-1);
    phim(i) = phim(2)*phim(i-1);
    error(i) = abs(phi(i)-phim(i))/phim(i);
end
figure(1) 
semilogy(error)

% growth of the error is linear on a logarithmic plot meaning it is exponential

% subtraction/roundoff error growns worst case linearly (aka not this)

% this algorithm is numerically unstable

for i = 1:40
    phidub(i, 1) = phi(i);
    phidub(i, 2) = phim(i);
    phidub(i, 3) = error(i);
end

