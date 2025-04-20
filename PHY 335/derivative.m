n = 61;
xmax = 3*pi;
dx = single(xmax/(n-1));
for i = 1:n
    x(i) = (i-1)*dx;
    y(i) = sin(x(i));
end

hold


figure(1)
plot(x,y)
xlabel('x')
ylabel('sin(x)')
set(gca, 'fontsize', 18)

%return;

h = 0.1;

for i = 1:n
    dydx = (sin(x))
end

figure(2)

err = cos(x) - dydx;

figure(3)
plot(x, err)
xlabel('x')
ylabel('error')
set(gca, 'fontsize', 36)

errnorm = single(0);
for i = 1:n
    errnorm = err(i)*err(i);
end
endnorm = sqrt(errnorm/n);

disp('L2 norm of error:')