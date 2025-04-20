% Example of integrating a Stiff ODE using
% explicit Euler vs implicit Forward Euler
%
clear u v x yvars % Clear the variable in Workspace!
%
% The following integrates the ODEs using the usual explicit Euler method.
%
h = 0.001;
x1 = 0;
x2 = 4;
n=ceil((x2-x1)/h);
x(1)=x1;
u(1)=1.0;
v(1)=0.0;
for i=1:n
    u(i+1) = u(i) + h*(998*u(i)+1998*v(i));
    v(i+1) = v(i) + h*(-999*u(i)-1999*v(i));
    x(i+1) = x1 + i*h;
end
% Analytic solution
ua = 2*exp(-x) - exp(-1000*x);
va = -exp(-x) + exp(-1000*x);
figure(1)
plot(x,u,x,ua,'black:','LineWidth',2)
title('Euler (explicit)');
xlabel('x');
ylabel('u');
axis([x1 x2 -2 4]);
figure(2)
plot(x,v,x,va,'black:','LineWidth',2)
title('Euler (explicit)');
xlabel('x');
ylabel('v');
axis([x1 x2 -4 2]);
%
% Now lets integrate the same equations with the same step size using the
% Implicit Forward Euler method.
%
x(1)=x1;
yvars(:,1)=[u(1); v(1)];
det = 1+h*1001+h^2*1000;
onepluschinv = (1/det)*[1+h*1999 h*1998; -h*999 1-h*998];
for i=1:n
    yvars(:,i+1)=onepluschinv*yvars(:,i);
    x(i+1) = x1 + i*h;
end
figure(3)
plot(x,yvars(1,:),x,ua,'black:','LineWidth',2)
title('Forward Euler (implicit)');
xlabel('x');
ylabel('u');
axis([x1 x2 -2 4]);
figure(4)
plot(x,yvars(2,:),x,va,'black:','LineWidth',2)
title('Forward Euler (implicit)');
xlabel('x');
ylabel('v');
axis([x1 x2 -4 2]);