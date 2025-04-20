clear
lambda=8;
h=0.0125; % stepsize
xfin=5; % end-point of integration
%
x=zeros(1); % reset the x and y arrays so we don't
y=zeros(1); % have to clear Workspace when we change h
%
i=1;
x(i)=0; % start at x=0
y(i)=1; % start with y(0)=1
while x(i)<xfin
    y(i+1)=y(i)-h*lambda*y(i); % Euler's method
    x(i+1)=i*h;
    i=i+1;
end
ya=y(1)*exp(-lambda*x); % compute analytic solution
yrelerr=abs(ya-y)./ya; % compute relative error
figure(1)
plot(x,ya,'k',x,y,'r') % plot numeric & analytic sol.
title('Analyic and Numerical solutions');
xlabel('x');
ylabel('y(x)');
figure(2)
semilogy(x,ya,'k',x,abs(y),'r') % plot numeric & analytic sol.
title('Analyic and |Numerical| solutions');
xlabel('x');
ylabel('y(x)');
figure(3)
loglog(x,yrelerr,'r') % plot error
title('|Analyic - Numerical|/Analytic');
xlabel('x');
ylabel('relativ error');
