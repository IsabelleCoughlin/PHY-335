%% Problem 2
%
% Integate Equation 1 with any method given $J_1$ at $x = 0$ to find $J_1(1)$.
% Dont use besselj or the recurrence relation to get initial condition.
% Instead use the series approximation for $J_v(x)$ for small x. 
%
%% Part 1
% Can you start the integrations at exactly $x = 0$?
%
% No, when the integrations start exactly at zero, there exists a
% singularity due to division by zero, and therefore the calculations
% cannot be done. We start with a small number near zero, since the series
% approximation allows us to use small x to approximate J_v(x).

clear
nvar = 2;
atol = 1.0e-9;
rtol = atol;
h1 = 0.00001;
hmin = 0.0;
x1 = 1e-3;
x2 = 1;
z = x1:0.1:x2;
counter = 1;

for i = x1:0.1:x2
    J(counter) = besselj(1,i);
    counter = counter+1;
end

%% Part 2
%
% If we integrate as we did in Problem 1, we see that the error difference between the 
% gets bigger as the x value increases. 

in = (x1/2) - (x1^3)/16 + (x1^5)/384;
in
de = 1/2 - (3*x1^2)/16 + (5*x1^4)/384;  
de
ys = [in, de];

sdp = NumericalRecipes.Output(9);
sbs = NumericalRecipes.Output(9);

m = HW6_P1Functor(0);

dedp = NumericalRecipes.Odeint(NumericalRecipes.StepperDopr5(),ys,x1,x2,atol,rtol,h1,hmin,sdp,m);
debs = NumericalRecipes.Odeint(NumericalRecipes.StepperBS(),ys,x1,x2,atol,rtol,h1,hmin,sbs,m);

yvlspd = dedp.integrate;
yvlsbs = debs.integrate;
axpd = sdp.xsave(1:sdp.count);
axbs = sbs.xsave(1:sbs.count);
aypd = sdp.ysave(1:sdp.count,1);
aybs = sbs.ysave(1:sbs.count,1);
azpd = sdp.ysave(1:sdp.count,2);
azbs = sbs.ysave(1:sbs.count,2);

% Plot all three lines
figure(1)
plot(z, J)
hold on
plot(axpd,aypd)
plot(axbs,aybs)
xlabel('x');
ylabel('y_2');
hold off
legend('bessel-true', 'apd','abs')

% How to find the differences
for i = 1:length(J)
    apddif(i) = J(i) - aypd(i);
    absdif(i) = J(i) - aybs(i);
end
% 
figure(2)
%plot(z, J)
hold on
plot(axpd,apddif)
plot(axbs,absdif)
xlabel('x');
ylabel('y_2');
title('Differences ')
hold off
legend('aRK','aBS')


%% Part 3
%
% Using this specific substitution of varaibles, w = (dy/dx)*x, we simplify
% the calculation of the derivative and aligns with the struction of the original
% equation due to the division by x. Therefore, at small x values when it would normally
% be really unstable, the scaling by x helps.
%
% Using the substitution w = (dy/dx)*x, we find the following

% y from power series solution
initial = (x1/2) - (x1^3)/16 + (x1^5)/384;
%w = dy/dx*x
deriv = (1/2)*x1 - ((3*x1^2)/16)*x1 + ((5*x1^4)/384)*x1;

% Initial Conditions
ystart = [initial, deriv];

outsdp = NumericalRecipes.Output(9);
outsbs = NumericalRecipes.Output(9);

d = HW6_P2Functor(0);

odedp = NumericalRecipes.Odeint(NumericalRecipes.StepperDopr5(),ystart,x1,x2,atol,rtol,h1,hmin,outsdp,d);
odebs = NumericalRecipes.Odeint(NumericalRecipes.StepperBS(),ystart,x1,x2,atol,rtol,h1,hmin,outsbs,d);

yvalspd = odedp.integrate;
yvalsbs = odebs.integrate;
xpd = outsdp.xsave(1:outsdp.count);
xbs = outsbs.xsave(1:outsbs.count);
ypd = outsdp.ysave(1:outsdp.count,1);
ybs = outsbs.ysave(1:outsbs.count,1);
zpd = outsdp.ysave(1:outsdp.count,2);
zbs = outsbs.ysave(1:outsbs.count,2);

figure(3)
plot(z, J)
hold on
plot(xpd,ypd)
plot(xbs,ybs)
xlabel('x');
ylabel('y_2');
hold off
legend('bessel-true', 'pd','bs')
for i = 1:length(J)
    pddif(i) = J(i) - ypd(i);
    bsdif(i) = J(i) - ybs(i);
end
%%
% The difference of the RK function in this case stays very close to zero
% compared the to BS function
figure(4)
%plot(z, J)
hold on
plot(xpd,pddif)
plot(xbs,bsdif)
xlabel('x');
ylabel('y_2');
title('Differences ')
hold off
legend('RK','BS')

%%
% Comparing the two functions we see a much more even difference, where it
% is getting larger as the value of x increases. Since the approximation fo
% series initial conditions is only valid for small x values, this makes
% sense. 

if outsdp.count == outsbs.count
    figure(5)
    plot(outsdp.xsave(1:outsdp.count),outsdp.ysave(1:outsdp.count,1)-outsbs.ysave(1:outsbs.count,1))
    xlabel('x');
    ylabel('\delta y_1');
end

