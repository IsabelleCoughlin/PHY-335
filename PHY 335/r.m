x_start = 1e-3; % Small x to avoid singularity
y0 = [x_start / 2; 1 / 2]; % Use series expansion

d = FunctionODE(0);

odedp = NumericalRecipes.Odeint(NumericalRecipes.StepperDopr5(),ystart,x1,x2,atol,rtol,h1,hmin,outsdp,d);
odebs = NumericalRecipes.Odeint(NumericalRecipes.StepperBS(),ystart,x1,x2,atol,rtol,h1,hmin,outsbs,d);

yvalspd = odedp.integrate;
yvalsbs = odebs.integrate;


[x, y] = NumericalRecipes.Odeint(@functor.dydx, y0, x_start, 1, NumericalRecipes.StepperDopr5());

% Compare with besselj
y_true = besselj(1, x);

figure;
plot(x, y(1,:), 'r', 'DisplayName', 'Numerical');
hold on;
plot(x, y_true, 'b', 'DisplayName', 'besselj(1,x)');
xlabel('x');
ylabel('J_1(x)');
legend;
title('Bessel Function J_1(x) Integration from x=0');
hold off;
