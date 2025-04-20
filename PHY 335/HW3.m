%% Homework 3
%
% Qromo is a function that implements Romberg integration using Midpnt or
% other extended versions. 
% function. It takes an instance of Midpnt class as input
% 
%% Problem 1: Integrate the following
%%
% First, we use the midpoint approximation with qromo to find the numeric
% solution to $\int_{0}^{\pi/2} \frac{\sqrt{x}}{sin(x)}dx$
t = NumericalRecipes.Midpnt(@(x) x^(1/2)/sin(x), 0, pi/2);
NumericalRecipes.qromo(t, 10^-5)

%%
% Now we integrate $\int_{\pi/2}^\infty \frac{sin(x)}{x^2}dx$
%
% This integral converges more slowly since sin(x) oscillates and $x^{-2}$ is not
% dominant enough to converge it quickly. 
s = NumericalRecipes.Midpnt(@(x) sin(x)/(x^2), pi/2, 10^99);
NumericalRecipes.qromo(s, 10^-4)
%%
% Obviously the solution is not correct, and we find that we can only
% acheive an accuracy of $10^{-4}$ without creating an error. Thus, there is an
% issue in the calculation using the original method. 
%
% We can use Midinf that will allow one of the bounds to be infinity and a
% 1/u variable change (very similar to the u-substitutio below). This
% improves the convergence and they are the same method in theory. For this
% method to work, either one limit can be going towards positive or
% negative infinity, but not both. This condition is satisfied for use in
% this example. 
%
r = NumericalRecipes.Midinf(@(x) sin(x)/(x^2), pi/2, 10^99);
NumericalRecipes.qromo(r, 10^-5)
%%
% Simple substitution can allow us to get finite limits which improves
% convergence when using numeric integration.
%
% We substitute x = 1/t into the equation so that when you integrate, 
% the bounds become t = 1/inf = 0, and t = 1/(pi/2) = 2/pi.
d = NumericalRecipes.Midpnt(@(t) sin(1/t), 0, 2/pi);
NumericalRecipes.qromo(d, 10^-5)
%%
% Now the function converges much quicker and produces the correct answer
% of 0.1658 and allows us to reach a higher accuracy of $10^{-5}$. 
%
% Lastly, we find $\int_0^\infty \frac{e^{-x}}{\sqrt{x}}dx$
%
u = NumericalRecipes.Midpnt(@(x) (exp(-x))/(x^(1/2)), 0, 10^99);
NumericalRecipes.qromo(u, 10^-5)

%% Problem 2: Integrate double integral over ellipse
%
% Integrate $\int \int e^{-xsin(y)}dxdy$ over the boundary given by $5x^{2} - 6xy +
% 5x^{2} = 2$
%
% We write a functor class that returns a function solely of x and creats
% the internal bounds for the y integral. First, we analytically solve for
% the outer bouds. 
%
% In order to solve for the bounds, we solve the equation for the ellipse
% as a function of y(x). We find $y(x) = \frac{(3x \pm \sqrt{(-16x^{2} +
% 10)}}{5}$. The two different solutions form our upper and lower bounds
% for the y integral. 
%
% Our bounds for the outer x integral occur when the functions for y become
% imaginary. This is because at those points of x, the y value cannot
% exist on the same plane that we're graphing on. Thus, solving for 
% Our bounds are found when the y(x) function becomes imaginary
% Solved for this analytically, we find $\pm \sqrt{\frac{5}{8}}$

boundtwo = (5/8)^(1/2);
boundone = -(5/8)^(1/2);

%%
% The computed integral is $\int_{-\sqrt{\frac{5}{8}}}^{\sqrt{\frac{5}{8}}} \int_{\frac{(3x - \sqrt{(-16x^{2} +
% 10)}}{5}}^{\frac{(3x + \sqrt{(-16x^{2} +
% 10)}}{5}} e^{-xsin(y)}dydx$
% As found in the "Profiler" this function runs 81 times. 
%%
l = NumericalRecipes.Midpnt(expfunc, boundone, boundtwo);
NumericalRecipes.qromo(l, 10^-5);
