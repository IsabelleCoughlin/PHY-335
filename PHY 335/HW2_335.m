%% HW 2: Interpolation
%% *Problem 1: Polynomial Interpolation with Numerical Recipes Library* 
%
% (a) Consider evaluating the function $y = frac{1}{(1 + 25*x.^2)}$ by polynomial 
% interpolation on the interval [−1, 1] using equally spaced points. Use Poly
% interp with 5 and 20 points. Compare values from Poly interp with the exact 
% values, especially near ±1. Comment on what you find—a graph would be nice.
%
%% Creating the datasets
% We use the function FuncGen, previously defined, to create random x, y
% pairs based on a passed function. We define parameters, such to create N
% number of points between (inclusive) A and B. The fourth input parameter
% allow you to pass coefficients for a polynomial function if f(x) is a
% polynomial. Since out original function, $y = frac{1}{(1 + 25*x.^2)}$, is
% rational, we pass an empty list. The last input parameter is for error to
% create  uniform noise within the data set, but that is left at zero. 
%%
clear
A = -1;
B = 1;
N = 5;
[x,y] = FuncGen(A,B,N,[],0);
[xp,yp] = FuncGen(A,B,4*N,[],0);
%%
% Now we have two sets of data points, one with 5 and one with 20 points. 

%%% Interpolation using an M-th order polynomial
%
% Define a 4th order polynomial and interpolate based on both datasets with
% varying number of points of the same function. We use the built-in
% Numerical Recipes function of Poly_interp().
%
M = 4;
p_5 = NumericalRecipes.Poly_interp(x,y,M);
p_20 = NumericalRecipes.Poly_interp(xp,yp,M);

%%%
% Generate finite set of x values over which to interpolate
%
h = (B - A)/99;
xi = A:h:B;


%%% Interpolation plots using original data
%
% We create an interpolated fit based on the five data points and plot them
% against the original function (approximately made with the 20-point data
% and the five original points). 
%

%%
for i=1:length(xi)
    [ypi(i),p_5] = p_5.interp(xi(i));
    errsp(i)=p_5.dy;
end
figure(1);
plot(xp,yp,'k',x,y,'r*',xi,ypi,'b');
xlabel('x');
ylabel('y');
legend('y(x)','exact','fit');
title('Interpolation fit using 5 data points');
set(gca,'fontsize', 15);

%%
% We can clearly see that this function is not a relatively good fit for
% our data. Especially near the end points, the function is much more
% curved and fit to a polynomial than our goal of the rational function.
% The limiting to only five data points does not give enough information to
% adequately fit to the correct function with this method. 
%
% We plot the interpolated fit based on the function generated with 20
% points, against the original data. 

%%
for i=1:length(xi)
    % ypi(i) = p.interp(xi(i));
    [yppi(i),p_20] = p_20.interp(xi(i)); % save error info
    errspp(i)=p_20.dy;
end
figure(2);
plot(xp,yp,'k',x,y,'r*',xi,yppi, 'b')
xlabel('x')
ylabel('y')
legend('y(x)','exact','fit')
title('Interpolation fit usign 20 data points')
set(gca,'fontsize', 15)

%%
% We can see that this function based off many more original data points
% gives a much better fit to the data. While many times this can lead to
% overfitting, it does not seem to be at that point. 
%
%% Error Data plots using original data
%
% Figure 3 shows plots with error bars for data shown in Figure 1, aka the
% data generated from 5-points. While it was obvious tha thte function was
% not a good fit in Figure 1, we also see that the error is more
% substantial in some areas than others. Especially closer to the end
% points, the error increases a lot. 

%%
figure(3);
errorbar(xi,ypi,errsp)
xlabel('x')
ylabel('y')
title('5-Point Interpolation Error')
legend('exact')
set(gca,'fontsize', 15)

%%
% Figure 4 plots error from Figure 2 of the 20-point interpolation. The
% error is much lower in this figure and shows that the interpolated fit is
% accurate for the data. 

%%
figure(4);
errorbar(xi,yppi,errspp)
xlabel('x')
ylabel('y')
title('20-Point Interpolation Error')
legend('error')
set(gca,'fontsize', 15)

%% Finding coefficients of interpolated function
%
% (b) Find the coefficients of the interpolating polynomials using polcof and/or polcoe. 
% Evaluate the function from the coefficients. Comment.
%
% Two different Numerical Recipes methods are compared to calculate the
% coefficients of the polynomial fitted using poly_interp(). polcof()uses
% the division method, with $\sigma(N.^3)$ to calculate the coefficients,
% while polcoe() uses vandermond matrix method with an error of
% $\sigma(N.^2)$. While polcof() generally gives better results, it has
% issues with division by zero at certain points. 
%
% Here, we compare the two methods and see similar results for the N = 5
% and N = 20 cases. 
%
% For the first graph, we see that the function is not entirely smooth, since it is plotted at
% only 20 points, although the data is generated from a smooth polynomial.
% This plot is similar to the fit found using the 5-point data, since this
% polynomial was generated to fit that plot. 

%%
format long
% We must flip the arrays since the coefficients are listed in the opposite
% way as polcof and polcoe read the input.
cofa = flip(NumericalRecipes.polcof(x,y));
figure(5);
plot(xp,polyval(cofa,xp), xi, ypi)
xlabel('x')
ylabel('y')
title('Polcof() Coefficient Polynomial')
subtitle('5-point data')
legend('polynomial', 'fit')
set(gca,'fontsize', 15)

%%

coea = flip(NumericalRecipes.polcoe(x,y));
figure(6);
plot(xp,polyval(coea,xp), xi, ypi)
xlabel('x')
ylabel('y')
title('Polcoe() Coefficient Polynomial')
subtitle('5-point data')
legend('polynomial', 'fit')
set(gca,'fontsize', 15)

%%
% We see no great difference between the two methods in for this function.
%
% The following graphs show the polynomials created for the N = 20 data
% values. These graphs fit much closer to the original interpolated fit,
% but with minor noticalble differences (probably majority due to it being
% only evaluated at 20 step sizes instead of 99 as the fit function was.
%%
cof = flip(NumericalRecipes.polcof(xp,yp));
figure(7);
plot(xp,polyval(cof,xp), xi, yppi)
xlabel('x')
ylabel('y')
title('Polcof() Coefficient Polynomial')
subtitle('20-point data')
legend('polynomial', 'fit')
set(gca,'fontsize', 15)

%%
coe = flip(NumericalRecipes.polcoe(xp,yp));
figure(8);
plot(xp,polyval(coe,xp), xi, yppi)
xlabel('x')
ylabel('y')
title('Polcoe() Coefficient Polynomial')
subtitle('20-point data')
legend('polynomial', 'fit')
set(gca,'fontsize', 15)

%%
% (c) Now use Rat interp with 5 and 20 points. Compare values from Rat interp with the 
% exact values. Comment on what you find—a graph would be nice
%
% We follow the same procedure as above, but with the Rat_interp()
% function from Numerical Recipes instead. We can immediately see in both
% graphs that the fit is much closer to our original equation for both
% trials of N = 5 and N = 20. Especially for N = 5, the fit was much closer
% and comparable to the fit at N = 20 for both the Rat_interp() and
% Poly_interp() fit. This shows that using the Rat_interp requires much
% less data to be accurate for rational functions as opposed to
% Poly_interp(). While this seems like a relatively obvious conclusion, it
% is shown here with the comparison. 

%%
% Using rat_interp function from Numerical Recipes
r_5 = NumericalRecipes.Rat_interp(x,y,M);
r_20 = NumericalRecipes.Rat_interp(xp,yp,M);

%%
% Plotting the rat-interp() interpolation from 5 points
%%
for i=1:length(xi)
    [yri(i),r_5] = r_5.interp(xi(i));
    err_r_5(i)=r_5.dy;
end
figure(9);
plot(xp,yp,'k',x,y,'r*',xi,yri,'b');
xlabel('x');
ylabel('y');
legend('y(x)','exact','fit');
title('Rat Interp fit using 5 data points');
set(gca,'fontsize', 15);

%%
% Plotting the rat-interp() interpolation from 20 points
%%
for i=1:length(xi)
    [yrri(i),r_20] = r_20.interp(xi(i));
    err_r_20(i)=r_20.dy;
end
figure(10);
plot(xp,yp,'k',x,y,'r*',xi,yrri,'b');
xlabel('x');
ylabel('y');
legend('y(x)','exact','fit');
title('Rat Interp fit using 20 data points');
set(gca,'fontsize', 15);


