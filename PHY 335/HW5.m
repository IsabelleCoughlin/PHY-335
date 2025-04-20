%% Homework 5

%% Problem 1 A
% Using your choice of methods from the Numerical Recipes library, find the roots of:
% cos{x}−0.8+ px^{2} = 0 for $p = 0.1, 0.2, 0.3, 0.4$.
%
% I chose to use Brent's Method, which uses 3-point interpolation to find
% the roots, and implements zbrak to find the intial brackets to search
% within. 
%%
ps = [0.1, 0.2, 0.3];
for l = ps
    myfuncd = hwfivefunctor(l);
    [x1b, x2b] = NumericalRecipes.zbrak(myfuncd,-10,10,20);
    fprintf(['\n--- Finding roots using Brent''s Method ---\n']);
    for k=1:length(x1b)
        newtroot = NumericalRecipes.zbrent(myfuncd,x1b(k),x2b(k),1.e-12);
        fprintf('Bracket %2d: [%f,%f]; Root at x = %.16f\n',k,x1b(k),x2b(k),newtroot);
    end
end
%%
% We find roots for 0.1, 0.2, 0.3 shown above. The method cannot find roots for x = 0.4. Looking at the graph
% below, we see there are no real roots for p = 0.4.
%%
myfuncd = hwfivefunctor(0.4);
x=-10:0.01:10;
figure(1);
% Plots the function
plot(x,myfuncd.func(x));
xlabel('x');
ylabel('f(x)');
% Plots derivative of the function
figure(2);
plot(x,myfuncd.df(x));
xlabel('x');
ylabel('f''(x)');

%% Problem 1 B
% To find whether a function has a double root, we usually begin by finding values where f(x) = 0 and f'(x) = 0.
% Based on the pattern of how the p value affects
% the equation, a double root would exist between p = 0.3 and p = 0.4.
% After some analysis, we reduce this range to 0.3 to 0.35. We find the
% point at which the four roots distinctively become 2--two of the roots
% approach the same value (within a detectable range). 
%%
b = false;
for i=0.3:0.0001:0.35
    myfuncd = hwfivefunctor(i);
    [x1b, x2b] = NumericalRecipes.zbrak(myfuncd,-10,10,150);
    for k=1:length(x1b)
        rooted(k) = NumericalRecipes.zbrent(myfuncd,x1b(k),x2b(k),1.e-12);
        deriva(k) = myfuncd.df(newroot);
        end
    for t = 1:(length(x1b)-1)
        if abs((rooted(t+1) - rooted(t))) < (0.2)
            fprintf('\nRoot at x = %.10f, and derivative d = %.10f\n',i, deriva(t));
            b = true;
        end
    end
    if b
        break;
    end
end
%%
% At about p = 0.32, we find a double root where the function grazes the
% x-axis. Therefore, after this point, the function does not cross the
% x-axis at all. To verify that this is a double root, we also see that the
% derivative of the value is also approaching zero. If we wanted to show
% that this is exactly zero, we would need to find a much more precise
% value of p--on the scale closer to 0.3231... Graphing this value, we see
% the function is near the x-axis and the value of the derivative at 0 is
% equal to zero (truncated from the true value for the purposes of the graph). 
%%
myfuncd = hwfivefunctor(0.3231);
x=-10:0.01:10;
figure(1);
% Plots the function
plot(x,myfuncd.func(x));
xlabel('x');
ylabel('f(x)');
% Plots derivative of the function
figure(2);
plot(x,myfuncd.df(x));
xlabel('x');
ylabel('f''(x)');
%% Problem 2
%
% Use zroots to find the roots of x^6 − 12.1*x^5 + 59.5*x^4 − 151.85*x^3 +
% 212.6625*x^2 − 156.6*x + 48.5625 = 0

a = [48.5625,  -156.6, 212.6625, -151.85, 59.5, -12.1, 1];
polish = false;
roots = NumericalRecipes.zroots(a, polish);
roots

%%
% Roots are the exact same whether you polish true or false