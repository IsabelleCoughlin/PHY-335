% Examples of Root Finding:
%

% Finding Brackets

%Zbrac
% Expand a range until a bracket is found:
[foundbrac, x1, x2] = NumericalRecipes.zbrac(@sin,pi/4,3*pi/4);
if foundbrac
    fprintf('Bracket [%f,%f] found\n',x1,x2);
end

% Zbrak
% Divide a range into N pieces and return all brackets
N=20;
[x1a, x2a] = NumericalRecipes.zbrak(@sin,-5*pi/2,5*pi/2,N);
fprintf('\nFound %d brackets\n',length(x1a));
for k=1:length(x1a)
    fprintf('Bracket %2d: [%f,%f]\n',k,x1a(k),x2a(k));
end

% Bisection (slow)
fprintf('\n--- Finding roots using Bisection ---\n');
for k=1:length(x1a)
    [foundbrac,x1,x2] = NumericalRecipes.zbrac(@sin,x1a(k),x2a(k));
    % need to expand bracket since we may end up with roots at bracket
    % limits.
    if foundbrac
        bisroot = NumericalRecipes.rtbis(@sin,x1,x2,1.e-12);
        fprintf('Bracket %2d: [%f,%f]; Root at x = %.16f\n',k,x1a(k),x2a(k),bisroot);
    end
end

% False Position Method
fprintf('\n--- Finding roots using False Position Method ---\n');
for k=1:length(x1a)
    fpsroot = NumericalRecipes.rtflsp(@sin,x1a(k),x2a(k),1.e-12);
    fprintf('Bracket %2d: [%f,%f]; Root at x = %.16f\n',k,x1a(k),x2a(k),fpsroot);
end

% Secant Method
fprintf('\n--- Finding roots using Secant Method ---\n');
for k=1:length(x1a)
    secroot = NumericalRecipes.rtsec(@sin,x1a(k),x2a(k),1.e-12);
    fprintf('Bracket %2d: [%f,%f]; Root at x = %.16f\n',k,x1a(k),x2a(k),secroot);
end

% Brent's Method
% Best method using just f(x)
fprintf('\n--- Finding roots using Brent''s Method ---\n');
for k=1:length(x1a)
    brentroot = NumericalRecipes.zbrent(@sin,x1a(k),x2a(k),1.e-12);
    fprintf('Bracket %2d: [%f,%f]; Root at x = %.16f\n',k,x1a(k),x2a(k),brentroot);
end

%return

% Newton's Method
% We are required to use a "FunctorD" class in order to define both a
% function and its first derivative.
%
% TestFunctionD.m defines a cubic polynomial with coefficients
% a x^3 + b x^2 + c^x +d

% Passes in the coefficients to the functor
myfuncd = hwfivefunctor(0.32);
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
% Find some brackets between -10 and 10
[x1b, x2b] = NumericalRecipes.zbrak(myfuncd,-10,10,150);
% 
% % Simple Newton's Method
% % Need sufficiently good guess at start
% fprintf('\n--- Finding roots using Newton''s Method ---\n');
% for k=1:length(x1b)
%     newtroot = NumericalRecipes.rtnewt(myfuncd,x1b(k),x2b(k),1.e-12);
%     fprintf('Bracket %2d: [%f,%f]; Root at x = %.16f\n',k,x1b(k),x2b(k),newtroot);
% end
% 
% % Safe Newton's Method
% % Leaves a slow/non-convergent cycle
% fprintf('\n--- Finding roots using Safe Newton''s Method ---\n');
% for k=1:length(x1b)
%     saferoot = NumericalRecipes.rtsafe(myfuncd,x1b(k),x2b(k),1.e-12);
%     fprintf('Bracket %2d: [%f,%f]; Root at x = %.16f\n',k,x1b(k),x2b(k),saferoot);
% end
% 
% % Compare to Brent's Method
% fprintf('\n--- Finding roots using Brent''s Method ---\n');
% for k=1:length(x1b)
%     newtroot = NumericalRecipes.zbrent(myfuncd,x1b(k),x2b(k),1.e-12);
%     fprintf('Bracket %2d: [%f,%f]; Root at x = %.16f\n',k,x1b(k),x2b(k),newtroot);
% end
