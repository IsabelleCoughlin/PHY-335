%% Homework 7

% We find the series of eigenvalues of Mathieus equation based on even and
% off functions. We use the functor class (HW7_Functor) to define our
% equation and its derivatives. The equation: $\frac{d^2y}{dt^2} +
% (a-2*qcos(2*t))*y = 0$, where q is a defined and passed parameter. We
% explore the eigenvalue sequence for both $q = 5$ and $q = 25$, for both
% even and odd function parameters of a. 
%
%  We perform shooting function to solve the boundary value problem. My
%  approach to this problem was to calculate values of a for the four
%  different cases, periodic even, antiperiodic even, periodic odd, and
%  antiperiodic odd. In order to do this, we must relate the "n" value for
%  each solution of eigenvalue to the value of a found. The value of n
%  corresponds to the number of roots in the solution the found value of a
%  represents. As such, the value of the eigenvalue increases as n
%  increases. Therefore, when n = 0, the lowest value of a is found and n
%  has no roots. We can also gleen information about the function from the
%  value of n/number of roots. For n = 0, the function cannot cross the
%  x-axis, and hence must be periodic. Then, n = 1 crosses the x-axis once
%  but does not return to where it began, so it must be anti-periodic.
%  Hence, we know n even values are periodic and n odd values are
%  antiperiodic. The approach to finding a_n() and b_n() series at certain
%  values of n included estimating the starting eigenvalue solution and
%  searching for a solution with the correct number of roots. The values
%  and graphs of these solutions are saved on a separate pdf, but were
%  found with the following methods. 
%
% For simplicity, all of the functions were saved separately for the four
% different cases. Hence, there us a Score and Load file for each of the
% four cases, described below. 
%%
% For shooting boundary problems, we instantiate with an initial guess, and
% the algorithm finds the closest solution value to that guess. 
clear
params = [-5]; % initial guess for the eigenvalue

%%
% Load function defines boundary conditions
Lload = MathieuLoad();
% Score evaluates when the boundary conditions are met
Lscore = MathieuScore();
% Functor defines the function
Lderivs = HW7_Functor(25);
% Shoot performs the shooting method until the score is met
MShoot = NumericalRecipes.Shoot(0,pi,Lload,Lderivs,Lscore);
[params,check] = NumericalRecipes.newt(params,MShoot);
if ~check
    fprintf('   Newt found lambda = %.9f\n',params(1));
else
    fprintf('   ??????  Newt found lambda = %.9f\n',params(1));
end
figure(1);
out=NumericalRecipes.Output(-1);
% We integrate the function to find the graph for y(t)
int=NumericalRecipes.Odeint(NumericalRecipes.StepperDopr5(),[1, 0, params(1)],0,pi,1.e-8,1.e-8,0.001,0,out,Lderivs);
int.integrate();
% Plot the function
plot(out.xsave(1:out.count),out.ysave(1:out.count,1));
title(['Periodic Even, n = 2, a(n)', num2str(params(1))])

%% Evaluate other versions
%
% Above, we evaluated the periodic even versions of this function. In order
% to find other solutions, we must slightly change the functions to find
% different conditions. Below we explore antiperiodic even, periodic odd,
% and antiperiodic odd. 

params1 = [220]; % initial guess for the eigenvalue
%even = true; % true->Even : false->odd
Lload = MathieuLoad_b();
Lscore = MathieuScore_b();
Lderivs = HW7_Functor(25);
MShoot = NumericalRecipes.Shoot(0,pi,Lload,Lderivs,Lscore);
[params1,check] = NumericalRecipes.newt(params1,MShoot);
if ~check
    fprintf('   Newt found lambda = %.9f\n',params1(1));
else
    fprintf('   ??????  Newt found lambda = %.9f\n',params1(1));
end
figure(2);
out=NumericalRecipes.Output(-1);
int=NumericalRecipes.Odeint(NumericalRecipes.StepperDopr5(),[1, 0, params1(1)],0,pi,1.e-8,1.e-8,0.001,0,out,Lderivs);
int.integrate();
plot(out.xsave(1:out.count),out.ysave(1:out.count,1));
title(['AntiPeriodic Even, n = 15, a(n)', num2str(params1(1))])



params2 = [100]; % initial guess for the eigenvalue
%even = true; % true->Even : false->odd
Lload = MathieuLoad_c();
Lscore = MathieuScore_c();
Lderivs = HW7_Functor(5);
MShoot = NumericalRecipes.Shoot(0,pi,Lload,Lderivs,Lscore);
[params2,check] = NumericalRecipes.newt(params2,MShoot);
if ~check
    fprintf('   Newt found lambda = %.9f\n',params2(1));
else
    fprintf('   ??????  Newt found lambda = %.9f\n',params2(1));
end
figure(3);
out=NumericalRecipes.Output(-1);
int=NumericalRecipes.Odeint(NumericalRecipes.StepperDopr5(),[0,1, params2(1)],0,pi,1.e-8,1.e-8,0.001,0,out,Lderivs);
int.integrate();
plot(out.xsave(1:out.count),out.ysave(1:out.count,1));
title(['Periodic Odd, n = 10, b(n)', num2str(params2(1))])



params3 = [225]; % initial guess for the eigenvalue
Lload = MathieuLoad_d();
Lscore = MathieuScore_d();
Lderivs = HW7_Functor(5);
MShoot = NumericalRecipes.Shoot(0,pi,Lload,Lderivs,Lscore);
[params3,check] = NumericalRecipes.newt(params3,MShoot);
if ~check
    fprintf('   Newt found lambda = %.9f\n',params3(1));
else
    fprintf('   ??????  Newt found lambda = %.9f\n',params3(1));
end
figure(4);
out=NumericalRecipes.Output(-1);
int=NumericalRecipes.Odeint(NumericalRecipes.StepperDopr5(),[0,1, params3(1)],0,pi,1.e-8,1.e-8,0.001,0,out,Lderivs);
int.integrate();
plot(out.xsave(1:out.count),out.ysave(1:out.count,1));
title(['Anti-Periodic Odd, n = 15, b(n)', num2str(params3(1))])