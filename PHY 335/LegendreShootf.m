clear
params = [1,3.0,3.0]; % guesses: [P(-1),lambda(-1),lambda(1)]
Lload1 = LegendreLoadf1();
Lload2 = LegendreLoadf2();
Lscoref = LegendreScoref();
Lderivs = LegendreDerivs();
LegendraShootf = NumericalRecipes.Shootf(-1.0,1.0,0.0,Lload1,Lload2,Lderivs,Lscoref);
[params,check] = NumericalRecipes.newt(params,LegendraShootf);
if ~check
    fprintf('   Newt found lambda = %.9f\n',params(2));
    fprintf('             lambda2 = %.9f\n',params(3));
    fprintf('               P(-1) = %.9f\n',params(1));
else
    fprintf('  ?????? Newt found lambda = %.9f\n',lambda(1));
end
% Now that we have found the parameters of a solution, we need to 
% integrate the equations again using these parameters so that we can
% have the actual solution of the ODE.
% Note that for the initial conditions we must construct by hand the 
% same conditions coded within LegendreLoadf1 and LegendreLoadf2.
figure(1);
out1=NumericalRecipes.Output(-1);
out2=NumericalRecipes.Output(-1);
int1=NumericalRecipes.Odeint(NumericalRecipes.StepperDopr5(),[params(1),-params(2)*params(1)/2,params(2)],-1,0,1.e-8,1.e-8,0.001,0,out1,Lderivs);
int1.integrate();
int2=NumericalRecipes.Odeint(NumericalRecipes.StepperDopr5(),[1,params(3)/2,params(3)],1,0,1.e-8,1.e-8,0.001,0,out2,Lderivs);
int2.integrate();
plot(out1.xsave(1:out1.count),out1.ysave(1:out1.count,1),out2.xsave(1:out2.count),out2.ysave(1:out2.count,1));
