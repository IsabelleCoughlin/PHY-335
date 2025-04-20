clear
lambda(1) = 3; % initial guess for the eigenvalue
even = true; % true->Even : false->odd
Lload = LegendreLoad();
Lscore = LegendreScore(even);
Lderivs = LegendreDerivs();
LegendraShoot = NumericalRecipes.Shoot(1.0,0.0,Lload,Lderivs,Lscore);
[lambda,check] = NumericalRecipes.newt(lambda,LegendraShoot);
if ~check
    fprintf('   Newt found lambda = %.9f\n',lambda(1));
else
    fprintf('   ??????  Newt found lambda = %.9f\n',lambda(1));
end
figure(1);
out=NumericalRecipes.Output(100);
int=NumericalRecipes.Odeint(NumericalRecipes.StepperBS(),[1,lambda(1)/2,lambda(1)],1,0,1.e-8,1.e-8,0.001,0,out,Lderivs);
int.integrate();
plot(out.xsave(1:out.count),out.ysave(1:out.count,1));
hold on
flip=1;
if ~even
    flip=-1;
end
plot(-out.xsave(1:out.count),flip*out.ysave(1:out.count,1));
hold off