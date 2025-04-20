%% Problem 1
%
%% Part 1
%
% We integrate the Bessel function from x = 1 to x = 30. In addition, we
% compare the integrations using two methods, StepperDopr5 and StepperBS.
% It was found that StepperBS is much more accurate for this problem, and
% requires much fewer steps in order to get to that approximation. While BS
% only used 28 steps, RK (StepperDopr5) used nearly 4000 when given range
% in the Output function. This makes sense, since BS is a much more
% efficient routine and works well for smooth functions, which this is. 


clear
% Define problem parameters

% Number of variables
nvar = 2;

% Tolerances
atol = 1.0e-9;
rtol = atol;
% Step size parameters
h1 = 0.1;
hmin = 0.0;

% Initial Conditions
x1 = 1;
x2 = 30;
z = x1:0.1:x2;
counter = 1;

%Plotting True-Values of Bessel Function
for i = x1:0.1:x2
    J(counter) = besselj(1,i);
    counter = counter+1;
end
    
% Initital Conditions
J1 = besselj(1,1)  % Should be about 0.4401
% recurrence relation
J1_prime = (besselj(0,1) - besselj(1,1)/x1) % Should be close to 0.3251
ystart = [J1, J1_prime];

% Calculates number of steps (defined to 291 to match J-true function
outsdp = NumericalRecipes.Output(290);
outsbs = NumericalRecipes.Output(290);

% Call functor
d = HW6_P1Functor(0)

% Create odeint class
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

% Prints number of derivative calls
%d.count

%% Part 2
%
% We compare the results of the different methods using a variety of plots.
% 

%% 
% We plot the three series on the same graph, and find very good
% approximations. When zooming into very small scales, we see that there
% are differences. Highlighed much more clearly in Figure 4.

% Plot vs the True Values
figure(1)
plot(z, J)
hold on
plot(xpd,ypd)
plot(xbs,ybs)
xlabel('x');
ylabel('y_2');
hold off
legend('bessel-true', 'RK','BS')

%% 
% We calculate the differences between the true values and each of the
% series. We see in Figure 4 that the BS method is much more accurate at
% these scales. There is minute fluxuations when looked very closely, but
% they remain constant the entire scale from 1-30. On the other hand, for
% the RK function, the errors compound and get larger as the x-value
% increases, but fluxuates sinusoidally. 

% How to find the differences
for i = 1:length(J)
    pddif(i) = J(i) - ypd(i);
    bsdif(i) = J(i) - ybs(i);
end
% 
figure(2)
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
% Similarly, we can see the difference between the two approximations, and
% a similar pattern occurs, since BS approximation is so close to the true,
% it is like comparing the RK approximation to the true values again. 

%Differences between the two 
if outsdp.count == outsbs.count
    figure(3)
    plot(outsdp.xsave(1:outsdp.count),outsdp.ysave(1:outsdp.count,1)-outsbs.ysave(1:outsbs.count,1))
    xlabel('x');
    ylabel('\delta y_1');
end

%% Part 3
% 
% We compare the number of times the derivative is called for each routine.
% While in this finalized code they are computed simultaneously, when done
% separately for their most efficient step sizes (argument for .Output was (-1)
% it was found that RK required 22243 functor calls and BS required
% 1988 functor calls. This again shows how much more efficient BS is when
% calcualting the values, and is as much of a factor of 10 difference. 

% Derivative counts for BOTH functions at 290 steps.
d.count