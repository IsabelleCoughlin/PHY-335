%% Homework 4
%% Problem 1
%
% Do by Monte Carlo integration of $\int \int e^{-xsin(y)}dxdy$ over the 
% boundary given by $5x^{2} - 6xy + 5x^{2} = 2$. Therefore, we create a
% rectangular box to surround out bounds of integration, namely the
% ellipse, and generate random points within that rectangle. If
% we conclude the point lies inside the ellipse, we evaluate the point with
% our function and include it in our integral. 
%
% In order to make the function more efficient, we transform our ellipse to
% a different coordinate system in order to make the rectangle fit more
% closely to the bounds. In this case, we transform to the (u,v)
% coordinates. After evaluating whether or not they are in the integral,
% we either add or do not add the point into the calculation of the
% integral by transforming them back into (x, y) coordinates and using our
% original equation $e^{-xsin(y)}$.
%%
clear;
ud=NumericalRecipes.Ran(454533);
j=0;
int = 0;
Ntot=50000;
area = 2*1;
u_min = -1;
u_max = 1;
v_min = -0.5;
v_max = 0.5;
for i=1:Ntot
    ut=-1 + 2*ud.doub;
    vt=(-1/2) + 1*ud.doub;
    if(ut^2 + (vt^2)*4 <= 1)
        xt = (1/sqrt(2))*(ut-vt);
        yt = (1/sqrt(2))*(ut+vt);
        int=int+exp(-xt*sin(yt));
        j = j+1;
        x(j)=xt;
        y(j)=yt;
    end
    if mod(i,100)==0
        figure(1)
        scatter(x,y)
        axis equal
        pause(0.001)
    end
end
figure(1)
scatter(x,y)
set(gca,'fontsize', 36)
axis equal

%% Plots
% The plot shows how we fill in the data points slowly until we have an
% adequate coverage of our integral shape. We see that it is also the same
% shape as our elliptical bounds, which suggests the transformations
% correctly calculated the bounds. We now find the total integral value,
% which we intermediately calculated through a summation during the
% process. 

%% Total Integral
% We found the total integral to be 1.4478 using 50000 points. The true
% value is 1.4490, which relates to a 0.0828% error percentage. At only 500
% points, we get a value of 1.4553, which is a 0.4348% error, also falling
% below the 0.5% threshhold. Hence, we don't need a ton of points in order
% to get a relatively good error percentage. In addition, the error
% percentage get increasingly harder to reduce as we increase the number of
% points. This suggests that error of the data converges very slowly
% compared to the number of points we integrate with. 
total = area*int/Ntot
j