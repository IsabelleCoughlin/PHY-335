%% Problem 3
%
% We compare the solution to integration of the differential equation $y' =
% -\lambda*y$, using both explicit and implicit Euler's method. We keep the
% value of $\lambda$ constant, but vary the step size to investigat the
% accuracy and stability of the solutions. 
%
%% Part A
% We use 1000 equal steps to integrate the solution, and compare it to the
% analytical solution, which is known to be $y = Ce^{-\lambda x}$

% Initial Conditions
x0 = 0;
y0 = 1;
% Test n = 1000
n = 1000;
h = 1/n;
lambda = 100;

% Allocate arrays and set intial conditions
x_val = zeros(1,n);
x_val(1,1) = x0;

y_val_e = zeros(1,n);
y_val_e(1,1) = y0;

y_val_i = zeros(1,n);
y_val_i(1,1) = y0;
% Counted for array indexing
c = 2;
% Iterate and calculate the values of the function
for i = 0:h:1
    % Explicit Calcualtion
    y_val_e(1,c) = y_val_e(1,c-1) + h*(-lambda*y_val_e(1,c-1));
    % Implicit Calculation
    y_val_i(1,c) = y_val_i(1,c-1)/(1 + lambda*h);
    % X is the same for both since the step is the same
    x_val(1,c) = x_val(1,c-1) + h;
    c = c+1;
end
%%

% Plot explicit vs implicit vs analytical
figure(1)
plot(x_val,y_val_e)
hold on
plot(x_val,y_val_i)
plot(x_val, y0*exp(-lambda*x_val))
xlabel('x_val');
ylabel('y_val');
hold off
legend('explicit-1', 'implicit-1', 'analytical-1')
%%
% In the graph above, we can see that the solutions are relatively similar,
% but the scale makes it a big difficult to see differences, therefore we
% transform it below. 
%%

% log(y) plot
x_val_log = log(x_val);
figure(2)
plot(x_val_log,y_val_e)
hold on
plot(x_val_log,y_val_i)
plot(x_val_log, y0*exp(-lambda*x_val))
xlabel('x_val');
ylabel('y_val');
hold off
legend('log_explicit-1', 'log_implicit-1', 'log_analytical-1')
%%
% In this grpah, we can see that the two solutions are relatively close in
% behavior to the analytic solution. It is also noticable that each
% solution is about the same amount off on each side o the analytic
% solution, one overestimating and one underestimating. This suggests that
% the $n = 100$ steps allows both solutions to be accurate and stable in
% this case. 
%

%% Part B
% In order for the explicit Euler's method to be stable, $h \ge 2/\lambda$.
% Since we set $\lambda = 100$ this implies that if $h \ge 0.02$ that the
% method is stable, but if $h > 0.02$ then the method is unstable. We find
% the minimum steps for the explicit Euler method to be stable is $n = 50$.
% 
%% Part C
% Correspondingly, to test this boundary, we change the conditions to test
% $n = 49$ ($h = 0.204...$), $n = 50$ ($h = 0.02$), and $n = 51$ ($h =
% 0.196...$). Hence, we expect for $n = 50$ that the function will remain
% constantly erraneous (not decay or grow), but for $n = 51$ we expect
% exponential decay. Yet, if $n = 49$, the function will be unstable and
% will grow exponentially. We see the graphs below in order $n = 49, 50, 51$.
% 

figure_count = 3;

for f = 49:1:51
    x_val = zeros(1,f);
    x_val(1,1) = x0;
    
    y_val_e = zeros(1,f);
    y_val_e(1,1) = y0;
    
    y_val_i = zeros(1,f);
    y_val_i(1,1) = y0;

    h = 1/f;
    p = 2;
    for i = 0:h:1
        y_val_e(1,p) = y_val_e(1,p-1) + h*(-lambda*y_val_e(1,p-1));
        y_val_i(1,p) = y_val_i(1,p-1)/(1 + lambda*h);
        x_val(1,p) = x_val(1,p-1) + h;
        p = p+1;
    end
    figure(figure_count)
    plot(x_val,y_val_e)
    hold on
    plot(x_val,y_val_i)
    plot(x_val, y0*exp(-lambda*x_val))
    xlabel('x_val');
    ylabel('y_val');
    hold off
    legend('explicit 3', 'implicit 3', 'analytical')
    figure_count = figure_count+1;
end

%% 
% As expected, we see an exponental increase in $n = 49$ representing an
% unstable method implementation. Interestingly for this problem, at
% exactly $n = 50$ is the boundary of stability for step size changes.
% Therefore, we can see the behavior at exactly the boundary and see that
% it is not increasing or decreasing exponentially. While it is not a
% "good" solution, it is stable. 