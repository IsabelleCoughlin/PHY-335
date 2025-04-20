%% Homework 1

%% Problem 1: Taylor Series
% *Explore truncation error in series expansion of sin(x)*

%%% Introduction
% In this file, we explore the error due to value truncation, by comparing
% an approximation of sin(x) to the actual value. We calculate the
% taylor series expansion of sin(x) at varying values of N, the number of
% steps for which we approximate. As n -> inf, the denominator of the
% expression approaches zero. These slight adjustments bring the total
% approximation closer to the true valule--except for when the value
% becomes too small for the computer to handle (aka the datatype in which
% the number is stored does not have the capacity to store as many decimal
% points as needed). Thus, at some value of N the accuracy levels off. At 
% this point though, the value of each term being added or subtracted would
% be so low to make any effects on the overall approximation due to
% roundoff. Thus, we can assume that this will not be a large source of error.
% On the other hand, when the values of x become very low (approaching
% zero),  their corresponding values of sin(x) become similarly low (as
% sin(x) = 0). Thus, this source of error is also unlilkely to create a
% large problem. 
%
% Generally, the trend we find is that the error is noticeably higher
% for the combination of low values of N and high values of X. This reason
% for such will be explained in more detail throughout.


%%% Defining bounds and Evaluating Error
% The main way to locate and analyze error for this function is to
% experiment with a variety of bounds for x and N. I found that the
% following are adequate for showing the patterns in error. In this
% section, we also calculate the error between the approximation and the
% real values of sin(x). This is done in a few ways to explore the
% underlying relationships. First, in the error array we store the
% numerical difference between the approximation and the actual value. In
% the following two arrays, we store the l2_norm of that error. The first
% one is calculated over N, meaning we can analyze the behavior of N as it
% relates to the error. In the other, we calculate the same for X. 
%
x = single(0:pi/20:pi); % Values of x to evaluate sin(x)
n = single(2:10); % Values of N

% Create an array to store approximated sin(x) values
sinSeries_array = single(zeros(length(x), length(n)));
% Create an array to store error between approximation and actual
error = single(zeros(length(x), length(n)));
% Arrays to store different L2-Norm values
l2_error_overn = zeros(1, length(n));
l2_error_overx = zeros(1, length(x));

% Loop through all values of x and n
for x_val  = 1:length(x)
    for n_val = 1:length(n)
        % Call sinseries approximation
        sinSeries_array(x_val, n_val) = sinSeries(x(x_val),n(n_val));
        % Calculate direct error
        error(x_val, n_val) = sinSeries_array(x_val, n_val) - sin(x(x_val));
    end
end

% Calculating the L2_Norm for Error over N
for n_val = 1:length(n)
    sum_over_N = 0;
    for x_val = 1:length(x)
        sum_over_N = sum_over_N + error(x_val, n_val)^2;
    end
    l2_error_overn(n_val) = sqrt(sum_over_N / length(x)); % Normalize
end

% Calculating the L2_Norm for Error over X
for x_val = 1:length(x)
    sum_over_X = 0;
    for n_val = 1:length(n)
        sum_over_X = sum_over_X + error(x_val, n_val)^2;
    end
    l2_error_overx(x_val) = sqrt(sum_over_X / length(n)); % Normalize
end

%%% Graphs
% We define many graphs that all provide useful information for analyzing
% the trends of approximation and error for the calculations. While many
% are commented out since they provide similar information, the most useful
% ones are kept and shown below. 
%

% Display sin(x)_approx vs x for various N
% figure(1);
% hold on;
% for i = 1:length(n)
%     plot(x, sinSeries_array(:, i), 'DisplayName', ['n = ', num2str(n(i))]);
% end
% hold off;
% 
% % Add labels to plots
% xlabel('x');
% ylabel('sin(x) Approximation');
% title('sinSeries Approximation vs x for Various N');
% legend('show');
% grid on;

% Plots error vs x for various n
figure(2);
hold on;
for i = 1:length(n)
    plot(x, error(:, i), 'DisplayName', ['n = ', num2str(n(i))]);
end
hold off;

% Make plot labels
xlabel('x');
ylabel('Error');
title('Error of sinSeries vs X');
legend('show');
grid on;
% 
% % Plot L2-norm vs N
% figure(3);
% plot(n, l2_error_overn);
% xlabel('N');
% ylabel('L2-Norm of Error over N');
% title('L2-Norm of Error (N) vs N');
% grid on;
% 
% % Plot L2-norm vs N
% figure(4);
% plot(x, l2_error_overx);
% xlabel('X');
% ylabel('L2-Norm of Error over X');
% title('L2-Norm of Error (X) vs X');
% grid on;
% 
% % Plot log10(l_2-error)
% figure(5);
% plot(n, log10(l2_error_overn));
% xlabel('N');
% ylabel('log10(L2-Norm (N))');
% title('log10(L2-Norm (N)) vs N');
% grid on;
% 
% % Plot log10(l_2-error)
% figure(6);
% plot(x, log10(l2_error_overx));
% xlabel('X');
% ylabel('log10(L2-Norm (X))');
% title('log10(L2-Norm (X)) vs X');
% grid on;
% 
% 
% % Plot log10(l_2-error) vs log10(n)
figure(6);
plot(log10(n), log10(l2_error_overn));
xlabel('log10(N)');
ylabel('log10(L2-Norm)');
title('log10(L2-Norm) vs log10(N)');
grid on;

% Plot log10(l_2-error) vs log10(n)
figure(7);
plot(log10(x), log10(l2_error_overx));
xlabel('log10(x)');
ylabel('log10(L2-Norm (X))');
title('log10(L2-Norm (X)) vs log10(X)');
grid on;


%Plot error vs n for different series of x
figure(8);
hold on;
for i = 1:length(x)
    plot(n, error(i, :), 'DisplayName', ['x = ', num2str(x(i))]);
end
hold off;
xlabel('N');
ylabel('Error');
title('Error vs N');
legend('show');
grid on;

%%% Meanings for graphs displayed
%
% Figure 2:
% Error of sinSeries vs X
%
% This graph displays the relationship for the numerical error between the
% approximated and true value of sin(x). Visually, we can see that the
% magnitude of error gets larger as the value of x increases. While
% slightly harder to see due to the multitude of series, the values with
% the highest error at constant values of x correspond to lower values of
% n. Both relationships will be more clearly shown in other graphs. 
%
% Figure 6:
% $log10(L2-Norm (N))$ vs $log10(N)$
%
% This graph shows a logarithmically scaled plot of L2-normalization of the error over N vs values of N.
% It shows that as N increases, the error of the function decreases. 
%
% Figure 7:
% $log10(L2-Norm (X))$ vs $log10(X)$
%
% This plot shows the increasing error in approximation as x increases. It
% is almost a linear graph and represents a direct relationship.
%
% Figure 8:
% Error vs N for series of x
%
% This figure shows the most intuitive and understandable patterns of the
% data. It is obvious how the error (calculated directly from subtracting
% sin(x) from the approximation of sin(x))decreases as the value of N
% increases. Similarly, we can see a pattern in the series of x, as x
% increases the error of the approximation of sin(x) increases as well.
% From the graph we can also see that the increase in error is regular and
% increases at specific intervals for certain values of x. 

%%% SinSeries Approximation Function
%
% Here, we compute the taylor series approximation using the summation
% formula. It takes the inputs from x and n to calculate the approximation.
% 
function [approx] = sinSeries(x, N)
    approx = single(0);
    for k = 0:N
        approx = approx + (((-1).^k)/factorial(2.*k+1)).*(x.^(2.*k +1));
    end
end

%%% Conclusion
%
% From the data, we can see that there exists certain ranges of x and N
% values for which the error in approximation using taylor series
% expansions increases dramatically. In this case, we can see a significant
% increase in error for small values of N and large value of X. Related to
% N, this comes from the general approximation of a taylor series being
% more accurate at higher values. Since it is meant to extent to infinity
% to adequately match the approximated function, smaller values do not have
% enough terms to "chip away" at the approximation. At small values
% of X, this is not very noticable since small values of x converge much
% faster than large values of X. For larger values of X, they take much
% longer to converge, thus is N is not sufficiently large, the series will
% be truncated before the approximation is adequate. 
%
%% Problem 2: Heaviside Calculus

% <include>prob2.m</include>

