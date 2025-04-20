x = pi/4:pi/2:4*pi; % Range of x values
n = 2:15; % Range of N values
sinSeries_array = zeros(length(x), length(n)); % Store approximations
error = zeros(length(x), length(n)); % Store errors

% Compute approximations and errors
for x_val = 1:length(x)
    for n_val = 1:length(n)
        sinSeries_array(x_val, n_val) = sinSeries(x(x_val), n(n_val)); % Approximation
        error(x_val, n_val) = sinSeries_array(x_val, n_val) - sin(x(x_val)); % Error
    end
end

% Plot approximations vs x for various N
figure(3);
hold on;
for i = 1:length(n)
    plot(x, sinSeries_array(:, i), 'DisplayName', ['N = ', num2str(n(i))]);
end
plot(x, sin(x), 'k--', 'DisplayName', 'True sin(x)'); % True sine function
hold off;
xlabel('x');
ylabel('sin(x) Approximation');
title('sinSeries Approximation vs x for Various N');
legend('show');
grid on;

% Plot error vs x for various N
figure(1);
hold on;
for i = 1:length(n)
    plot(x, error(:, i), 'DisplayName', ['N = ', num2str(n(i))]);
end
hold off;
xlabel('x');
ylabel('Error');
title('Error of sinSeries vs sin(x)');
legend('show');
grid on;

% Compute L2-norm of the error for each N
L2_norm = sqrt(sum(error.^2, 1)); % Sum over x dimension

% Plot L2-norm vs N
figure(4);
plot(n, L2_norm, '-o');
xlabel('N');
ylabel('L2-Norm of Error');
title('L2-Norm of Error vs N');
grid on;

% Plot error vs N for different x values
figure(2);
hold on;
for i = 1:length(x)
    plot(n, error(i, :), 'DisplayName', ['x = ', num2str(x(i))]);
end
hold off;
xlabel('N');
ylabel('Error');
title('Error vs N for Various x');
legend('show');
grid on;

% Function to compute the truncated series expansion of sin(x)
function [approx] = sinSeries(x, N)
    approx = 0;
    for k = 0:N
        approx = approx + ((-1)^k / factorial(2*k + 1)) * (x^(2*k + 1));
    end
end
