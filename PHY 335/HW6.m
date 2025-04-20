%% Homework 6
%
%
x1 = 1 % Starting value
x2 = 30

J_0 = besselj(0,1)
J_1 = besselj(1,x1)
J_1_prime = J_0 - J_1/x1

J_1_0 = besselj(1,0)


J1 = besselj(1,1)  % Should be about 0.4401
J1_prime = (besselj(0,1) - besselj(1,1)/x1) % Should be close to 0.3251


