%% Problem 2: Heaviside Calculus

% The second approximations for the numerical derivative of a function f(x) is found using a 2-point centered difference.

% $$\frac{df(x)}{dx} = Af(x + \frac{h}{2}) + Bf(x - \frac{h}{2}) + Err[f(x)]$$

% We apply the $\hat{D}$ operator to the equation.

% $$\hat{D}[\frac{df(x)}{dx}] = Ae^{\frac{1}{2}h\hat{D}}f(x) + Be^{-\frac{1}{2}h\hat{D}}f(x) + Err[f(x)]$$

% Here, we can divide both sides of the equation by f(x). Also, we
% recognize the $e^{\frac{1}{2}h\hat{D}}$ and $e^{-\frac{1}{2}h\hat{D}}$
% terms with Taylor-Series expansions, as shown below.

% $$e^{\frac{1}{2}h\hat{D}} = 1 + \frac{1}{2}h\hat{D} +
% \frac{1}{2}(\frac{1}{2}h)^{2}(\hat{D})^{2} +
% \frac{1}{6}(\frac{1}{2}h)^{3}(\hat{D})^{3} + ...

% $$e^{-\frac{1}{2}h\hat{D}} = 1 - \frac{1}{2}h\hat{D} +
% \frac{1}{2}(\frac{1}{2}h)^{2}(\hat{D})^{2} -
% \frac{1}{6}(\frac{1}{2}h)^{3}(\hat{D})^{3} + ...

% With this equation, we can solve for A and B. 

% Since we have a $\hat{D}$ term on the left hand side of the equation, we
% must solve for A and B such that the D^0 term is zero, and the
% D^1 term is 1.

%%% Solving for D^0

% $$ A + B = 0
% $$ A = -B

%%% Solving for D^1

% $$ (A-B)\frac{h}{2} = 1
% $$ (A-(-A))\frac{h}{2} = 1
% $$ 2A\frac{h}{2} = 1
% $$ A = \frac{1}{h}
% $$ B = -\frac{1}{h}

%%% Error Factor 

% The next term in the expansion relates to the leading error term.
% We can exlude the factors used in solving for coefficients. 

% To find the error estimate, we plug coefficients back into the original
% function, other than the first two terms we used to solve A,B. 

% $$ \hat{D} = (\frac{1}{h})(\frac{1}{2}(\frac{1}{2}h)^{2}(\hat{D})^{2} + \frac{1}{6}(\frac{1}{2}h)^{3}(\hat{D})^{3}) + (-\frac{1}{h})(
% \frac{1}{2}(\frac{1}{2}h)^{2}(\hat{D})^{2} - \frac{1}{6}(\frac{1}{2}h)^{3}(\hat{D})^{3}) + Err

% $$ \hat{D} = \frac{1}{8}(h)(\hat{D})^{2} + \frac{1}{48}(h)^{2}(\hat{D})^{3}) - 
% \frac{1}{8}h(\hat{D})^{2} + \frac{1}{48}(h)^{2}(\hat{D})^{3}) + Err

% $$ \hat{D} = \frac{1}{48}(h)^{2}(\hat{D})^{3}) + \frac{1}{48}(h)^{2}(\hat{D})^{3}) + Err

% Since the second order factor is zero, the leading error term of this
% function is determined by the $O(h^{4})$ factor. 

% $$ \hat{D} = \frac{1}{24}(h)^{2}(\hat{D})^{3}) + Err

% We find the leading order in the truncation error to be:

% $$ Error = -\frac{1}{24}(h)^{2}(f(x))^{3})

%%% Conclusion

% Thus, the derivative can be approximated by

% f(x)^{1} = ((\frac{1}{h})(f(x + \frac{1}{2}h) - f(x - \frac{1}{2}h)))
% -\frac{1}{24}(h)^{2}(f(x))^{3}) + O(h^{4})




