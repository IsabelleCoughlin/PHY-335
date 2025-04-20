%% Problem 2: Heaviside Calculus
% Use Heaviside calculus to obtain a formula of maximum order for
% $\int_{0}^{2h}f(x)dx$ in terms of $f(0), f(h), f'(h), f''(h),$ 
% and $f(2h)$. Compare your formula with Simpson’s rule, Numerical Recipes (4.1.4). 
% If it is different, explain why.
%
% This is equivalent to performing the 3-point scheme integration method.
% Our goal is to generate an approximation to the integral using a linear
% combination of the values given, primarily $f(0), f(h), f'(h), f''(h),$ 
% and $f(2h)$. We write this as:
%
% $\int_{0}^{2h}f(x)dx = Af(x-h)+Bf(x)+Cf(x+h)+Err[f(x)]$
%
% In this case, our x-variable is h, so that x-h = 0, x = h, and x + h =
% 2h.
%
% First, we define the heaviside operator, $\hat{D}$, which is a derivative
% operator such that
%
% $\hat{D}f(x) = f'(x)$
%
% Thus, it's inverse represents integration as
%
% $\frac{1}{\hat{D}}f(x) = \int f(x)$
%
% Yet, this represents in indefinite integral. In order to apply bounds, we
% must apply another operator to create the definite integral from $x=0$ to
% $x=2h$. In order to apply this, we must use the underlying principle of 
% using the taylor series expansions to define what $f(x+h)$ and $f(x-h)$ are. 
%
% $f(x+h) = f(x) + hf'(x) + \frac{1}{2}h^2f''(x)+...$
%
% $f(x-h)= f(x) - hf'(x) + \frac{1}{2}h^2f''(x) + ...$
%
% Acknowledging that $f'(x) = \hat{D}f(x)$, $f''(x) = {\hat{D}}^{2}f(x)$,
% we can write the same equations as
%
% $f(x+h) = f(x) + h\hat{D}f(x) + \frac{1}{2}h.^2{\hat{D}}^{2}f(x) + ...$ 
% 
% $f(x+h) = f(x)[1 + h\hat{D} + \frac{1}{2}h.^2{\hat{D}}^{2} + ... $
%
% $f(x-h) = f(x) - h\hat{D}f(x) + \frac{1}{2}h.^2f{\hat{D}}^{2}f(x) - ... $
% 
% $f(x-h) = f(x)[1 - h\hat{D} + \frac{1}{2}h.^2{\hat{D}}^{2} - ... $
%
% Respectively, these equations correspond to the taylor series expansions of 
% for e^{h\hat{D}}*f(x), e^{-h\hat{D}}*f(x). As such, we can rewrite our
% original equation using the replacements above. 
% 
% $\int_{0}^{2h}f(x)dx = Ae^{-h\hat{D}}*f(x)+Bf(x)+Ce^{h\hat{D}}*f(x)+Err[f(x)]$
%
% Now, we can use the same taylor series approximations to add a "shift
% operator" to the integral. When solving a definite integral, we first
% solve the indefinite integral. Then we evaluate the integral at the end
% bound, and then subtract the integral evaluated at the start bound. Thus,
% for this equation, we subtract f(x+h) from f(x-h), or 
%
% $(e^{h\hat{D}} - e^{-h\hat{D}})\frac{1}{\hat{D}}f(x) =
% \int_{0}^{2h}f(x)dx$

