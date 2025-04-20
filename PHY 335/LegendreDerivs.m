classdef LegendreDerivs < NumericalRecipes.FunctorODE
    methods(Static)
        function derivs = dydx(x,y)
            % y(1) = P, y(2)=Q, y(3)=lambda
            derivs(1) = y(2);
            derivs(2) = (2.0*x*y(2) - y(3)*y(1))/(1.0-x*x);
            derivs(3) = 0.0;
        end
    end
end