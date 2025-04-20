classdef HW7_Functor < NumericalRecipes.FunctorODE
    properties
        q
    end
    methods
        function obj = HW7_Functor(q)
            %obj@NumericalRecipes.FunctorODE();
            obj.q = q;
        end
        function derivs = dydx(obj, x, y)
            derivs(1) = y(2);
            derivs(2) = -(y(3) - 2*obj.q*cos(2*x))*y(1);
            derivs(3) = 0;
        end
    end
end
