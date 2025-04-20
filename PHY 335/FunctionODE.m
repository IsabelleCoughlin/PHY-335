classdef FunctionODE < NumericalRecipes.FunctorODE
    properties
        count = 0;
    end
    methods
        function obj = FunctionODE(count)
            obj = obj@NumericalRecipes.FunctorODE();
            obj.count = count;
        end
        function derivs = dydx(obj,x,y)
            obj.count = obj.count + 1;
            derivs(1) = y(2)/x; 
            derivs(2) = -(x-(1/x))*y(1);
        end
    end
end