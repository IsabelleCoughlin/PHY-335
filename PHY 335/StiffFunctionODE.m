classdef StiffFunctionODE < NumericalRecipes.FunctorODE
    methods
        function obj = StiffFunctionODE()
            obj = obj@NumericalRecipes.FunctorODE();
        end
        function derivs = dydx(obj,x,y)
            derivs(1) = 998*y(1)+1998*y(2);
            derivs(2) = -999*y(1)-1999*y(2);
        end
    end
end