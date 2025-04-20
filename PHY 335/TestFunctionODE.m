classdef TestFunctionODE < NumericalRecipes.FunctorODE
    properties
        count = 0;
    end
    methods
        function obj = TestFunctionODE(count)
            obj = obj@NumericalRecipes.FunctorODE();
            obj.count = count;
        end
        function derivs = dydx(obj,x,y)
            obj.count = obj.count + 1;
            derivs(1) = y(2);
            derivs(2) = -(1/x)*y(2) - (1- (1/x^2))*y(1);
            %derivs(2) = -(1/x)*y(2) - (1 - (1/x^2)) * y(1);
   
        end
    end
end