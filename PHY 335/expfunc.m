classdef expfunc < NumericalRecipes.Functor
    methods
        % Constructor
        function obj = expfunc()
            % e^(-x*siny)
            obj = obj@NumericalRecipes.Functor();
        end
        function val = func(obj, x)
                % Write the bounds
            
            lowerbound = (3*x - (-16*x^2 + 10)^(1/2))/5
            upperbound = (3*x + (-16*x^2 + 10)^(1/2))/5
            p = NumericalRecipes.Midpnt(@(y) exp(-x*sin(y)), lowerbound, upperbound)
            val = NumericalRecipes.qromo(p, 10^-5)
        end
    end
end