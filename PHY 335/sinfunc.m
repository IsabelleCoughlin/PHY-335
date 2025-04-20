classdef sinfunc < NumericalRecipes.Functor
    methods
        % Constructor
        function obj = sinfunc()
            % e^(-x*siny)
            obj = obj@NumericalRecipes.Functor();
        end
        function val = func(obj, x)
            val = exp(sin(x));
        end
    end
end