%Functor for a x^3 + b x^2 + c^x +d
classdef hwfivefunctor < NumericalRecipes.FunctorD
    properties
        p
    end
    methods
        function obj = hwfivefunctor(p)
            obj = obj@NumericalRecipes.FunctorD();
            obj.p = p;
        end
        function val = func(obj,x)
            val = cos(x)-0.8+obj.p*x.^2;
        end
        function val = df(obj,x)
            val = -sin(x) + 2.0*obj.p*x;
        end
    end
end