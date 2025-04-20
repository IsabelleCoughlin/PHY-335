%Functor for a x^3 + b x^2 + c^x +d
classdef TestFunctionD < NumericalRecipes.FunctorD
    properties
        a
        b
        c
        d
    end
    methods
        function obj = TestFunctionD(a,b,c,d)
            obj = obj@NumericalRecipes.FunctorD();
            obj.a = a;
            obj.b = b;
            obj.c = c;
            obj.d = d;
        end
        function val = func(obj,x)
            val = ((obj.a.*x + obj.b).*x + obj.c).*x +obj.d;
        end
        function val = df(obj,x)
            val = (3.0*obj.a.*x + 2.0*obj.b).*x + obj.c;
        end
    end
end