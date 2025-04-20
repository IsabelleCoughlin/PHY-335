classdef MathieuScore < NumericalRecipes.FunctorShoot
    methods
        function y = vector(obj,x,v)
            %Which one is unique
            y(1,1) = 1-v(1);
        end
    end
end