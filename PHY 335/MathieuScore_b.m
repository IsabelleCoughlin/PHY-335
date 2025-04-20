classdef MathieuScore_b < NumericalRecipes.FunctorShoot
    methods (Static)
        function y = vector(x,v)
            %Which one is unique
            y(1,1) = v(2);
        end
    end
end