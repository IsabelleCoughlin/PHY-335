classdef MathieuScore_d < NumericalRecipes.FunctorShoot
    methods (Static)
        function y = vector(x,v)
            %Which one is unique
            y(1,1) = v(1)
        end
    end
end