classdef MathieuLoad < NumericalRecipes.FunctorShoot
    methods (Static)
        function y = vector(x,v)
            y(1) = 1;        
            y(2) = 0;    
            y(3) = v(1);        
        end
    end
end