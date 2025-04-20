classdef LegendreScoref < NumericalRecipes.FunctorShoot
    methods(Static)
        function y = vector(x,v)
            % v(1) = P, v(2)=Q, v(3)=lambda
            % y(1,1) = P, y(2,1)=Q, y(3,1)=lambda
            % y must be a column vector!
            y(1,1) = v(1);
            y(2,1) = v(2);
            y(3,1) = v(3);
        end
    end
end