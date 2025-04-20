classdef MathieuLoad_b < NumericalRecipes.FunctorShoot
    methods (Static)
        function y = vector(x,v)
            y(1) = 1;         % P
            y(2) = 0;    % Q
            y(3) = v(1);        % lambda
        end
    end
end