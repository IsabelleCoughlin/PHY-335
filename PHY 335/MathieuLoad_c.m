classdef MathieuLoad_c < NumericalRecipes.FunctorShoot
    methods (Static)
        function y = vector(x,v)
            y(1) = 0;         % P
            y(2) = 1;    % Q
            y(3) = v(1);        % lambda
        end
    end
end