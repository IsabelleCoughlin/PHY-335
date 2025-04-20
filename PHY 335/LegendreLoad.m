classdef LegendreLoad < NumericalRecipes.FunctorShoot
    methods(Static)
        function y = vector(x,v)
            % v(1) =  lambda
            % y(1) = P, y(2)=Q, y(3)=lambda
            y(1) = 1.0;         % P
            y(2) = v(1)/2.0;    % Q
            y(3) = v(1);        % lambda
        end
    end
end