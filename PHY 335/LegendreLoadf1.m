classdef LegendreLoadf1 < NumericalRecipes.FunctorShoot
    methods(Static)
        function y = vector(x,v)
            % v(1) =  P_1, v(2)=lambda1, v(3)=lambda2
            % y(1) = P, y(2)=Q, y(3)=lambda
            y(1) = v(1);                % P (parameter)
            y(2) = -v(2)*v(1)/2.0;      % Q (boundary condition)
            y(3) = v(2);                % lambda (parameter)
        end
    end
end