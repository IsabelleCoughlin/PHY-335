classdef LegendreScore < NumericalRecipes.FunctorShoot
    properties
        even
    end
    methods
        function obj = LegendreScore(even)
            obj.even = even;
        end
        function y = vector(obj,x,v)
            % v(1) = P, v(2)=Q, v(3)=lambda
            % y(1,1) = P or Q
            if obj.even
                y(1,1) = v(2);
            else
                y(1,1) = v(1);
            end
        end
    end
end