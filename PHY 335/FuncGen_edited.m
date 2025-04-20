% FuncGen is a simple example function that can be used to
% generate the x and y values representing some function.
% This function is written as an example of a function that
% can be called by any other function or script within MATLAB
% so long as this file is within the MATLAB search path.
% Because the name of the file is the same as the name of the
% function "FuncGen".m
%
% This function returns 2 arrays.  The first is the values of 
% x at which the function is evaluated.  The second is the 
% values y=f(x) at each x value.
% 
% The input values a, b, and n are used to determing the
% values of x at which the function will be evaluated.
% The input value p is an array that is used if the function
% f(x) is chosen to be a polynomial.  See the MATLAB polyval
% for more information on this array.  The input value err
% is used to add random error to the returned function
% y = f(x) + error(x)
%
function [x,y] = FuncGen(a,b,n,p,err)
    % Set up the values of x
    x = a:(b-a)/(n-1):b;
    %
    % Choose the function to be returned
    if isequal(p,[])
        y = (1+(25*(x.^2))).^(-1);
    else
        y = polyval(p,x);
    end
%     y = sin(x);
%     y = exp(-4*x.^2);
%     y = exp(-2*abs(x));
    %
    % Create a list of random values [one for each (x,y) point]
    % with a value between -err..err.
    %del = (2*rand(1,length(y))-1)*err;
    %del = random('Uniform',-err,err,1,length(y));
    %
    % Add the random error to each function value;
    %y = y + del;
end