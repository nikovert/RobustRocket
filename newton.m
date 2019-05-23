function x = newton(f, df, x0, atol, rtol, nmax)
%
% NEWTON Newton's Method
%   Newton's method for finding successively better approximations to the 
%   zeroes of a real-valued function.
%
% Input:
%   f - input funtion
%   df - derived input function
%   x0 - inicial aproximation
%   tol - tolerance
%   nmax - maximum number of iterations
%
% Output:
%   x - aproximation to root
%
% Example:
%	x = newton( @(x) exp(x)+x, @(x) exp(x)+1, 0, 1e-4, 1e-4, 10 )
%
    if nargin == 3
        atol = 1e-7;
        rtol = 1e-7;
        nmax = 1000;
    elseif nargin == 4
        rtol = 1e-7;
        nmax = 1000;
    elseif nargin == 5
        nmax = 1000;
    elseif nargin ~= 6
        error('newton: invalid input parameters');
    end

    n = 1;
    x = x0;
    dx = f(x)/df(x);

    while n < nmax && abs(dx) > atol + rtol*abs(x)
        x = x - dx;
        n = n + 1;
        dx = f(x)/df(x);
    end
end