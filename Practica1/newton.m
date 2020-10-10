function [xk, res, it] = newton(x0, tol, itmax, fun, dfun)
    
    fNewton = @(x) x - fun(x)/dfun(x);
    [xk, res, it] = iteracio_simple(x0, tol, itmax, fNewton);

end
