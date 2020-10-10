function [xk, res, it] = secant(a, b, tol, itmax, fun)
    
    % Initial conditions
    xk = [a, b];
    it = 0;
    res = [abs(b - a)];

    % Tolerance to iterate over different successions
    tolk = res(1);

    while it < itmax & tolk > tol
        fk = xk(end) - fun(xk(end))*(xk(end) - xk(end - 1))/(fun(xk(end)) - fun(xk(end - 1)));
        tolk = abs(xk(end) - fk);
        xk = [xk, fk];
        res = [res, tolk];
        it = it + 1;
    end

end
