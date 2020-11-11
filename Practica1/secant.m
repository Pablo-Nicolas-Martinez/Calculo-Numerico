function [xk, res, it] = secant(a, b, tol, itmax, fun)
    
    % Initial conditions
    xk = [a, b];
    it = 0;
    res = [fun(a), fun(b)];

    % Tolerance to iterate over different successions
    tolk = tol + 1;

    while tolk > tol
        if it == itmax
            printf("Warning: Reached maximum number of iterations, itmax = \n. Stopping iterative method", itmax);
            break;
        end
        fk = xk(end) - fun(xk(end))*(xk(end) - xk(end - 1))/(fun(xk(end)) - fun(xk(end - 1)));
        tolk = abs(xk(end) - fk);
        xk = [xk, fk];
        res = [res, fun(fk)];
        it = it + 1;
    end

end
