function [xk, res, it] = newton(x0, tol, itmax, fun, dfun)

    % Initialization of the output parameters
    it = 0;
    tolk = abs(fun(x0)/dfun(x0));
    xk = [x0];
    res = [fun(x0)];
    
    % Loop for computing consecutive iterates
    while tolk > tol
        if it == itmax
            printf("Warning: Reached maximum number of iterations, itmax = %d. Stopping iterative method", itmax);
            break;
        end
        
        % Computing the next step
        fk = xk(end) - fun(xk(end))/dfun(xk(end));
        tolk = abs(xk(end) - fk);
        xk = [xk, fk];
        res = [res, fun(fk)];
        it = it + 1;
    end 
end
