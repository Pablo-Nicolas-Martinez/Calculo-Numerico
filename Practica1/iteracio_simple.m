function [xk, res, it] = iteracio_simple(x0, tol, itmax, fun)
    
    % Initialization of the output parameters
    it = 0;
    xk = [x0];
    res = [fun(x0) - x0];
    tolk = abs(res(1));
    
    % Loop for computing consecutive iterates
    while tolk > tol
        if it == itmax
            printf("Warning: Reached maximum number of iterations, itmax = %d. Stopping iterative method", itmax);
            break;
        end
        
        % Computing the next step
        fk = fun(xk(end));
        tolk = abs(xk(end) - fk);
        xk = [xk, fk];
        res = [res, fun(xk(end)) - xk(end)];
        it = it + 1;
    end 
end
