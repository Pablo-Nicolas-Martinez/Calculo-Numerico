function [xk, res, it] = iteracio_simple(x0, tol, itmax, fun)
    
    % Initialization of the output parameters
    it = 0;
    tolk = abs(fun(x0) - x0);
    xk = [x0];
    res = [];
    
    % Loop for computing consecutive iterates
    while tolk > tol
        if it == itmax
            printf("Warning: Reached maximum number of iterations, itmax = \n. Stopping iterative method", itmax);
            break;
        end
        
        % Computing the next step
        fk = fun(xk(end));
        tolk = abs(xk(end) - fk);
        xk = [xk, fk];
        res = [res, tolk];
        it = it + 1;
    end 
end
