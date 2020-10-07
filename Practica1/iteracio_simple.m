function [xk, res, it] = iteracio_simple(x0, tol, itmax, f)
    
    % Initialization of the output parameters
    it = 0;
    res = [f(x0)];
    xk = [x0];
    
    current = x0;
    error = f(x0);
    
    % Loop for computing consecutive iterates
    while ( abs(current - error) > tol) 
        
        % Checking the number of iterations
        if (it > itmax)
            printf("El mètode no ha convergit en prou temps");
            break;
        end
        
        % Computing the next step
        current = f(current);
        error = f(error);
        xk = [xk; current];
        res = [res; error];
    end
end