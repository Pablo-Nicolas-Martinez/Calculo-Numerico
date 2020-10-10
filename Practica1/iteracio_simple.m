function [xk, res, it] = iteracio_simple(x0, tol, itmax, fun)
    
    % Initialization of the output parameters
    it = 0;
    res = [abs(fun(x0) - x0)];
    xk = [x0];
    
    tolk = res(1);
    
    % Loop for computing consecutive iterates
    while tolk > tol & it< itmax
        % Computing the next step
        fk = fun(xk(end));
        tolk = abs(xk(end) - fk);
        xk = [xk, fk];
        res = [res, tolk];
        it = it + 1;
    end
end
