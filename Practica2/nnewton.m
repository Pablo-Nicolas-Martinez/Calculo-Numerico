function [XK, DFk, res, it] = nnewton(x0, tol, itmax, fun, dfun)
   
    % Preliminary declarations
    XK = x0;
    res = norm(fun(x0));
    it = 0;
    tolk = norm(dfun(x0)\x0);
    
    % Loop for iterative search of solutions
    while tolk > tol
        if it == itmax
            fprintf("Maximum number of iterations achieved. Stopping iterative procedure.\n");
            break;
        end
        dk = dfun(XK(:,end))\fun(XK(:,end));
        tolk = norm(dk);
        XK = [XK, XK(:, end) - dk];
        res = [res, norm(fun(XK(:, end)))];
        it = it + 1;
    end
    
    DFk = dfun(XK(:, end));
   
end