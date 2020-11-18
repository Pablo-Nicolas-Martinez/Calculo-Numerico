function X = contzvc(x0, ds, smax, sgn, f, Df)
    
    % Preliminary declarations: finding a point in the curve
    %[XK, ~, ~, ~] = nnewton(x0, 1e-12, 100, f, Df);
    %X = XK(:, end);
    
    X = x0;
    
    % Choosing an appropriate initial submatrix and vector
    n = length(x0); % Not consistent with the problem definition
    Df0 = Df(x0);
    Dtmp = Df0(:, 2:end);
    for i = 1:n
        if (det(Dtmp) < 1e-15)
            Dtmp(:, i) = Df0(:, i);
        else
            vkp = Dtmp\(-Df0(:, i));
            vk = [vkp(1:i - 1); 1; vkp(i:end)];
            break
        end
    end
    
    % Normalization of the vector and final definition
    vk = sgn*vk/norm(vk);
    
    % Loop to iterate on points of the curve
    for i = 0:ds:smax
        
        % Prediction vector
        vk = [Df(X(:, end)); vk']\[zeros(n - 1, 1); 1];
        
        % Correction to the estimate
        fun = @(x) [f(x); dot(x - X(:, end), vk) - ds];
        dfun = @(x) [Df(x); vk'];
        [Xnew, ~, ~, ~] = nnewton(X(:, end) + ds*vk, 1e-12, 100, fun, dfun);
        
        X = [X, Xnew(:, end)];
    end
    
end