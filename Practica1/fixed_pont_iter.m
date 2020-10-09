% Fixed point iteration. x = f(x).
% Input:   1. x0: initial guess
%          2. tol: tolerance so that abs(x_k+1 - x_k) < tol
%          3. itmax: maximum number of iterations allowed
%          4. fun: function''s name 
% Output:  1. xk: resulting sequence 
%          2. res: resulting residuals 
%          3. it: number of required iterations  

function [xk,res,it] = fixed_point_iter(x0,tol,itmax,fun)

  xk = [x0] ; it = 0 ; res=[abs(x0-fun(x0))];

  tolk=abs(res(1));

  while it < itmax & tolk > tol
      fk = fun(xk(end)) ; tolk = abs(xk(end)-fk) ;
      xk = [xk  fk];
      res = [res tolk] ;
      it = it + 1 ;
  end

