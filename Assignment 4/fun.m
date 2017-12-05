function [F,G,H] = fun(x)

syms y1 y2
f = @(y) 10*(y2-y1^2)^2 +(y1-1)^2;

F = double(subs(f,[y1;y2],x));
grad = gradient(f,[y1;y2]);
G = double(subs(grad,[y1;y2],x));
hess = hessian(f,[y1;y2]);
H = double(subs(hess,[y1;y2],x));

end