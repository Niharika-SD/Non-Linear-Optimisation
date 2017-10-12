function [F,G,H] = fun(x)
syms y1 y2
my_func = 10*(y2-y1^2)^2 +(y1-1)^2;

grad = gradient(my_func,[y1;y2]);
hess = hessian(my_func,[y1;y2]);
F = double(subs(my_func,[y1;y2],x));
G = double(subs(grad,[y1;y2],x));
H = double(subs(hess,[y1;y2],x));
end