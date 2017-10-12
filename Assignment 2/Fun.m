function [F,J] = Fun(x)

syms y

my_func = @(y) (y-3).^2;

F = my_func(x);

J_1 = jacobian(my_func,y);
J = subs(J_1,y,x);

end