function [F,J] = eigen_val_compute(x)

syms y

A =[4,2,1;2,3,0;1,0,1];
%my_func = @(y) y.^2;

k = ones(size(A,1),1);
my_func = @(y) vertcat(A*y(1:end-1),(y(1:end-1)'.^2)*k) - [y(end)*y(1:end-1);1];

F = my_func(x);

J_1 = @(y) horzcat(vertcat(A -y(end)*eye(size(A)),2*y(1:end-1)'),vertcat(-y(1:end-1),0));
J = J_1(x);

end