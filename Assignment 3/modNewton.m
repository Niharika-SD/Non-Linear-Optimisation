function [B, flag] = modNewton(H,beta)

flag = 0; 
if (norm(H-H','fro')~=0)
    fprintf('The given hessian is not symmetric: error \n')
    return;
end

if (beta<=1)
    fprintf(' The value of beta needs to be greater than 1: error \n')
    return;
end


if (norm(H,'fro')==0)
    epsi = 1;
else
    epsi = norm(H,2)/beta;
    [V,D] = eig(H);
end

D_bar = max(D,epsi*eye(size(H)));

if norm(D_bar-D,'fro') ~=0
    flag =1;
    B = V*D_bar*V';
else
    B = H;
end


end