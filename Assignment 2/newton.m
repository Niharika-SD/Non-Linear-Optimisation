function [x,F,J,iter,status] = newton(Fun,x0,maxit,printlevel,tol)

x = x0;
count = 1;
status =0;
fh = str2func(Fun);

for i = 1:maxit
    
    [F_iter,J_iter] = fh(x);
    
    if (count < printlevel)
        fprintf('\n');
        fprintf('Iteration %d || Function Value : %f', i,norm(F_iter,'fro'));
    end
    
    if (i==1)
        F_x0 = F_iter;
    end
    
    x = x - J_iter\F_iter;
    
    if((norm(F_iter,'fro')/norm(F_x0,'fro')) < tol)
        status =1;
        break;
    end
   
    count = count + 1;
    
end

F = F_iter;
J = J_iter;
iter = i;
end