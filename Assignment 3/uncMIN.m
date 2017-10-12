function [x,F,G,H,iter,status] = uncMIN(fun,x0,step,maxit,printlevel,tol)

x = x0;
count = 1;
status =0;
fh = str2func(fun);
tau = 0.5;
nu = 0.7;
alpha_init =1;

for i = 1:maxit
    
    [F_iter,G_iter,H_iter] = fh(x);
    alpha_iter= alpha_init;
  
    if (i==1)
        F_x0 = F_iter;
    end
    
    if((norm(F_iter,'fro')/norm(F_x0,'fro')) < tol)
        status =1;
        break;
    end
    
    if (count < printlevel)
        fprintf('\n');
        fprintf('Iteration %d || Function Value : %f \n', i,norm(F_iter,'fro'));
    end
    
  
    
    if (step == 0)
        p_iter = -G_iter;
        alpha_iter =0.01;
    else
        [B_iter,~] = modNewton(H_iter,10);
        p_iter = -B_iter\G_iter;
        
        for l = 1:100
        
        [F_iter_l,~,~] = fh(x+ alpha_iter*p_iter);
        
        if (F_iter_l > F_iter+nu*alpha_iter*G_iter'*p_iter)
        alpha_iter = alpha_iter*tau;
        end
        
        end
    end
    
    x = x + alpha_iter*p_iter;
    count = count + 1;
    
end

iter = i;
F = F_iter;
G= G_iter;
H= H_iter;

end