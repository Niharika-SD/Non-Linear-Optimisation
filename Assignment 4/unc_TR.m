function [x,F,G,H,iter,status] = unc_TR(fun,x0,maxit,printlevel,tol)

x_k = x0;
eta_vs = 0.9;
eta_s =0.1;
gam_d = 0.5;
gam_i =2; 
status = 1;
radius = 0.5; 

for iter = 1: maxit
    
    fh = str2func(fun);
    [F_k,G_k,H_k]  = fh(x_k);
    
    if (iter ==1)
        F_0 = F_k;
        G_0 =G_k;
    end
    
    if (printlevel)
        fprintf('\n Iter : % d || function value : %f ',iter, double(norm(F_k,'fro')/norm(F_0,'fro')))
    end
    
    % trial step using steihaug_CG and second derivative information
    [s_k, ~, ~] = steihaug_CG(H_k,G_k,radius,tol);

    [F_k_s,~,~] = fh(x_k+s_k);
    
    rho_k = -(F_k-F_k_s)/(G_k'*s_k+0.5*s_k'*H_k*s_k);
    
    if (rho_k>=eta_vs)
        % very successful iteration
        x_k =x_k +s_k;
        radius = gam_i*radius;
    
    elseif(rho_k>=eta_s)
        %  successful iteration
        
        x_k =x_k +s_k;
        
    else
       % unsuccessful iteration
         
        radius = gam_d*radius;
        
    end
        
    if(norm(G_k,2) <= 10e-06*max(1 ,norm(G_0,2)))
        x = x_k;
        F = F_k;
        G = G_k;
        H = H_k;
        status = 0; 
        break;
    end
end

end