function [p, iters, flag] = steihaug_CG(B,g,radius,tol)

if (norm(B-B')~=0)
    fprintf('Matrix B is non symmetric: error ')
end

p_0 = zeros(size(g));
r_0 = g;
s_0 = -g;

for k = 1:1000
    
    if (k==1)
        p_k =p_0;
        r_k=r_0;
        s_k = s_0;
    end
    
    if(norm(r_k,2)> tol* norm(r_0,2))
    
        if (s_k'*B*s_k>0)
            alpha_k = (r_k'*r_k)/(s_k'*B*s_k);
        else
            tau = sqrt((radius^2 - norm(p_k,2)^2)/norm(s_k,2)^2); 
            p_k = p_k +tau*s_k;
            p = p_k;
            iters =k;
            flag =-1;
            break;
        end
    
        if (norm(p_k+alpha_k*s_k,2)< radius)
           p_k = p_k +alpha_k*s_k;
        else 
           tau = sqrt((radius^2 - norm(p_k,2)^2)/norm(s_k,2)^2); 
           p_k = p_k +tau*s_k;
           p = p_k;
           iters =k;
           flag = 1;
           break;
        end
        
        r_k_1 = r_k + alpha_k*B*s_k;
        Beta = (r_k_1'*r_k_1)/(r_k'*r_k) ;
        s_k =  -r_k_1 + Beta*s_k;
        r_k = r_k_1;
    
    else
        
        iters = k;
        p = p_k;
        flag =0;
        break
    
    end
end
end
