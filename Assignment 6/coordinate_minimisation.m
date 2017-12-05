 function [x,iter] = coordinate_minimisation(H,type)

n= size(H,1);
x0 = 0.001*rand(n,1);
%x0 = ones(n,1);

alpha = 1/normest(H,2);
for iter = 0:999
        
        if (iter == 0)
            x = x0;
        end
        
        fprintf(' Iteration %d || Function Value: %f \n',iter+1,(0.5*x'*H*x))
        ik = mod(iter,n)+1;
        grad_quad = H*x;
        
        if(type == 0)
            alpha = x(ik)- (H(ik,:)*x - H(ik,ik)*x(ik)/(2*H(ik,ik)));
        elseif(type == 2)
           ik = randi(n);
        elseif(type == 3)          
           ik = find(abs(grad_quad) == max(abs(grad_quad)));
        end
        
        vect = zeros(n,1);
        vect(ik) =1;
        
        
        if(type==0)
            x = x- vect*alpha;
        else
            x = x- vect*alpha*grad_quad(ik);
        
        end
        if(norm(H*x,2)<= 10e-06*max(1,norm(H*x0,2)))
            break;
        end
        
end   
    
end

 
 

