close all
clear all

type =3;

for  n = [10,100,1000]
    
    a = 1;
    for rc = [1/10,1/100,1/1000,1/10000]
        
        R = sprandsym(n,0.1,rc,1);
        [x,iter(a)] = coordinate_minimisation(R,type);
        func_mag(a) = 0.5*x'*R*x;
        a = a+1;
    end
  
    figure;subplot(1,2,1)
    
    plot([10,100,1000,10000],func_mag,'r');
    title(strcat('Performance of algorithm  ', num2str(type+1), ' for dimsension: ',num2str(n)))
    xlabel('condition number')
    ylabel('Magnitude of function')
    
    subplot(1,2,2)
    plot([10,100,1000,10000],iter+1,'g')
    title(strcat('Performance of algorithm  ',num2str(type+1), ' for dimsension: ',num2str(n)))
    xlabel('condition number')
    ylabel('number of iterations')

end
    

