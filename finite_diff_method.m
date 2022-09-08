%% Central Difference Method
% Daniele Kauctz Monteiro (2022)
% danielekauctz@hotmail.com

% M, K, C: mass, stiffness, and damping matrices of the dynamic model
% x, v, a: displacement, velocity, and acceleration matrices
           ... rows = degrees of freedom
           ... columns = number of time points
% x0, v0: initial displacement and velocity vectors
% F: force matrix
     ... rows = degrees of freedom
     ... columns = number of time points
% DOF: degrees of freedom
% dt: step (time interval)
% n: number of time points


function[x,v,a] = finite_diff_method(M,K,C,x0,v0,F,DOF,dt,n)
             
    x = zeros(DOF,n);
    x(:,1)= x0;       

    v = zeros(DOF,n-1);
    v(:,1)= v0;       

    a = zeros(DOF,n-1);
    a(:,1) = inv(M)*(F(:,1) -(C*v(:,1))-(K*x(:,1)));  
    
    x_1 = (dt^2/2*a(:,1)) - (dt*v(:,1)) + x(:,1);
    
    C1 = inv((M/dt^2)+(C/(2*dt)));
    C2 = K - (2/dt^2*M);
    C3 = M/dt^2 - (C/(2*dt));

    for i = 1:n-1
        
        if i == 1
            x(:,i+1) = C1*(F(:,i)-(C2*x(:,i))-(C3*x_1));
        else
            x(:,i+1) = C1*(F(:,i)-(C2*x(:,i))-(C3*x(:,i-1))); 
        end
        
    end
    
    for i = 1:n-2
                       
        v(:,i+1) = (x(:,i+2) - x(:,i))/(2*dt);
        a(:,i+1) = (1/dt^2)*(x(:,i+2) - (2*x(:,i+1)) + x(:,i));
        
    end

end