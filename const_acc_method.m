%% Constant Acceleration Method
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


function[x,v,a] = const_acc_method(M,K,C,x0,v0,F,DOF,dt,n)
             
    x = zeros(DOF,n);
    x(:,1)= x0;       

    v = zeros(DOF,n);
    v(:,1)= v0;       

    a = zeros(DOF,n);
    a(:,1) = inv(M)*(F(:,1) -(C*v(:,1))-(K*x(:,1)));  

    for i = 1:n-1

        v(:,i+1) = v(:,i) + (dt*a(:,i));
        x(:,i+1) = x(:,i) + (dt*v(:,i)) + (dt^2/2*a(:,i));
        
        a(:,i+1) = inv(M)*(F(:,i+1) -(C*v(:,i+1))-(K*x(:,i+1)));

    end

end