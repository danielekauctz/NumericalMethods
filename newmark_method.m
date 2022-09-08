%% Newmark's Method
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


function[x,v,a] = newmark_method(M,K,C,x0,v0,F,DOF,dt,n)
             
    x = zeros(DOF,n);
    x(:,1)= x0;       

    v = zeros(DOF,n);
    v(:,1)= v0;       

    a = zeros(DOF,n);
    a(:,1) = inv(M)*(F(:,1) -(C*v0)-(K*x0));  
   

    alpha = 0.25;     % alpha (= 1/4 or 1/6)
    delta  = 0.50;    % delta
    
    b0 = 1/(alpha*dt^2);
    b1 = 1/(alpha*dt);
    b2 = (1/(2*alpha)) - 1;
    b5 = delta/(alpha*dt);
    b6 = (delta/alpha) - 1;
    b7 = (dt/2)*((delta/alpha) - 2);
    c  = inv((b0*M) + (b5*C) + K);

    for i = 1:n-1

        x(:,i+1) = c*(F(:,i+1) + (M*((b0*x(:,i))+(b1*v(:,i))+(b2*a(:,i)))) + (C*((b5*x(:,i))+(b6*v(:,i))+(b7*a(:,i)))));
        v(:,i+1) = (b5*(x(:,i+1)-x(:,i))) - (b6*v(:,i)) - (b7*a(:,i));
        a(:,i+1) = (b0*(x(:,i+1)-x(:,i))) - (b1*v(:,i)) - (b2*a(:,i));

    end

end