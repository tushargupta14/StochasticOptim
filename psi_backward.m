function psi_f = psi_backward(psi_i,fi,delta_t,G,w,y,z,theta,T,t_curr,t_next)
   
    %rho = 2.66*10^-12;
    %kv = 0.54;
    %b =1.45;
    %kb = 285;
    %E_b = 7517;
    
    tspan = [t_curr t_next];
    options = odeset('RelTol',1e-9,'AbsTol',1e-12);
    [t,psi] = ode45(@(t,psi)psiODE(t,psi,fi,G,w,y,z,theta,T),tspan,psi_i,options);
    psi_f = psi(end,:);
    
end
