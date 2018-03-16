function zf = z_backward(zi,y,G,T,S,t_curr,t_next)
   
    %kv = 1.5;
  
    
    growth_rate = G;
    tspan = [t_curr t_next];
    options = odeset('RelTol',1e-9,'AbsTol',1e-12);
    [t,z] = ode45(@(t,z)zODE(t,z,y,growth_rate,T,S),[t_curr t_next],zi,options);
    
    zf = z(end,:);
    %zf = [z1 z2 z3 z4 z5 z6 z7 z8 z9];
end