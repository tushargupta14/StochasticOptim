function wf = w_backward(wi,z,y,G,T,S,t_curr,t_next)
   
    %kv = 1.5;
  
    
    growth_rate = G;
    tspan = [t_curr t_next];
    options = odeset('RelTol',1e-9,'AbsTol',1e-12);
    [t,w] = ode45(@(t,w)wODE(t,w,z,y,growth_rate,T,S),[t_curr t_next],wi,options);
    
    wf = w(end,:);
    %zf = [z1 z2 z3 z4 z5 z6 z7 z8 z9];
end