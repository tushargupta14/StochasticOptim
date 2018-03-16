function yf = y_forward(yi,diffusion,delta_t,G,B,t_curr,t_next)

% g is the diffusion function 

%yf = zeros(1,9);
h = delta_t;
%or h = delta_t/10;
% y_ode is the drift function 
t_vector = t_curr:h:t_next;

opts = sdeset('SDEType','Ito');
y = sde_euler(@(t,y_ode)yODE(t,y_ode,G,B),diffusion,t_vector,yi,opts);
yf = y(end,:);

end