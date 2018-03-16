function dw_dt = wODE(t,w,z,y,G,T,S)


b =1.45;
kb = 285;
E_b = 7517;
kv = 0.54;
g = 1.5;
kg = 1.44 *10^8;
E_g = 4859;
rho = 2.66*10^-12;
DelG_Dely = DelG_dy(S,T);
DelB_Dely = DelB_dy(S,y(4),y(8),T);


d2G_dy2 = Del2G_DY2(S,T);
d2B_dy2 = Del2B_DY2(S,y(4),y(8),T);

A = -(2*DelG_Dely*(w(1)*-3*rho*kv*(y(4)+y(8))+y(2)*w(3)+2*y(3)*w(4)+3*y(4)*w(5)+y(6)*w(7)+2*y(7)*w(8)+3*y(8)*w(9))+2*w(6)*DelB_Dely);
B  = -(d2G_dy2*(-3*rho*kv*z(1)*(y(4)+y(8))+z(3)*y(2)+2*z(4)*y(3)+3*z(5)*y(4)+z(7)*y(6)+2*z(8)*y(7)+3*z(9)*y(8))+d2B_dy2*z(6));

dw1 = A+B;

dw2 = -2*w(3)*G;
dw3 = -2*w(4)*2*G;
dw4 = -(-3*rho*kv*G*2*w(1)+3*G*2*w(5));
dw5 = -2*w(6)*kb*exp(-E_b/T)*S^b;
dw6 = -2*G*w(7);
dw7 = -2*G*2*w(8);
dw8 = -(-3*rho*kv*G*2*w(1)+3*G*2*w(9));
dw9 = - 2*w(6)*kb*exp(-E_b/T)*S^b;

dw_dt = [dw1;dw2;dw3;dw4;dw5;dw6;dw7;dw8;dw9];




end
