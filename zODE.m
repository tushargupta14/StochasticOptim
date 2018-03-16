function dz_dt = zODE(t,z,y,G,T,S)

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

    dz1 = -(-3*rho*kv*DelG_Dely*(y(4)+y(8))*z(1)+ DelG_Dely*y(2)*z(3)+2*DelG_Dely*y(3)*z(4)+3*DelG_Dely*y(4)*z(5)+DelB_Dely*z(6)+DelG_Dely*y(6)*z(7)+2*DelG_Dely*y(7)*z(8)+3*DelG_Dely*y(8)*z(9));
    
    dz2 = -G*z(3);
    dz3 = -2*G*z(4);
    dz4 =  -(-3*rho*kv*G*z(1) + 3*G*z(5));
    dz5 = -z(6)*kb*exp(-E_b/T)*S^b;
    dz6 = -G*z(7);
    dz7 = -2*G*z(8);
    dz8 = -(-3*rho*kv*G*z(1) + 3*G*z(9));
    dz9 = -z(6)*kb*exp(-E_b/T)*S^b;

    dz_dt = [dz1;dz2;dz3;dz4;dz5;dz6;dz7;dz8;dz9];

end
