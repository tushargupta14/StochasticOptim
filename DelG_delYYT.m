function G_YYT = DelG_delYYT(T,S,Cs,exper,dCS_dT)

E_g = 4859;
kg = 1.44 *10^8;
g = 1.5;
    

G_YYT = kg*g*(g-1)*exp(-E_g/T)*(E_g/(T^2)*S^(g-2)/Cs^2 + (g-2)*S^(g-3)*exper/Cs^2 - 2*dCS_dT/Cs^3);



end
