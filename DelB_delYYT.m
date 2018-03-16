function B_YYT = DelB_delYYT(T,S,Cs,exper,dCS_dT,theta4,theta8,y4,y8)
b =1.45;
kb = 285;
E_b = 7517;

B_YYT = kb*b*(b-1)*exp(-E_b/T)*((E_b/T^2)*(S^(b-2)/Cs^2)*(y4+y8) + (b-2)*S^(b-3)*exper*(y4+y8)/Cs^2 -2* S^(b-2)*dCS_dT*(y4+y8)/Cs^3 +S^(b-2)*(theta4+theta8)/Cs^2); 


end
