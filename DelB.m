function DB_dy = DelB(S,y4,y8,T)

b =1.45;
kb = 285;
E_b = 7517;
Cs = 6.29 * 10^-2 + 2.46*10^-3 * (T-273) - 7.14 * 10^-6 * (T-273)^2 ;
%Cm = 7.76 * 10^-2 + 2.46*10^-3 * (T-273) - 8.1*10^-6 * (T-273)^2 ;
DB_dy = kb*b*exp(-E_b/T)*S^(b-1)*(y4+y8)/Cs;

end 