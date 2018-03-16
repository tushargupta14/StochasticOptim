function dpsi_dt = psiODE(t,psi,fi,G,w,y,z,theta,T)
    rho = 2.66*10^-12;
    kv = 0.54;
    b =1.45;
    kb = 285;
    E_b = 7517;
    
    Cs = 6.29 * 10^-2 + 2.46*10^-3 * (T-273) - 7.14 * 10^-6 * (T-273)^2 ;
    %Cm = 7.76 * 10^-2 + 2.46*10^-3 * (T-273) - 8.1*10^-6 * (T-273)^2 ;
  
    S = (y(1) -Cs)/ Cs; % Relative Saturation 
    G_YT = DelG_y_T(T,S,y,theta);
    B_YT = DelB_y_T(T,S,theta,y);
    G_Y = DelG_dy(S,T);
    B_Y = DelB_dy(S,y(4),y(8),T);
    
    G_YY =  Del2G_DY2(S,T);
    B_YY = Del2B_DY2(S,y(4),y(8),T);
    G_T = DelG_T(S,theta,T,y(1));

    DCs_dT = 2.46*10^-3 -14.28*10^-6*(T-273);
    exper = (Cs*(theta(1) - DCs_dT) - DCs_dT * (y(1) - Cs))/Cs^2;
    G_YYT = DelG_delYYT(T,S,Cs,exper,DCs_dT);
    B_YYT = DelB_delYYT(T,S,Cs,exper,DCs_dT,theta(4),theta(8),y(4),y(8));
    
    A = -(6*rho*kv*(psi(1)*G_Y*(y(4)+y(8))+w(1)*G_YT*(y(4)+y(8))+w(1)*G_Y*(theta(4)+theta(8))));
    B = 2*(psi(3)*y(2)*G_Y+w(3)*theta(2)*G_Y+w(3)*y(2)*G_YT)+4*(psi(4)*G_Y*y(3)+w(4)*theta(3)*G_Y+w(4)*y(3)*G_YT)+6*(psi(5)*G_Y*y(4)+w(5)*theta(4)*G_Y+w(5)*y(4)*G_YT);
    C = 2*(B_YT*w(6)+B_Y*psi(6));
    D = 2*(psi(7)*G_Y*y(6)+w(7)*theta(6)*G_Y+w(7)*y(6)*G_YT)+4*(psi(8)*G_Y*y(7)+w(8)*G_Y*theta(7)+w(8)*y(7)*G_YT)+6*(psi(9)*G_Y*y(8)+w(9)*G_Y*theta(8)+w(9)*y(8)*G_YT);
    E = -3*rho*kv*(fi(1)*G_YY*(y(4)+y(8))+z(1)*G_YY*(theta(4)+theta(8))+z(1)*G_YYT*(y(4)+y(8)));
    F = fi(3)*G_YY*y(2)+ z(3)*G_YY*theta(2)+z(3)*G_YYT*y(2) + 2*(fi(4)*y(3)*G_YY+z(4)*G_YYT*y(3)+z(4)*theta(3)*G_YY)+ 3*(fi(5)*G_YY*y(4)+z(5)*G_YYT*y(4)+z(5)*G_YY*theta(4));
    % G reserved for the growth rate
    H = B_YYT*z(6)+B_YY*fi(6) + fi(7)*G_YY*y(6)+z(7)*G_YYT*y(6)+z(7)*G_YY*theta(6)+ 2*(fi(8)*G_YY*y(7)+z(8)*G_YYT*y(7)+z(8)*G_YY*theta(7)) + 3*(fi(9)*G_YY*y(8)+z(9)*G_YYT*y(8)+z(9)*G_YY*theta(8));
  
    dpsi1 = -(A+B+C+D+E+F+H);
    dpsi2 = -2*(psi(3)*G +w(3)*G_T);
    dpsi3 = -4*(psi(4)*G +w(4)*G_T);
    dpsi4 = -(-6*rho*kv*(psi(1)*G+w(1)*G_T)+6*(psi(5)*G+w(5)*G_T));
    
    dpsi5 = -2*kb*exp(-E_b/T)*(psi(6)*S^b + w(6)*(E_b/(T^2))*S^b + w(6)*b*S^(b-1)*exper);
    dpsi6 = -2*(psi(7)*G+w(7)*G_T);
    dpsi7 = -4*(psi(8)*G+w(8)*G_T);
    dpsi8 = -(-6*rho*kv*(psi(1)*G+w(1)*G_T)+ 6*(psi(9)*G+w(9)*G_T));
    dpsi9 = -2*kb*exp(-E_b/T)*(psi(6)*S^b + w(6)*(E_b/(T^2))*S^b + w(6)*b*S^(b-1)*exper);
    
    dpsi_dt = [dpsi1;dpsi2;dpsi3;dpsi4;dpsi5;dpsi6;dpsi7;dpsi8;dpsi9];

end