%%% Main file for simulating stochastic part of the paper 

b =1.45;
kb = 285;
E_b = 7517;
kv = 0.54;
g = 1.5;
kg = 1.44 *10^8;
E_g = 4859;
rho = 2.66*10^-12;

t0 = 0;
delta_t = 1;
batch_time = 1800;
tf = batch_time;
dt = 1; %% same as delta_t

n = length(t0:dt:tf);

y0 = [0.1743 66.66 1.83*10^4 5.05*10^6 1.93*10^9 0.867 0 0 0];
zf = [0 0 0 0 1 0 0 0 -1];
wf= [0 0 0 0 1 0 0 0 -1];
theta0 = [0 0 0 0 0 0 0 0 0];
fi_f = [0 0 0 0 0 0 0 0 0];
psi_f = [0 0 0 0 0 0 0 0 0];

N = 3; % number of iterations 
tolerance = 10^4;
M = -10^-7;
T_vec = ones(1,length(t0:delta_t:batch_time))*323;
DH_vec = zeros(N,length(t0:delta_t:batch_time));

iteration = 1;
diffusion = [2.659*10^-5;0;25.882;1.517*10^4;6.57*10^6;0.5486;25.91599;1382.3464;8.7453*10^4];
g1 = 2.659*10^-5;
g2 = 0;
g3 = 25.882;
g4 = 1.517*10^4;
g5 = 6.57*10^6;
g6 = 0.5486;
g7 = 25.91599;
g8 = 1382.3464;
g9 = 8.7453*10^4;

iteration
cmap = hsv(N);

figure;hold on;
title('Hamiltonian Derivative')

while iteration < N
   
    y_mat = zeros(length(t0:delta_t:tf),9);
    z_mat = zeros(length(t0:delta_t:tf),9);
    w_mat = zeros(length(t0:delta_t:tf),9);
    theta_mat = zeros(length(t0:delta_t:tf),9);
    fi_mat = zeros(length(t0:delta_t:tf),9);
    psi_mat = zeros(length(t0:delta_t:tf),9);
    y_mat(1,:) = y0;
    z_mat(end,:) = zf;
    w_mat(end,:) = wf;
    theta_mat(1,:) = theta0;
    fi_mat(end,:) = fi_f;
    psi_mat(end,:) = psi_f;
    DelH_dy_mat = zeros(length(t0:delta_t:tf),9); %1801*9 
    DelH_dz_mat = zeros(length(t0:delta_t:tf),9);

    
    
    % y forward integration 
    for t = 1:1:(n-1)
        t_curr = t;
        t_next = t+1;
        Cs = 6.29 * 10^-2 + 2.46*10^-3 * (T_vec(1,t_curr)-273) - 7.14 * 10^-6 * (T_vec(1,t_curr)-273)^2 ;
        S  = (y_mat(t_curr,1) - Cs)/Cs;
        G = (kg* exp(-E_g/T_vec(1,t_curr)))*S^g;
        %T_vec(1,t_curr)
        B = (kb*exp(-E_b/T_vec(1,t_curr)))*S^b*(y_mat(t_curr,5)+y_mat(t_curr,9));    
        y_mat(t_next,:) = y_forward(y_mat(t_curr,:),diffusion,delta_t,G,B,t_curr,t_next);   
    end
    %y_mat
    display('Calculated y');
    for t = n:-1:2
        t_curr = t;
        t_next = t-1;
        Cs = 6.29 * 10^-2 + 2.46*10^-3 * (T_vec(1,t_curr)-273) - 7.14 * 10^-6 * (T_vec(1,t_curr)-273)^2 ;%%%%
        S  = (y_mat(t_curr,1) - Cs)/Cs; %%%%%%
        G = (kg* exp(-E_g/T_vec(1,t_curr)))*S^g; %%%%%
        z_mat(t_next,:) = z_backward(z_mat(t_curr,:),y_mat(t_curr,:),G,T_vec(1,t_curr),S,t_curr,t_next);  %%%%%
    end
    display('Calculated z');
    % backward integration for w
    for t = n:-1:2
        t_curr = t;
        t_next = t-1;
        Cs = 6.29 * 10^-2 + 2.46*10^-3 * (T_vec(1,t_curr)-273) - 7.14 * 10^-6 * (T_vec(1,t_curr)-273)^2 ;%%%%
        S  = (y_mat(t_curr,1) - Cs)/Cs; %%%%%%
        G = (kg* exp(-E_g/T_vec(1,t_curr)))*S^g; %%%%%
        w_mat(t_next,:) = w_backward(w_mat(t_curr,:),z_mat(t_curr,:),y_mat(t_curr,:),G,T_vec(1,t_curr),S,t_curr,t_next);  %%%%%
    end
    display('Calculated w');
    for t_curr = 1:1:(n-1) 
        
        Cs = 6.29 * 10^-2 + 2.46*10^-3 * (T_vec(1,t_curr)-273) - 7.14 * 10^-6 * (T_vec(1,t_curr)-273)^2 ;
        S  = (y_mat(t_curr,1) - Cs)/Cs;
        G = (kg* exp(-E_g/T_vec(1,t_curr)))*S^g;
        DelG_dT = DelG_T(S,theta_mat(t_curr,:),T_vec(1,t_curr),y_mat(t_curr,1));
        DelB_dT = DelB_T(S,T_vec(1,t_curr),theta_mat(t_curr,:),y_mat(t_curr,:));
        t_next = t_curr+1;
        theta_mat(t_next,:) = theta_forward(theta_mat(t_curr,:),delta_t,G,DelG_dT,DelB_dT,y_mat(t_curr,:),t_curr,t_next);
        
    end 
    display('Calculated theta');
    for t = n:-1:2
        t_curr = t;
        t_next = t-1;
        Cs = 6.29 * 10^-2 + 2.46*10^-3 * (T_vec(1,t_curr)-273) - 7.14 * 10^-6 * (T_vec(1,t_curr)-273)^2 ; %%%%
        S  = (y_mat(t_curr,1) - Cs)/Cs; %%%%
        G = (kg* exp(-E_g/T_vec(1,t_curr)))*S^g;%%%
        
        fi_mat(t_next,:) = fi_backward(fi_mat(t_curr,:),delta_t,G,y_mat(t_curr,:),z_mat(t_curr,:),theta_mat(t_curr,:),T_vec(1,t_curr),t_curr,t_next);
   
    end
    
   display('Calculated fi');
    % Backward integration for psi
    
    for t = n:-1:2
        t_curr = t;
        t_next = t-1;
        Cs = 6.29 * 10^-2 + 2.46*10^-3 * (T_vec(1,t_curr)-273) - 7.14 * 10^-6 * (T_vec(1,t_curr)-273)^2 ; %%%%
        S  = (y_mat(t_curr,1) - Cs)/Cs; %%%%
        G = (kg* exp(-E_g/T_vec(1,t_curr)))*S^g;
       
        psi_mat(t_next,:) = psi_backward(psi_mat(t_curr,:),fi_mat(t_curr,:),delta_t,G,w_mat(t_curr,:),y_mat(t_curr,:),z_mat(t_curr,:),theta_mat(t_curr,:),T_vec(1,t_curr),t_curr,t_next);
   
    end
    display('Calculated psi');
    %%%% Calculating the value of DH_vec
    for t_curr = 1:1:n
        
        Cs = 6.29 * 10^-2 + 2.46*10^-3 * (T_vec(1,t_curr)-273) - 7.14 * 10^-6 * (T_vec(1,t_curr)-273)^2 ;
        S  = (y_mat(t_curr,1) - Cs)/Cs;
        G = (kg* exp(-E_g/T_vec(1,t_curr)))*S^g;
        B = (kb*exp(-E_b/T_vec(1,t_curr)))*S^b*(y_mat(t_curr,5)+y_mat(t_curr,9));                   
        DelH_dy_mat(t_curr,:) = DelH_y(G,z_mat(t_curr,:),y_mat(t_curr,:),S,T_vec(1,t_curr));
        DelH_dz_mat(t_curr,:) = DelH_z(G,B,rho,kv,y_mat(t_curr,:))*10^-6;
   
        %display('done')     
    end
    
    DelH_dw = [g1^2/2 g2^2/2 g3^2/2 g4^2/2 g5^2/2 g6^2/2 g7^2/2 g8^2/2 g9^2/2];
    for t = 1:1:n
        sum = 0;
        for i=1:1:9
            sum = sum + DelH_dy_mat(t,i)*theta_mat(t,i) + DelH_dz_mat(n-t+1,i)*fi_mat(n-t+1,i)+DelH_dw(i)*psi_mat(n-t+1,i);    %%$$$$$$$$$$$
            %theta_mat(t,i)
        end
        DH_vec(iteration,t)=sum;
    end
    
    for t = 1:1:n
       
        if abs(DH_vec(iteration,t)) > tolerance
            T_vec(1,t)  = check_constraint1(y_mat(t,1),T_vec(1,t),M,DH_vec(iteration,t));
            %T_vec(t) = T_vec(t) + M*DH_vec(iteration,t);
        end 
    end
    
    plot(t0:delta_t:tf,DH_vec(iteration,:),'Color',cmap(iteration,:));
    iteration = iteration+1
    
    
end

hold off;

%plot_T_vec(T_vec,n);
figure
plot(t0:delta_t:tf,T_vec)
title('Teperature Profile');


