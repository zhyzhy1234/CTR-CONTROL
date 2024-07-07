%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project name: CTR-control
%Programer   : zhy
%Finish date : 
%Records     : 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
close all;

load('data_n_w.mat', 'n_w'); load('data_epsilon.mat', 'epsilon'); load('data_r_w.mat', 'r_w');load('data_L.mat', 'L'); 
load('data_eta_s0.mat', 'eta_s0');load('data_eta_l0.mat', 'eta_l0'); load('data_l_w0.mat', 'l_w0'); 
load('data_eta_sl.mat', 'eta_sl');load('data_eta_ll.mat', 'eta_ll'); load('data_l_wl.mat', 'l_wl');
load('data_eta_sr.mat', 'eta_sr');load('data_eta_lr.mat', 'eta_lr'); load('data_l_wr.mat', 'l_wr')

global Step_simu  Time_simu 
global Num_d newSt_d Exp_d 
global Exp_a newSt_a v_a
global ref_point w dot_w   
global error s_d delt_s 

Step_simu = 5000;                        
Time_simu = zeros(Step_simu,1);          

Num_d = 3;                                
eta_point = zeros(Step_simu+1,6,3);      
 
newSt_d = zeros(Step_simu+1,6,Num_d);     

newSt_d(1,:,1) = [0 0 0 0 0.5 1.2]; 
newSt_d(1,:,2) = [0 0 0 0 0.8 2.1];  
newSt_d(1,:,3) = [0 0 0 0 1.0 0.9];  


alpha = [0.7 0.6 0.6];

Exp_d = zeros(Step_simu+1,2,Num_d);      

newSt_a        = zeros(Step_simu+1,6);   
newSt_a(1,5:6) = [1.27,1.5];              
newSt_a(1,4)   = 0;                       
Exp_a          = zeros(Step_simu,3);     
v_a0 = 0.6/3;                             

u_a      = zeros(Step_simu+1,1);          
u_d      = zeros(Step_simu+1,3);          
kappa_a     = zeros(Step_simu+1,1);        
point_af = zeros(Step_simu+1,2);         
point_df = zeros(Step_simu+1,2,3);        
index_w  = zeros(Step_simu+1,3);          
l_y      = zeros(Step_simu+1,3);          
fai_p    = zeros(Step_simu+1,3);          

dot_w = zeros(Step_simu+1,3);             

w = ones(Step_simu+1,3);

error  = zeros(Step_simu,2,Num_d);       
s_d    = zeros(Step_simu,Num_d);         
delt_s = zeros(Step_simu,Num_d);
      
Sum_sav = zeros(Step_simu,Num_d);

k =1;     
gama = 1; 
kd = [0 2 2];   

c_c1 = 0;

newSt_d(1,:,1) = [0 0 0 0 0.5 1.2];  
newSt_d(1,:,2) = [0 0 0 0 0.8 2.1];  
newSt_d(1,:,3) = [0 0 0 0 1.0 0.9]; 

v_a = zeros(Step_simu,2); 
theta = zeros(Step_simu,1);      
old_theta = pi/2; old_rhod = 0;    
dot_theta =  zeros(Step_simu,1); 
dot_rhod =  zeros(Step_simu,1);    
beta  =  zeros(Step_simu,1); 
rho_as =  zeros(Step_simu,1); 
lh_d   =  zeros(Step_simu,1);
mid_var0 =  zeros(Step_simu,1); 
gamma_beta  = zeros(Step_simu,1); 
gamma_rho  = zeros(Step_simu,1);    

gamma_p  = zeros(Step_simu,1); 
gamma_s  = zeros(Step_simu,1); 


l_set = zeros(Step_simu,1);
rho_d  = zeros(Step_simu,1);
beta_d =  zeros(Step_simu,1);

e_beta = zeros(Step_simu,1);
e_rho  = zeros(Step_simu,1);

z_beta = zeros(Step_simu,1);
z_rho   = zeros(Step_simu,1);

sigma_beta = zeros(Step_simu,1);
sigma_rho    = zeros(Step_simu,1); 

chi_d = zeros(Step_simu,2,Num_d); 


alphaf = 0.1; initial_value = 0;
filtered_speed = zeros(Step_simu,3); 
        
  
color_h = [78/255 098/255 171/255];
color_er = [70/255 158/255 180/255];
color_el = [135/255 207/255 164/255];

figure;   

mainAx = axes; axis([0 10 0 10]);
plot(mainAx,eta_l0(:,1),eta_l0(:,2),'-.');hold on;
plot(mainAx,eta_ll(:,1),eta_ll(:,2));hold on;
plot(mainAx,eta_lr(:,1),eta_lr(:,2));hold on;

delta_s = 8*ones(Step_simu,1);
for m = 1:Step_simu
    
    count = m;
    Time_simu(m) = m/10;

    waitforbuttonpress; 
    key = get(gcf, 'CurrentCharacter');

    switch key
        case 'w' 
            Exp_a(m,1) = v_a0; Exp_a(m,3) = 0;    
        case 's' 
            Exp_a(m,1) = -v_a0; Exp_a(m,3) = 0;   
        case 'a' 
            Exp_a(m,1) = v_a0; Exp_a(m,3) = 0.8;  
        case 'd' 
            Exp_a(m,1) = v_a0; Exp_a(m,3) = -0.8;
        case 'stop' 
            Exp_a(m,1) = 0; Exp_a(m,3) = 0;
     end
    

    for i = 1:Num_d 
        [rho_da(m,i), theta_da(m,i)] = Fun_CalDistanAndAngle(newSt_a(m,5:6),newSt_d(m,5:6,i));  
    end
    
    newSt_a(m+1,4) = 0.1*Exp_a(m,3) + newSt_a(m,4);
    
    if rho_da(m,1)< 0.25 || rho_da(m,2)< 0.25 || rho_da(m,3)< 0.25
        [min_difference, index] = sort(rho_da(m,:));
        [aaa,newSt_a(m+1,4)] = Fun_CalDistanAndAngle([cos(theta_da(m,index(1))), sin(theta_da(m,index(1)))] + [cos(theta_da(m,index(2))), sin(theta_da(m,index(2)))],[0,0]);
    end
    newSt_a(m+1,5) = 0.1*Exp_a(m,1)*cos(newSt_a(m+1,4)) + newSt_a(m,5);
    newSt_a(m+1,6) = 0.1*Exp_a(m,1)*sin(newSt_a(m+1,4)) + newSt_a(m,6);                           
    
    v_a(m,1) = Exp_a(1)*cos(newSt_a(m,4));                                                     
    v_a(m,2) = Exp_a(1)*sin(newSt_a(m,4));                                                      
    
     
     center = [9.5, 8.5]; rho_at(m,1) = norm(newSt_a(m,5:6)-center,2);
     
     if rho_at(m,1) <= epsilon
        break; 
     end
     
    [rho_dd(m,1), theta_dd(m,1)] = Fun_CalDistanAndAngle(newSt_d(m,5:6,2),newSt_d(m,5:6,3));   
    [rho_dd(m,2), theta_dd(m,2)] = Fun_CalDistanAndAngle(newSt_d(m,5:6,1),newSt_d(m,5:6,3));   

    gamma_p_0 = 0.1; gamma_p_inf = 0.05; k_p = 0.05; 
    gamma_s_0 = 0.2; gamma_s_inf = 0.05; k_s = 0.05;
    
    gamma_p(m) = (gamma_p_0 - gamma_p_inf).*exp(-k_p*m/10) + gamma_p_inf;
    gamma_s(m) = (gamma_s_0 - gamma_s_inf).*exp(-k_s*m/10) + gamma_s_inf;
    
   %% 
   
    point_df(m,:,2) = eta_ll(w(m,2),1:2);  fai_p(m,2) = eta_ll(w(m,2),3);                              
    point_df(m,:,3) = eta_lr(w(m,3),1:2);  fai_p(m,3) = eta_lr(w(m,3),3);                              
    
    
    [index_w(m+1,1),l_y(m+1,1),point_af(m+1,:)] = Fun_CartesianToSf(newSt_a(m+1,5:6), eta_l0(:,1:2));  
    [index_w(m,2),l_y(m,2),point_df(m,:,2)] = Fun_CartesianToSf(newSt_d(m,5:6,2), eta_ll(:,1:2));      
    [index_w(m,3),l_y(m,3),point_df(m,:,3)] = Fun_CartesianToSf(newSt_d(m,5:6,3), eta_lr(:,1:2));     
    
    fai_i(m,2) = eta_ll(index_w(m,2),3);  fai_i(m,3) = eta_lr(index_w(m,3),3);                        
    
    kappa_a(m+1,1) = eta_l0(index_w(m+1,1),4);  
    
    s_a(m+1,1) = index_w(m+1,1)/100;
    u_a(m+1,1)  = (index_w(m+1,1)/100 - index_w(m,1)/100) /0.1; 
     
    theta_wi = 0.826; c_a = 1.2;
    delta_s0 = min([2*epsilon - rho_da(m,3)*sin(theta_da(m,3)-fai_i(m,3)), rho_da(m,3)*sin(theta_da(m,3)-fai_i(m,3)), epsilon]./tan(c_a*theta_wi));
    
    mid_s = eta_l0(index_w(m+1,1),10);
    
    if mid_s  > delta_s0
        mid_s = delta_s0;
    end
        
    delta_s(m+1,1) = delta_s0 + l_y(m+1,1)*kappa_a(m+1,1)*eta_l0(index_w(m+1,1),9)*mid_s;
    s_av(m+1,1)    = floor(index_w(m+1,1) - delta_s(m+1,1)*100);  
    
    if s_av(m+1,1)  < 1
        s_av(m+1,1)  = 1;
    end
    
    dot_sav(m+1,1)  = (s_av(m+1,1)/100 - s_av(m,1)/100 )/0.1;
    kappa_av(m+1,1) = eta_l0(s_av(m+1,1),4);                      
    piont_c         = eta_l0(s_av(m+1,1),5:6);                    
    piont_c2        = eta_l0(s_av(m+1,1),1:2);                    
    
    Sum_sav(m+1,1)   = -0.1*eta_l0(s_av(m+1,1),9) * kappa_av(m+1,1)* dot_sav(m+1,1) + Sum_sav(m,1);                        
    
    s_d(m,1) = s_av(m+1,1)/100;                                   

    
    for i = 2:Num_d
        
        s_d(m,i)    = s_d(m,1) + L(i,2)*Sum_sav(m+1,1);           
        delt_s(m,i) = s_d(m,i) - index_w(m,i)/100;                
        
        z_s(m,i) = delt_s(m,i)/gamma_s(m);
        r_s(m,i) = 1 - z_s(m,i)^2;
        sigma_s(m,i)  = atanh(z_s(m,i));
        
        u_d(m+1,i) = (1 - L(i,2)*eta_l0(s_av(m+1,1),9) * kappa_av(m+1,1))* dot_sav(m+1,1) + 2*sigma_s(m,i);
        
        if u_d(m+1,i) < 0 
             u_d(m+1,i) = 0;
        end
                
        if u_d(m+1,i) > 1/3         
             u_d(m+1,i) = 1/3; 
        end
        

        error(m,:,i) = [cos(fai_p(m,i)) sin(fai_p(m,i));sin(fai_p(m,i)) -cos(fai_p(m,i))]...
                          *(newSt_d(m,5:6,i) - point_df(m,1:2,i))';                                   
                      
        error_p(m,i) = norm(error(m,:,i),2);
         
        error_t(m,i) = norm(error(m,1,i),2);
        error_n(m,i) = norm(error(m,2,i),2);
        
        z_n(m,i) = error(m,2,i)/gamma_p(m);
        r_n(m,i) = 1 - z_n(m,i)^2;
        sigma_n(m,i)  = atanh(z_n(m,i));
        
        z_t(m,i) = error(m,1,i)/gamma_p(m);
        r_t(m,i) = 1 - z_t(m,i)^2;
        sigma_t(m,i)  = atanh(z_t(m,i));
        
        Exp_d(m,2,i) = fai_p(m,i) + atan(k*r_n(m,i)*sigma_n(m,i));                                              

        if m>1
            dot_fai_p(m,i) =( fai_p(m,i) - fai_p(m-1,i))/0.1;
        else
            dot_fai_p(m,i) =0;
        end
        if abs(sigma_t(m,i))>0.01 && abs(r_n(m,i))>0.01
            mmm(m,i) = r_t(m,i)*sigma_n(m,i)*(error(m,1,i)*dot_fai_p(m,i)- z_n(m,i)*k_p*(gamma_p_0 - gamma_p_inf).*exp(-k_p*m/10)/10)/(r_n(m,i)*sigma_t(m,i));
        else
           mmm(m,i) =0;
        end
        
        u_p(m,i) = u_d(m+1,i)*cos(Exp_d(m,2,i) - fai_p(m,i)) + 2*r_t(m,i)*sigma_t(m,i) - ...
                   error(m,2,i)*dot_fai_p(m,i) + z_t(m,i)*k_p*(gamma_p_0 - gamma_p_inf).*exp(-k_p*m/10)/10 + mmm(m,i);
        
        
        if u_p(m,i) < 0
            u_p(m,i) = 0;
        end

        w(m+1,i) = floor(0.1*u_p(m,i)*100 + w(m,i));  
        if w(m+1,i)>size(eta_l0,1)
            fin_count = m;
            break;
        end
        
        
        filtered_speed(m+1,i) = firstOrderFilter(alphaf, u_d(m+1,i), filtered_speed(m,i));

        Exp_d(m,1,i) = tanh(filtered_speed(m+1,i));

        newSt_d(m+1,5:6,i) = 0.1*Exp_d(m,1,i)*[cos(Exp_d(m,2,i)), sin(Exp_d(m,2,i))] + newSt_d(m,5:6,i);   
        newSt_d(m+1,1,i)   = Exp_d(m,1,i);
        newSt_d(m+1,4,i)   = Exp_d(m,2,i);
 
    end
    
   %% 
    
    theta(m) = theta_dd(m,1) - fai_i(m,3) ; 
    dot_theta(m) = ( theta(m) - old_theta)/0.1;
    old_theta = theta(m);
    
    [rho(m), beta(m)] = Fun_CalDistanAndAngle(newSt_d(m,5:6,1),newSt_a(m,5:6));     
    
    c_rho = 5; c_beta = 2;
    l_0 = 2; ld = 0.5;
    
    gamma_beta_0 = 0.7; gamma_beta_inf = 0.1; k_beta = 0.05; 
    gamma_rho_0 = 5;  gamma_rho_inf = 0.1; k_rho = 0.05;
    
    gamma_beta(m) = (gamma_beta_0 - gamma_beta_inf).*exp(-k_beta*m/10) + gamma_beta_inf;
    gamma_rho(m)  = (gamma_rho_0 - gamma_rho_inf).*exp(-k_rho*m/10) + gamma_rho_inf;
    
    ee = rho(m);
    l_01 = (1-alpha(3)^2)*gamma_rho_0/(cos(gamma_beta_0) - alpha(3)) - 0.3;
    if l_01 > (rho(m) - ee)
        l_set(m) = rho(m) - ee;
    else
        l_set(m) = rho(m) - l_01;
    end
     
    mid_var0(m) = abs(sin(beta(m) - fai_i(m,3) - theta(m)));
    rho_as(m)   = rho_da(m,3)*sin(theta(m) - theta_da(m,3) + fai_i(m,3))/mid_var0(m);
    lh_d(m)  = l_set(m)*(mid_var0(m) - alpha(1))/(1 - alpha(1)^2);
    rho_d(m) = lh_d(m)/mid_var0(m) + rho_as(m) ;
    beta_d(m) = pi/2 + theta(m) + fai_i(m,3);
    
    dot_rhod(m) = ( rho_d(m) - old_rhod )./0.1;
    old_rhod = rho_d(m);
    
    e_beta(m) = beta_d(m) - beta(m);
    e_rho(m)  = rho_d(m) -  rho(m);
    
    z_beta(m) = e_beta(m)/gamma_beta(m);
    z_rho(m)  = e_rho(m)/gamma_rho(m);
    
    sigma_beta(m) = atanh(z_beta(m));
    sigma_rho(m)  = atanh(z_rho(m));
    
    chi_d(m,:,i) = [cos(beta(m)) -sin(beta(m)); sin(beta(m)) cos(beta(m))]*...
                   [ c_rho*sigma_rho(m) + dot_rhod(m) + cos(beta(m))*v_a(m,1) + sin(beta(m))*v_a(m,2);...
                     rho(m)*(c_beta*sigma_beta(m) + dot_theta(m)) - sin(beta(m))*v_a(m,1) + cos(beta(m))*v_a(m,2)];
                 
    u_d1 = 0.6*tanh(norm(chi_d(m,:,i),2))/3/0.2;
    
    fai1 = atan2(chi_d(m,2,i),chi_d(m,1,i)); 

    Exp_d(m,:,1) =  [u_d1,fai1];    
    newSt_d(m+1,5,1) = 0.1*u_d1*cos(fai1) + newSt_d(m,5,1); 
    newSt_d(m+1,6,1) = 0.1*u_d1*sin(fai1) + newSt_d(m,6,1);
               
    rho_as(m)   = rho_da(m,3)*sin(theta(m) - theta_da(m,3) + fai_i(m,3))/mid_var0(m);
    
   
    cla(mainAx);

    axis(mainAx,[0 10 0 10]);

  
%% 
    vector_length = 0.5;                                       
    x_component = cos(newSt_a(m,4))*vector_length; y_component = sin(newSt_a(m,4))*vector_length;
    quiver(mainAx,newSt_a(m,5), newSt_a(m,6), x_component, y_component, 'r', 'LineWidth', 1);  hold on;    

    plot(mainAx,newSt_a(1:m,5), newSt_a(1:m,6), 'k-'); hold on;  
    
%%
    scatter(mainAx,newSt_a(m,5),   newSt_a(m,6),   'SizeData', 20,'MarkerEdgeColor', 'none', 'MarkerFaceColor',[220/255 0/255 60/255]);hold on;% p_a 
    scatter(mainAx,newSt_d(m,5,1), newSt_d(m,6,1), 'SizeData', 20,'MarkerEdgeColor', 'none', 'MarkerFaceColor',[0/255 128/255 0/255]); hold on;% p_h
    scatter(mainAx,newSt_d(m,5,2), newSt_d(m,6,2), 'SizeData', 20,'MarkerEdgeColor', 'none', 'MarkerFaceColor',[0/255 128/255 0/255]); hold on;% p_el
    scatter(mainAx,newSt_d(m,5,3), newSt_d(m,6,3), 'SizeData', 20,'MarkerEdgeColor', 'none', 'MarkerFaceColor',[0/255 128/255 0/255]); hold on;% p_er
    
    
    plot(mainAx,newSt_d(1:m+1,5,1), newSt_d(1:m+1,6,1), 'g-'); hold on;   
    plot(mainAx,newSt_d(1:m+1,5,2), newSt_d(1:m+1,6,2), 'g-'); hold on;      
    plot(mainAx,newSt_d(1:m+1,5,3), newSt_d(1:m+1,6,3), 'k-'); hold on;  
    

    h2 = plot(mainAx,newSt_d((1:m+1),5,1),newSt_d((1:m+1),6,1),'color',[color_el,0.2],'linewidth',10,'lineStyle','-');hold on;
    h3 = plot(mainAx,newSt_d((1:m+1),5,2),newSt_d((1:m+1),6,2),'color',[color_h,0.2],'linewidth',10,'lineStyle','-');hold on;
    h4 = plot(mainAx,newSt_d((1:m+1),5,3),newSt_d((1:m+1),6,3),'color',[color_er,0.2],'linewidth',10,'lineStyle','-');hold on;


    h1 = plot(mainAx,newSt_a((1:m+1),5),newSt_a((1:m+1),6),'color','r','linewidth',1,'lineStyle','-.');hold on;
    h2 = plot(mainAx,newSt_d((1:m+1),5,1),newSt_d((1:m+1),6,1),'color',[color_el,1],'linewidth',1,'lineStyle','-');hold on;
    h3 = plot(mainAx,newSt_d((1:m+1),5,2),newSt_d((1:m+1),6,2),'color',[color_h,1],'linewidth',1,'lineStyle','-');hold on;
    h4 = plot(mainAx,newSt_d((1:m+1),5,3),newSt_d((1:m+1),6,3),'color',[color_er,1],'linewidth',1,'lineStyle','-');hold on;
    
    %viscircles(mainAx,newSt_d(m+1,5:6,1), 0.2, 'EdgeColor', 'k', 'LineStyle', '-.','LineWidth', 1); hold on;
    viscircles(mainAx,newSt_d(m+1,5:6,2), 0.2, 'EdgeColor', 'k', 'LineStyle', '-.','LineWidth', 1); hold on;
    viscircles(mainAx,newSt_d(m+1,5:6,3), 0.2, 'EdgeColor', 'k', 'LineStyle', '-.','LineWidth', 1); hold on;
    axis([0 10 0 10]);  
    
    pause(0.1);
    
    
end