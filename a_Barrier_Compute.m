%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project name: barrier
%Programer   : zhy
%Finish date : 
%Records     : 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
clc
clear all
close all

color1 = [255/256 47/256 47/256];
color2 = [0/256 205/256 0];       
%load('data_epsilon.mat', 'epsilon'); 
load('data_r_w.mat', 'r_w');
load('data_L.mat', 'L'); 

alpha = 0.6;           
r_out = r_w + L(3,2);  
r_inn = r_w + L(2,2);  
l_y = (r_out-r_inn)/2; 


Cr_a = linspace(r_inn,r_out, 1000); 
r_i = r_out; Cf_ar_out = -asin(alpha*r_i./Cr_a) + asin(alpha) + sqrt(1/(alpha^2)-1) - sqrt((Cr_a./(alpha*r_i)).^2 -1); 
r_i = r_inn; Cf_ar_inn =  asin(alpha*r_i./Cr_a) - asin(alpha) - sqrt(1/(alpha^2)-1) + sqrt((Cr_a./(alpha*r_i)).^2 -1); 

rho_ai   = linspace(0,r_out, 1000); 
r_i = r_out; Lr_a_out = sqrt(r_i^2 - 2*alpha*r_i.*rho_ai + rho_ai.^2); 
Lf_ar_out = asin(alpha) + asin(sqrt(1-alpha^2)*r_i./Lr_a_out) -pi/2;    

r_i = r_inn; Lr_a_inn = sqrt(r_i^2 + 2*alpha*r_i.*rho_ai + rho_ai.^2);
Lf_ar_inn = -asin(alpha) - asin(sqrt(1-alpha^2)*r_i./Lr_a_inn) + pi/2;  


initial_guess = 2.5; 
ra_solution = fsolve(@(ra) equation_to_solve(ra, alpha, r_w,r_out), initial_guess);
theta_ai = -asin(alpha*r_out./ra_solution) + asin(alpha) + sqrt(1/(alpha^2)-1) - sqrt((ra_solution./(alpha*r_out)).^2 -1);
 
[min_difference, index] = min(abs(Cf_ar_out - theta_ai));


alpha = l_y/sqrt(l_y^2 + (ra_solution*sin(theta_ai))^2);
rho_ai   = linspace(0,sqrt(l_y^2 + (ra_solution*sin(theta_ai))^2), 1000); 
r_i = r_out; LLr_a_out = sqrt(r_i^2 - 2*alpha*r_i.*rho_ai + rho_ai.^2); 
LLf_ar_out = asin(alpha) + asin(sqrt(1-alpha^2)*r_i./LLr_a_out) -pi/2;    

r_i = r_inn; LLr_a_inn = sqrt(r_i^2 + 2*alpha*r_i.*rho_ai + rho_ai.^2);
LLf_ar_inn = -asin(alpha) - asin(sqrt(1-alpha^2)*r_i./LLr_a_inn) + pi/2;  


theta = linspace(0, pi/4,1000); 
x_outer = r_out .* cos(theta);y_outer = r_out .* sin(theta); 
x_inner = r_inn .* cos(theta);y_inner = r_inn .* sin(theta); 


Cx_outer =  Cr_a .* cos(Cf_ar_out);  Cy_outer = Cr_a .* sin(Cf_ar_out);  
Cx_inner = Cr_a .* cos(Cf_ar_inn);   Cy_inner = Cr_a .* sin(Cf_ar_inn); 

CCx_outer = -Cr_a .* cos(Cf_ar_out) + r_w*2; CCy_outer = Cr_a .* sin(Cf_ar_out);  
CCr_a = sqrt(CCx_outer.^2+CCy_outer.^2); 
CCf_ar_out = atan2(CCy_outer,CCx_outer); 


Lx_outer = Lr_a_out .* cos(Lf_ar_out);    Ly_outer = Lr_a_out .* sin(Lf_ar_out);  
Lx_inner = Lr_a_inn .* cos(Lf_ar_inn);    Ly_inner = Lr_a_out .* sin(Lf_ar_inn);  


LLx_outer = LLr_a_out .* cos(LLf_ar_out);    LLy_outer = LLr_a_out .* sin(LLf_ar_out);  
LLx_inner = LLr_a_inn .* cos(LLf_ar_inn);    LLy_inner = LLr_a_out .* sin(LLf_ar_inn);  

figure(1)
width = 440; height = 270; 
set(gcf, 'Position', [10, 10, width, height]);

y_zero1 = zeros(1,1000); idx1 = 1; idx2 = 1000;
%%
h1 = plot(Cr_a,   Cf_ar_out,'color','k','linewidth',1.5,'lineStyle','-');hold on;    
     plot(Cr_a,   Cf_ar_inn,'color','k','linewidth',1.5,'lineStyle','-');hold on;    
h2 = plot(Lr_a_out, Lf_ar_out,'color','k','linewidth',1.5,'lineStyle','--');hold on; 
     plot(Lr_a_inn, Lf_ar_inn,'color','k','linewidth',1.5,'lineStyle','--');hold on; 
%%      
fill([Cr_a(idx1:1000), fliplr(Cr_a(idx1:1000))], [Cf_ar_out(idx1:1000), fliplr(y_zero1(idx1:1000))], [255/256 235/256 205/256], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
fill([Lr_a_inn(1:idx2), fliplr(Lr_a_inn(1:idx2))], [Lf_ar_inn(1:idx2), fliplr(y_zero1(1:idx2))], [255/256 235/256 205/256], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
%% 
h1 = plot(Cr_a,   Cf_ar_out,'color','k','linewidth',1.5,'lineStyle','-');hold on;    
     plot(Cr_a,   Cf_ar_inn,'color','k','linewidth',1.5,'lineStyle','-');hold on;    
%% 
h2 = plot(Lr_a_out, Lf_ar_out,'color','k','linewidth',1.5,'lineStyle','--');hold on; 
     plot(Lr_a_inn, Lf_ar_inn,'color','k','linewidth',1.5,'lineStyle','--');hold on; 
%% 
h3 = plot(CCr_a, CCf_ar_out,'color','k','linewidth',1.5,'lineStyle','-.');hold on;  
%% 
h6 = plot([r_inn,r_out], [0,0],'color',[210/256 105/256 30/256],'linewidth',1.5,'lineStyle','-');hold on;
     plot([r_out,r_out], [0,0.406],'color',[210/256 105/256 30/256],'linewidth',1.5,'lineStyle','-');hold on;
     plot([r_inn,r_inn], [0,0.406],'color',[210/256 105/256 30/256],'linewidth',1.5,'lineStyle','-');hold on;
%% 
h4 = plot(LLr_a_out, LLf_ar_out,'color','g','linewidth',1.5,'lineStyle','-');hold on; 
     plot(LLr_a_inn, LLf_ar_inn,'color','g','linewidth',1.5,'lineStyle','-');hold on; 
     plot([r_inn r_out], [0.0005 0.0005],'color','g','linewidth',0.5,'lineStyle','-');hold on;
%% 
h5 = plot(CCr_a(551:1000), CCf_ar_out(551:1000),'color','b','linewidth',1.5,'lineStyle','-');hold on; 
     plot(Cr_a(551:1000), Cf_ar_out(551:1000),'color','b','linewidth',1.5,'lineStyle','-');hold on; 
%%
%legend([h1,h3,h2,h4,h5,h6],'B-ALL', 'B-ALR','B-SL','BCH','BERG','TL','FontSize',10,'Interpreter','latex','Location', 'north','NumColumns',4);
ylabel('$\theta_{ai}$','Interpreter','latex','FontSize',12);hold on; 
xlabel('$r_a$','Interpreter','latex','FontSize',12);hold on; 
axis([r_inn r_out 0 0.406]);