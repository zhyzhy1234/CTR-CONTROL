%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project name: Path smoothing
%Programer   : zhy
%Finish date : 
%Records     : Given waypoints, find the smoothed path (arc length accuracy: 0.01m)

%               -----------  eta_el  eta_point(m,:,2)
%               -.-.-.-.-.-  eta_0   eta_point(m,:,1)
%               -----------  eta_er  eta_point(m,:,3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
close all;


obstacleVertices_3 = [0, 6; 0, 7; 3, 7; 3, 6];
obstacleVertices_4 = [3, 5; 3, 6; 4, 6; 4, 5];
obstacleVertices_5 = [1, 2.5; 1, 3.5; 2, 3.5; 2, 2.5];

obstacleVertices_2 = [5.5, 8; 5.5, 10; 7.5, 10; 7.5, 7; 6.5, 7; 6.5, 8];
obstacleVertices_6 = [5.5, 4+0.5; 5.5, 6+0.5; 7.5, 6+0.5; 7.5, 5+0.5; 6.5, 5+0.5; 6.5, 4+0.5];
obstacleVertices_7 = [5.5, 0.3; 3.5, 0.3; 3.5, 1.3; 4.5, 1.3; 4.5, 2.3; 6.5, 2.3; 6.5, 1.3; 5.5, 1.3];

obstacleVertices_8 = [7.5, 2; 8.5, 2; 8.5, 3; 7.5, 3];
obstacleVertices_1 = [9, 3.5; 10, 3.5; 10, 4.5; 9, 4.5];

% A*  waypoints
p_0 = [0.5 1.5; ...
       2.5 1.5; ...
       4.5 3.5; ...
       6.5 3.5; ...
       8.5 4.5; ...
       8.5 7.5; ...
       9.5 8.5
];

epsilon = 0.6;  
r_w = 2.5;  

p_w = p_0;
n_w = size(p_w,1);

L = [0 0;...        
     0 -epsilon;... 
     0 epsilon];   
   
start_0 = 0; start_r = 0; start_l = 0;                                         
end_0   = 0; end_r   = 0; end_l   = 0;                                         
l_w0    = zeros(n_w-2,2);    l_wr = zeros(n_w-2,2);    l_wl = zeros(n_w-2,2); 
eta_s0  = zeros(n_w-2,15); eta_sr = zeros(n_w-2,15); eta_sl = zeros(n_w-2,15); 


for i = 1:n_w - 2 
    
    if i > 1
        p_w(i,1:2) = eta_s0(i-1,5:6);     
        eta_s0(i,1:2) = eta_s0(i-1,5:6);  
        eta_sr(i,1:2) = eta_sr(i-1,5:6);  
        eta_sl(i,1:2) = eta_sl(i-1,5:6); 
    end
    
    eta_s0(i,:) = Fun_DubinSmooth(p_w(i,:),p_w(i+1,:),p_w(i+2,:),r_w,L(1,:));
    eta_sr(i,:) = Fun_DubinSmooth(p_w(i,:),p_w(i+1,:),p_w(i+2,:),r_w,L(3,:));
    eta_sl(i,:) = Fun_DubinSmooth(p_w(i,:),p_w(i+1,:),p_w(i+2,:),r_w,L(2,:)); 
    

    if i == n_w - 2 
        end_0 = start_0 + eta_s0(i,14) + norm( eta_s0(i,5:6)-p_0(n_w,:),2); l_w0(i,:) = [start_0, end_0];
        end_r = start_r + eta_sr(i,14) + norm( eta_sr(i,5:6)-p_0(n_w,:),2); l_wr(i,:) = [start_r, end_r];
        end_l = start_l + eta_sl(i,14) + norm( eta_sl(i,5:6)-p_0(n_w,:),2); l_wl(i,:) = [start_l, end_l];
    else
        end_0 = start_0 + eta_s0(i,14); l_w0(i,:) = [start_0, end_0]; start_0 = end_0;
        end_r = start_r + eta_sr(i,14); l_wr(i,:) = [start_r, end_r]; start_r = end_r;
        end_l = start_l + eta_sl(i,14); l_wl(i,:) = [start_l, end_l]; start_l = end_l;
    end   
    
end


Iend_0 = floor(end_0*100);  
varpi = linspace(0,end_0, Iend_0); 
for m = 1:Iend_0
    number_w0 = varpi(m); 
    for i = 1:n_w-2
        if number_w0 >= l_w0(i, 1) && number_w0 <= l_w0(i, 2)
            index_w0 = i; 
            break;
        end 
    end
    eta_l0(m,:) = Fun_SfToCartesian(eta_s0(index_w0,:),number_w0 - l_w0(index_w0,1),[0 0],index_w0);   
end


Iend_l = floor(end_l*100);  
varpi = linspace(0,end_l, Iend_l); 
for m = 1:Iend_l
    number_wl = varpi(m); 
    for i = 1:n_w-2
        if number_wl >= l_wl(i, 1) && number_wl <= l_wl(i, 2)
            index_wl = i;  
            break;
        end 
    end
    eta_ll(m,:) = Fun_SfToCartesian(eta_sl(index_wl,:),number_wl - l_wl(index_wl,1),[0 0],index_wl);  
end


Iend_r = floor(end_r*100);  
varpi = linspace(0,end_r, Iend_r); 
for m = 1:Iend_r
    number_wr = varpi(m); 
    for i = 1:n_w-2
        if number_wr >= l_wr(i, 1) && number_wr <= l_wr(i, 2)
            index_wr = i;  
            break;
        end 
    end
    eta_lr(m,:) = Fun_SfToCartesian(eta_sr(index_wr,:),number_wr - l_wr(index_wr,1),[0 0],index_wr);   
end

save('data_p_0.mat', 'p_0');
save('data_n_w.mat', 'n_w'); save('data_epsilon.mat', 'epsilon'); save('data_r_w.mat', 'r_w');save('data_L.mat', 'L'); 
save('data_eta_s0.mat', 'eta_s0');save('data_eta_l0.mat', 'eta_l0'); save('data_l_w0.mat', 'l_w0'); 
save('data_eta_sl.mat', 'eta_sl');save('data_eta_ll.mat', 'eta_ll'); save('data_l_wl.mat', 'l_wl');
save('data_eta_sr.mat', 'eta_sr');save('data_eta_lr.mat', 'eta_lr'); save('data_l_wr.mat', 'l_wr');



figure(2)

width = 500;  height = 400; 
set(gcf, 'Position', [10, 10, width, height]);

obsColor = [0, 0,0];

patch('Vertices', obstacleVertices_1, 'Faces', 1:size(obstacleVertices_1, 1), 'FaceColor', obsColor, 'EdgeColor', 'black');hold on;
patch('Vertices', obstacleVertices_2, 'Faces', 1:size(obstacleVertices_2, 1), 'FaceColor', obsColor, 'EdgeColor', 'black');hold on;
patch('Vertices', obstacleVertices_3, 'Faces', 1:size(obstacleVertices_3, 1), 'FaceColor', obsColor, 'EdgeColor', 'black');
patch('Vertices', obstacleVertices_4, 'Faces', 1:size(obstacleVertices_4, 1), 'FaceColor', obsColor, 'EdgeColor', 'black');
patch('Vertices', obstacleVertices_5, 'Faces', 1:size(obstacleVertices_5, 1), 'FaceColor', obsColor, 'EdgeColor', 'black');
patch('Vertices', obstacleVertices_6, 'Faces', 1:size(obstacleVertices_6, 1), 'FaceColor', obsColor, 'EdgeColor', 'black');
patch('Vertices', obstacleVertices_7, 'Faces', 1:size(obstacleVertices_7, 1), 'FaceColor', obsColor, 'EdgeColor', 'black');
patch('Vertices', obstacleVertices_8, 'Faces', 1:size(obstacleVertices_8, 1), 'FaceColor', obsColor, 'EdgeColor', 'black');

h1 = plot(p_0(:,1),p_0(:,2)); hold on;

scatter(9.5, 8.5, 'SizeData', 40, 'MarkerEdgeColor', 'blue', 'MarkerFaceColor', 'red');
center = [9.5, 8.5];
viscircles(center, epsilon, 'EdgeColor', 'green', 'LineStyle', '--','LineWidth', 1);


h2 = plot(eta_l0(:,1),eta_l0(:,2),'-.');hold on; 
h3 = plot(eta_ll(:,1),eta_ll(:,2),'--');hold on;
h4 = plot(eta_lr(:,1),eta_lr(:,2),'--');hold on;
legend([h1,h2,h3,h4],'$ \eta_L^0 $','$ \eta_S^0 $','$ \eta_S^l $','$ \eta_S^r $','Interpreter','latex','Location', 'northwest','FontSize',12, 'NumColumns', 2); %
box on; axis([0 10.5 0 10.5]);  %axis equal;         




