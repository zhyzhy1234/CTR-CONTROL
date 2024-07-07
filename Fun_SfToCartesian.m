%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project name: 
%Programer   : 
%Finish date :
%Records     :                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = Fun_SfToCartesian(p_s,w,ref_point0,index_w)
    
    waypoint_s = p_s(1:2);
    waypoint_q = p_s(3:4);
    waypoint_k = p_s(5:6);
    waypoint_t = p_s(7:8);
    waypoint_o = p_s(9:10);
    theta_c    = p_s(11);
    R = p_s(12);
    epsilon  = p_s(13); 

    
    [rho_sq,fai_sq] = Fun_CalDistanAndAngle(waypoint_q,waypoint_s); % waypoint_p - waypoint_s
    [rho_kt,fai_kt] = Fun_CalDistanAndAngle(waypoint_t,waypoint_k); % waypoint_t - waypoint_k

%     vector_cp = waypoint_q - waypoint_o;
%     fai_cp = atan2(vector_cp(2), vector_cp(1));
    
    w_1 = rho_sq;     
    w_2 = theta_c*R;  

    if w <= w_1    
        ref_point = [waypoint_s + w.*[cos(fai_sq) sin(fai_sq)],fai_sq,0,waypoint_o,index_w,w_1,epsilon,0];                                      % 01 waypoint_s -- waypoint_p
    elseif w_1 < w && w < w_1 + w_2
%         theta_circle = fai_cp - flag * (w - w_1)/R;
        theta_circle = fai_sq + epsilon*( pi/2 - (w - w_1)/R );
        ref_point = [waypoint_o + R.*[cos(theta_circle) sin(theta_circle)],theta_circle - epsilon * pi/2, 1/R,waypoint_o,index_w,w_1,epsilon,w-w_1]; % 02 waypoint_c -- waypoint_p -- waypoint_k 
    elseif w >= w_1 + w_2 
        ref_point = [waypoint_k + (w - w_1 - w_2).*[cos(fai_kt) sin(fai_kt)],fai_kt, 0,waypoint_o,index_w,w_1,epsilon,0];                        % 03 waypoint_k -- waypoint_t 
    else
        ref_point = ref_point0;
    end
    result = ref_point;
end

