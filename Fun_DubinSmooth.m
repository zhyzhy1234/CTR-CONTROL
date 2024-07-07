%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project name: 
%Programer   : 
%Finish date : 
%Records     :        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = Fun_DubinSmooth(waypoint_s,waypoint_i,waypoint_t,r,l_y)
    
    vector_sq = waypoint_i - waypoint_s;
    vector_kt = waypoint_t - waypoint_i;

    [rho_sq,fai_sq] = Fun_CalDistanAndAngle(waypoint_i,waypoint_s); % waypoint_p - waypoint_s    0 - 2*pi
    [rho_kt,fai_kt] = Fun_CalDistanAndAngle(waypoint_t,waypoint_i); % waypoint_t - waypoint_k
    
    bias = fai_sq - fai_kt;
    if bias >= pi
        fai_sq = fai_sq - 2*pi;
    elseif bias <= -pi
        fai_sq = fai_sq + 2*pi;
    end

    mu = sign(fai_sq - fai_kt);

    theta_c = acos( dot(vector_sq,vector_kt)/(rho_sq*rho_kt) ); 

    d = abs(r*tan(theta_c/2));

    waypoint_q = waypoint_i - d*[cos(fai_sq), sin(fai_sq)]; 
    waypoint_k = waypoint_i + d*[cos(fai_kt), sin(fai_kt)];

    waypoint_o = waypoint_q + r * mu * [sin(fai_sq), -cos(fai_sq)];
    
    r1 = r - mu*l_y(2);
    waypoint_s1 = waypoint_s + ([cos(fai_sq) sin(fai_sq);sin(fai_sq) -cos(fai_sq)]*l_y(1:2)')';
    waypoint_t1 = waypoint_t + ([cos(fai_kt) sin(fai_kt);sin(fai_kt) -cos(fai_kt)]*l_y(1:2)')';
    waypoint_q1 = waypoint_q + ([cos(fai_sq) sin(fai_sq);sin(fai_sq) -cos(fai_sq)]*l_y(1:2)')';
    waypoint_k1 = waypoint_k + ([cos(fai_kt) sin(fai_kt);sin(fai_kt) -cos(fai_kt)]*l_y(1:2)')';
    
    [rho_sq,fai_sq] = Fun_CalDistanAndAngle(waypoint_q1,waypoint_s1); % waypoint_p - waypoint_s

    w_1 = rho_sq;     
    w_2 = theta_c*r1; 
    
    Sum_w = w_1 + w_2;
    
    result = [waypoint_s1,waypoint_q1,waypoint_k1,waypoint_t1,waypoint_o,theta_c,r1,mu,Sum_w,w_1];
end

