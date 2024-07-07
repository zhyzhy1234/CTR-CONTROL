%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project name:  
%Programer   : 
%Finish date : 
%Records     : 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rho,theta] = Fun_CalDistanAndAngle(newSt_d1,newSt_d2)

    q_d  = newSt_d1 - newSt_d2;               
    rho  = norm(q_d,2);
    
    x = q_d(1);
    y = q_d(2);
    
    if y < 0
        theta = 2*pi + atan2(real(y),real(x)); 
    else 
        theta = atan2(real(y),real(x));        
    end
    
end