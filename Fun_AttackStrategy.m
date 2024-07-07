%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Project name:  attacker policy
%Programer   : zhy
%Finish date : 
%Records     : 
%                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Exp_a = Fun_AttackStrategy(newSt_d,newSt_a)

    rhod_a = 0.3;
    c_esc = 1; 
    k = 0;
    
    q_ad(:,:) = reshape( (- newSt_d(1,5:6,:) + newSt_a(1,5:6)),[2,3])'; 
    rho_ad = vecnorm(q_ad,2,2)';                                   

    Neighbor_a = zeros(1,2);
    for i = 1:3  
        if rho_ad(i) < rhod_a                                        
            Neighbor_a = q_ad(i,:)/rho_ad(i) + Neighbor_a;           
            k = 1;
        end
    end
    
    if k == 1   
        Exp_a =  Fun_AngleConversion(c_esc.*Neighbor_a);
    else
        Exp_a =  pi*rand; 
%         Exp_a(2) =  0.2*rand; 
    end

     
end