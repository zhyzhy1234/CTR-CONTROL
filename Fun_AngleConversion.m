function alpha = Fun_AngleConversion(X) 
    x = X(1);
    y = X(2);
    
    if y < 0
        alpha = 2*pi + atan2(real(y),real(x)); 
    else 
        alpha = atan2(real(y),real(x));     
    end    
end