function diff = equation_to_solve(ra, alpha, r_w,r_out)
    diff = -asin(alpha*r_out./ra) + asin(alpha) + sqrt(1/(alpha^2)-1) - sqrt((ra./(alpha*r_out)).^2 -1)  - acos(r_w/ra);
end