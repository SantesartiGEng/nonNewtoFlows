function visc = visc_modTPL_subcase1(p, Gamma)

 %%%- TPL Model  - sub case1
    
    gamma0_TPL = 1/p(2); 
    K = p(1).*( gamma0_TPL ).^(1-p(3));
    
    x_PL = [K p(3)];
    
    if     Gamma <= gamma0_TPL
            visc = p(1);
    elseif (Gamma > gamma0_TPL) 
            visc = visc_modPL(  x_PL, Gamma );
    end
    
end