function visc = visc_modTPL(p, Gamma)
% % %%%- truncated power-law model 
    
    gamma0_TPL = 1/p(2); 
    gammaInf_TPL = 1/p(3);
    K = p(1).*( (gamma0_TPL).^(1-p(4)) );
    muInf = p(1).*( (gamma0_TPL/gammaInf_TPL ).^(1-p(4)) );
    
    x_PL = [K p(4)];
    
    if     Gamma <= gamma0_TPL
            visc = p(1);
    elseif (Gamma > gamma0_TPL) & (Gamma <= gammaInf_TPL) 
            visc = visc_modPL(  x_PL, Gamma );
    elseif Gamma > gammaInf_TPL   
            visc = muInf;
    end

    
end