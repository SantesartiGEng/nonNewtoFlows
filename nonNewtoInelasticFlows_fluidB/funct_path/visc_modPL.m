function visc2 = visc_modPL(p,Gamma)
%%%- Power law model 
    visc2 = p(1).*(Gamma.^(p(2)-1));

end