function visc = visc_modSRB(p,Gamma)

% %%%- SRB model
    visc = p(1).*(  ( 1+ (Gamma.*p(3)).^p(5)  )./( 1+ (Gamma.*p(2)).^p(5) )  ).^( ( 1-p(4) )./p(5) );
    
end