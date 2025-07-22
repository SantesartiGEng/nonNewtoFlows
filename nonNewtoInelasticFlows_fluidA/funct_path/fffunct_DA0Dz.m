function DA0Dz= fffunct_DA0Dz(gradpTPL_aux, DgradpTPL_Dz_aux, Raux, theta, R0aux, RInfaux, mu0, K, n, muInf, alpha)  
    
    DA0Dz = Raux*theta*gradpTPL_aux/(2*muInf);
    
    DA0Dz = DA0Dz - DgradpTPL_Dz_aux/4*( (Raux^2+RInfaux^2)/muInf - R0aux^2/mu0 ); 
    
    DA0Dz = DA0Dz - DgradpTPL_Dz_aux/(gradpTPL_aux*alpha)*( (-gradpTPL_aux/(2*K))^(1/n) )*( RInfaux^alpha - R0aux^alpha);
   
    
end