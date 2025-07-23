function DA0Dz_MSR= fffunct_DA0Dz_MSR(gradpTPL_aux, DgradpTPL_Dz_aux, Raux, theta, R0aux, mu0, K, n, alpha)  
    
    
    DA0Dz_MSR = -( (-gradpTPL_aux/(2*K))^(1/n) )*Raux^(alpha-1)*theta;
    
    DA0Dz_MSR = DA0Dz_MSR + DgradpTPL_Dz_aux*(R0aux^2)/(4*mu0); 
    
    DA0Dz_MSR = DA0Dz_MSR + DgradpTPL_Dz_aux/(gradpTPL_aux*alpha)*( (-gradpTPL_aux/(2*K))^(1/n) )*( (Raux^alpha)/n + R0aux^alpha);
   
    
end