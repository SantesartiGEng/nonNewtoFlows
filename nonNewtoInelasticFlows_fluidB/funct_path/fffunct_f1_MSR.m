function f1_MSR = fffunct_f1_MSR(gradpTPL_aux, DgradpTPL_Dz_aux, Raux, theta, R0aux, mu0, K, n, alpha)
    
    
    f1_MSR = (Raux^alpha)/2 - (R0aux^alpha)/(alpha+2);
    f1_MSR = f1_MSR*DgradpTPL_Dz_aux/(gradpTPL_aux*alpha*n);
    f1_MSR = f1_MSR - (Raux^(alpha-1))*theta/2;
    f1_MSR = f1_MSR*( (-gradpTPL_aux/(2*K))^(1/n) );
    f1_MSR = +(R0aux^2)*f1_MSR;
   
    radvel_MSR_LSR_in_R0 = ffunct_radvel_MSR_LSR(R0aux, gradpTPL_aux, DgradpTPL_Dz_aux, Raux, theta, R0aux, mu0, K, n, alpha);                                
     
    f1_MSR = f1_MSR + R0aux*radvel_MSR_LSR_in_R0;
    
            

end