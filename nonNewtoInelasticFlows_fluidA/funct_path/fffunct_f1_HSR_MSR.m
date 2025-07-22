function f1_HSR_MSR = fffunct_f1_HSR_MSR(gradpTPL_aux, DgradpTPL_Dz_aux, Raux, theta, R0aux, RInfaux, mu0, K, n, muInf, alpha, beta)
    
    
    f1_HSR_MSR = R0aux^(alpha+2);
    f1_HSR_MSR = f1_HSR_MSR*DgradpTPL_Dz_aux/gradpTPL_aux;
    f1_HSR_MSR = f1_HSR_MSR*( (-gradpTPL_aux/(2*K))^(1/n) );
    f1_HSR_MSR = f1_HSR_MSR*(-1/(n*alpha*beta));
    
    DA1Dz_aux = fffunct_DA1Dz(gradpTPL_aux, DgradpTPL_Dz_aux, Raux, theta, RInfaux, K, n, muInf, alpha);
    f1_HSR_MSR = f1_HSR_MSR + DA1Dz_aux*(R0aux^2)/2;
    
    radvel_HSR_LSR_in_R0 = ffunct_radvel_HSR_LSR(R0aux, gradpTPL_aux, DgradpTPL_Dz_aux, Raux, theta, R0aux, RInfaux, mu0, K, n, muInf, alpha);                                
    f1_HSR_MSR = f1_HSR_MSR + R0aux*radvel_HSR_LSR_in_R0;
    
            

end