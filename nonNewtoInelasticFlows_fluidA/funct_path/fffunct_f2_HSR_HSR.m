function f2_HSR_HSR = fffunct_f2_HSR_HSR(gradpTPL_aux, DgradpTPL_Dz_aux, Raux, theta, R0aux, RInfaux, mu0, K, n, muInf, alpha, beta)
    
    
    f2_HSR_HSR = RInfaux^2*theta*gradpTPL_aux*Raux/(4*muInf);
    
    f2_HSR_HSR = f2_HSR_HSR - DgradpTPL_Dz_aux*RInfaux^2*(2*Raux^2-RInfaux^2)/(16*muInf);
    
    radvel_HSR_MSR_in_RInf = ffunct_radvel_HSR_MSR(RInfaux, gradpTPL_aux, DgradpTPL_Dz_aux, Raux, theta, R0aux, RInfaux, mu0, K, n, muInf, alpha, beta);
    f2_HSR_HSR = f2_HSR_HSR + RInfaux*radvel_HSR_MSR_in_RInf;
    
            

end