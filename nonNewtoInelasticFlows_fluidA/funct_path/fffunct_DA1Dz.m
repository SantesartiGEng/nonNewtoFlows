function DA1Dz= fffunct_DA1Dz(gradpTPL_aux, DgradpTPL_Dz_aux, Raux, theta, RInfaux, K, n, muInf, alpha)  
    
    DA1Dz = RInfaux^alpha;
    DA1Dz = DA1Dz*( (-gradpTPL_aux/(2*K))^(1/n) );
    DA1Dz = DA1Dz*(-DgradpTPL_Dz_aux/(gradpTPL_aux*alpha));

    DA1Dz = DA1Dz + gradpTPL_aux*Raux*theta/(2*muInf);
    
    DA1Dz = DA1Dz - DgradpTPL_Dz_aux/(4*muInf)*(RInfaux^2+Raux^2);

end