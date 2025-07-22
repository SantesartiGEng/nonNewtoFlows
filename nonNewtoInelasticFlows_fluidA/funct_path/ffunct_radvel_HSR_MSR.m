function radvel_radius_HSR_MSR = ffunct_radvel_HSR_MSR(r, gradpTPL_aux, DgradpTPL_Dz_aux, Raux, theta, R0aux, RInfaux, mu0, K, n, muInf, alpha, beta)
     
     radvel_radius_HSR_MSR = r^(alpha+1);
     radvel_radius_HSR_MSR = radvel_radius_HSR_MSR*DgradpTPL_Dz_aux/gradpTPL_aux;
     radvel_radius_HSR_MSR = radvel_radius_HSR_MSR*( (-gradpTPL_aux/(2*K))^(1/n) );
     radvel_radius_HSR_MSR = radvel_radius_HSR_MSR/(n*alpha*beta);
     
     DA1Dz = fffunct_DA1Dz(gradpTPL_aux, DgradpTPL_Dz_aux, Raux, theta, RInfaux, K, n, muInf, alpha);
     radvel_radius_HSR_MSR =  radvel_radius_HSR_MSR - DA1Dz*r/2;
     
     f1_HSR_MSR = fffunct_f1_HSR_MSR(gradpTPL_aux, DgradpTPL_Dz_aux, Raux, theta, R0aux, RInfaux, mu0, K, n, muInf, alpha, beta);
     radvel_radius_HSR_MSR = radvel_radius_HSR_MSR + f1_HSR_MSR/r;
     
  

end