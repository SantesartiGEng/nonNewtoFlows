function radvel_radius_HSR_HSR = ffunct_radvel_HSR_HSR(r, gradpTPL_aux, DgradpTPL_Dz_aux, Raux, theta, R0aux, RInfaux, mu0, K, n, muInf, alpha, beta)
     
     radvel_radius_HSR_HSR = DgradpTPL_Dz_aux*r/(16*muInf);
     radvel_radius_HSR_HSR = radvel_radius_HSR_HSR*(2*Raux^2 - r^2);
     
     radvel_radius_HSR_HSR = radvel_radius_HSR_HSR - r*theta*gradpTPL_aux*Raux/(4*muInf);
     
     f2_HSR_HSR = fffunct_f2_HSR_HSR(gradpTPL_aux, DgradpTPL_Dz_aux, Raux, theta, R0aux, RInfaux, mu0, K, n, muInf, alpha, beta);
     radvel_radius_HSR_HSR = radvel_radius_HSR_HSR + f2_HSR_HSR/r;
     
  

end