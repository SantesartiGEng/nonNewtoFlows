function radvel_radius_HSR_LSR = ffunct_radvel_HSR_LSR(r, gradpTPL_aux, DgradpTPL_Dz_aux, Raux, theta, R0aux, RInfaux, mu0, K, n, muInf, alpha)
     
     DA0Dz= fffunct_DA0Dz(gradpTPL_aux, DgradpTPL_Dz_aux, Raux, theta, R0aux, RInfaux, mu0, K, n, muInf, alpha);               
     radvel_radius_HSR_LSR = DA0Dz;
     radvel_radius_HSR_LSR = radvel_radius_HSR_LSR + DgradpTPL_Dz_aux*(r^2)/(8*mu0); 
     radvel_radius_HSR_LSR = -r/2*radvel_radius_HSR_LSR;

end