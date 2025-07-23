function radvel_radius_MSR_LSR = ffunct_radvel_MSR_LSR(r, gradpTPL_aux, DgradpTPL_Dz_aux, Raux, theta, R0aux, mu0, K, n, alpha)
     
     DA0Dz_MSR = fffunct_DA0Dz_MSR(gradpTPL_aux, DgradpTPL_Dz_aux, Raux, theta, R0aux, mu0, K, n, alpha);               
     radvel_radius_MSR_LSR = DA0Dz_MSR;
     radvel_radius_MSR_LSR = radvel_radius_MSR_LSR + DgradpTPL_Dz_aux*(r^2)/(8*mu0); 
     radvel_radius_MSR_LSR = -r/2*radvel_radius_MSR_LSR;

end