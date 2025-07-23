function radvel_radius_MSR_MSR = ffunct_radvel_MSR_MSR(r, gradpTPL_aux, DgradpTPL_Dz_aux, Raux, theta, R0aux, mu0, K, n, alpha)
     
     radvel_radius_MSR_MSR = (Raux^alpha)/2 - (r^alpha)/(alpha+2);
     radvel_radius_MSR_MSR = radvel_radius_MSR_MSR*DgradpTPL_Dz_aux/(gradpTPL_aux*alpha*n);
     radvel_radius_MSR_MSR = radvel_radius_MSR_MSR - (Raux^(alpha-1))*theta/2;
     radvel_radius_MSR_MSR = radvel_radius_MSR_MSR*( (-gradpTPL_aux/(2*K))^(1/n) );
     radvel_radius_MSR_MSR = -r*radvel_radius_MSR_MSR;
     
     f1_MSR = fffunct_f1_MSR(gradpTPL_aux, DgradpTPL_Dz_aux, Raux, theta, R0aux, mu0, K, n, alpha);
     radvel_radius_MSR_MSR = radvel_radius_MSR_MSR + f1_MSR/r;
     
  

end