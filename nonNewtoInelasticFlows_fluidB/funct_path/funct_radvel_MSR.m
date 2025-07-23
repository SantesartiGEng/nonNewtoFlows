function [radvel, radiusvect] = funct_radvel_MSR(index_z_axis, nodes_ray, R, theta, gradpTPL, Dgradp_TPL_Dz, mu0, tau0, K, n, alpha)
   %%%this function computes the radial velocity profile of a TPL fluid in
   %%%the MSR regime in a specific index_z_axis coordinate
    
    radvel = 0*ones(1,nodes_ray);
    Raux =  R(1,index_z_axis);
    radiusvect = linspace(0, Raux, nodes_ray);
    
    gradpTPL_aux = gradpTPL(1,index_z_axis);
    
    absgradpTPL = abs(gradpTPL);
    absgradpTPL_aux = absgradpTPL(1,index_z_axis);
    
    DgradpTPLDz_aux = Dgradp_TPL_Dz(1,index_z_axis);
    
    R0aux = funct_R0(absgradpTPL_aux, tau0);
    
    for i=1:nodes_ray
        
        if     radiusvect(1,i) <= R0aux
               radvel(1,i) = ffunct_radvel_MSR_LSR(radiusvect(1,i), gradpTPL_aux, DgradpTPLDz_aux, Raux, theta, R0aux, mu0, K, n, alpha);
                           
        elseif radiusvect(1,i) > R0aux
               radvel(1,i) = ffunct_radvel_MSR_MSR(radiusvect(1,i), gradpTPL_aux, DgradpTPLDz_aux, Raux, theta, R0aux, mu0, K, n, alpha);

        end
    end

end