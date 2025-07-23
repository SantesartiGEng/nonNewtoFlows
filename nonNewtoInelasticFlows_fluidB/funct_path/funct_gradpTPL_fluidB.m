function gradpTPL = funct_gradpTPL_fluidB(gradpNewto_mu0, gradpPowlaw, Q, R, mu0, K, n, tau0)
    % This function computes the gradp of the fluid B approximated by the TPL model flowing in a conical
    % nozzle following the Algorithm 1 reported in the reference article: 
    % Santesarti, G., Marino, M., Viola, F., Verzicco, R., & Vairo, G. (2025). A Quasi-Analytical Solution for "Carreau-Yasuda-like" Shear-thinning Fluids Flowing in Slightly Tapered Pipes. arXiv preprint 2502.14991 at https://doi.org/10.48550/arXiv.2502.14991
    
    %%%--- General arrays
    gradpTPL = gradpNewto_mu0*0;
    
    Abs_gradpNewto_mu0 = abs(gradpNewto_mu0);
    Abs_gradpPowlaw = abs(gradpPowlaw);
    
    %% gradp evaluation at the inlet
    mingradp = min( [Abs_gradpNewto_mu0(1) Abs_gradpPowlaw(1)] );
    
    if mingradp <= 2*tau0./R(1)
        gradpTPL(1) = gradpNewto_mu0(1);
    elseif mingradp > 2*tau0./R(1)
        gradpTPL(1) = ffunct_Qgradp_fzero(mingradp, Q, R(1), mu0, K, n, tau0);
    end
    
    
    %% gradp evaluation along the axis
    for i=2:( size(R,2) ) 
        
        back_gradp = abs(gradpTPL(i-1));
        
        if  abs( gradpTPL(i-1) ) > 2*tau0./R(i-1) % if in the previous section there is a MSR regime
                gradpTPL(i) = ffunct_Qgradp_fzero(back_gradp, Q, R(i), mu0, K, n, tau0);
            
        elseif gradpTPL(i-1) <= 2*tau0./R(i-1) % if in the previous section there is a LSR regime
                   gradpTPL(i) = gradpNewto_mu0(i); % an attempt to assume the LSR regime in the current section

                   if abs( gradpTPL(i) ) > 2*tau0./R(i) % check to verify if there is a MSR regime
                      gradpTPL(i) = ffunct_Qgradp_fzero(back_gradp, Q, R(i), mu0, K, n, tau0);
                   end
        end
        
        
    end
    
    
    
    
end