function gradpTPL = funct_gradpTPL_fluidA(gradpNewto_mu0, gradpPowlaw, gradpNewto_muInf, Q, R, mu0, tau0, K, n, muInf, tauInf)
    % This function computes the gradp of the fluid A approximated by the TPL model flowing in a conical
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
        
    elseif (mingradp > 2*tau0./R(1) && mingradp <= 2*tauInf./R(1))
        gradpTPL(1) = ffunct_Qgradp_MSR_fzero(mingradp, Q, R(1), mu0, tau0, K, n);
        
    elseif mingradp > 2*tauInf./R(1)
        gradpTPL(1) = ffunct_Qgradp_HSR_fzero(abs(gradpNewto_muInf(1,1)), Q, R(1), mu0, tau0, K, n, muInf, tauInf);
        
    end
   
    
    
    %% gradp evaluation along the axis
    for i=2:( size(R,2) ) 
        
        back_gradp = abs(gradpTPL(i-1));
        
        if       back_gradp <= 2*tau0./R(i-1) % if in the previous section there is a LSR regime
            
                    gradpTPL(i) = gradpNewto_mu0(i);% an attempt to assume the LSR regime in the current section

                    if abs( gradpTPL(i) ) > 2*tau0./R(i) % check to verify if there is a MSR regime
                        
                        gradpTPL(i) = ffunct_Qgradp_MSR_fzero(back_gradp, Q, R(i), mu0, tau0, K, n);
                        
                    end
            
        elseif   (back_gradp > 2*tau0./R(i-1) && back_gradp <= 2*tauInf./R(i-1) )% if in the previous section there is a MSR regime
            
                    gradpTPL(i) = ffunct_Qgradp_MSR_fzero(back_gradp, Q, R(i), mu0, tau0, K, n);% an attempt to assume the MSR regime in the current section
                   
                    if abs(gradpTPL(i)) > 2*tauInf./R(i) % check to verify if there is a HSR regime
                        
                        gradpTPL(i) = ffunct_Qgradp_HSR_fzero(back_gradp, Q, R(i), mu0, tau0, K, n, muInf, tauInf);
                        
                    end
                    
        elseif   back_gradp > 2*tauInf./R(i-1) % if in the previous section there is a HSR regime
    
                    gradpTPL(i) = ffunct_Qgradp_HSR_fzero(back_gradp, Q, R(i), mu0, tau0, K, n, muInf, tauInf);
                    
        end
        
        
    end
    
    
    
    
end