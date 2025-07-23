function AxialVelTPL = funct_AxialVel_TPL(z, gradp, R, R0, mu0, K, n)
    % this function returns the axial velocity field along the axis for a TPL fluid
  
    AxialVelTPL = z*0;
    
    alpha = (n+1)/n;
    
    for i=1:size(AxialVelTPL,2)
        if     R0(1,i) >= 100 %LSR regime
             AxialVelTPL(1,i) = -gradp(1,i).*( R(1,i).^2 )./(4*mu0);
             
        elseif R0(1,i) < 100  %MSR regime
             A0 = -gradp(1,i).*( ( (R0(1,i)./100).*R(1,i)).^2 )/(4*mu0) + ( ( -gradp(1,i)./(2*K) ).^(1/n) ).*(  ( R(1,i).^(alpha) ) - ( ( (R0(1,i)./100).*R(1,i)).^(alpha) ) )./alpha;
             AxialVelTPL(1,i) =  A0;
        end
    end
        
    
end