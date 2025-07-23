function AxVelTPL_radius = funct_AxVelTPL_radius(gradp, r, R, R0, mu0, K, n)
   % this function returns the axial velocity field along the radius for a
   % TPL fluid 
  
    AxVelTPL_radius = r*0;
    
    alpha = (n+1)/n;
   
  
    if R0 >= 100 
        for i=1:size(r,2)
            AxVelTPL_radius(1,i) = -gradp.*( R.^2 - r(i).^2 )./(4*mu0);
        end
        
    elseif R0 < 100 
        A0 = -gradp.*( ( (R0./100).*R).^2 )/(4*mu0) + ( ( -gradp./(2*K) ).^(1/n) ).*(  ( R.^(alpha) ) - ( ( (R0./100).*R).^(alpha) ) )./alpha;
        
        for i=1:size(r,2)
            if r(i) <= ( (R0./100).*R)
                AxVelTPL_radius(i) = (gradp.*( r(i).^2 )./(4*mu0)) + A0;
               
            elseif r(i) > ( (R0./100).*R)
                
                AxVelTPL_radius(i) = ( ( -gradp./(2*K) ).^(1/n) );
                AxVelTPL_radius(i) = AxVelTPL_radius(i).*( (R.^alpha) - ( r(i).^alpha ) )./alpha;
            end
        end
        
    end
    
end