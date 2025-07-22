function AxVelTPL_radius = funct_AxVelTPL_radius_HSR(gradp, r, R0, RInf, mu0, K, n, muInf)
    % this function returns the axial velocity field along the radius for a
    % TPL fluid in HSR regime
    
    R = r(1,end);
    alpha = (n+1)/n;
    
    AxVelTPL_radius = r*0;
    
    A0 = -gradp./4.*( ( (R0./100.*R).^2 )/mu0 + ( R.^2 - (RInf./100.*R).^2 )/muInf ) + ( ( -gradp./(2*K) ).^(1/n) ).*(  (RInf./100.*R).^alpha  - (R0./100.*R).^alpha )./alpha;
    
    A1 = -gradp./(4*muInf).*( R.^2 - (RInf./100.*R).^2 ) + ( ( -gradp./(2*K) ).^(1/n) )./alpha.*( (RInf./100.*R).^alpha ) ;
             
  
 
    for i=1:size(r,2)
        if r(i) <= ( (R0./100).*R)
            AxVelTPL_radius(i) = (gradp.*( r(i).^2 )./(4*mu0)) + A0;
            
        elseif r(i) > ( (R0./100).*R) && r(i) <= ( (RInf./100).*R)
            AxVelTPL_radius(i) = -( ( -gradp./(2*K) ).^(1/n) )./alpha.*( r(i).^alpha );
            AxVelTPL_radius(i) = AxVelTPL_radius(i) + A1;
            
        elseif r(i) > ( (RInf./100).*R)
            AxVelTPL_radius(i) = -gradp./(4*muInf).*( R.^2 - r(i).^2 );
            
        end
    end
    
end