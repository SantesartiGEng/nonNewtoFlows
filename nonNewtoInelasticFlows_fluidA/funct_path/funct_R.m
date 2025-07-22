function R = funct_R(z,theta,Rin)
    % this function returns the radius R(z) along the conical nozzle 
    
    R = Rin-z*tan(theta*pi/180);
    
end