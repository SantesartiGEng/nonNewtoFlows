function R = funct_R(z,theta,Rin)
    % this function returns the ray R(z) along the conical nozzle with
    % input the axial position z 
    %R0 = 4.35/2;
    %RL = 0.41/2;
    %alpha = 3.523;
    
    R = Rin-z*tan(theta*pi/180);
    
end