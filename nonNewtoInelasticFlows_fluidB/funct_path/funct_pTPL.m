function pTPL = funct_pTPL(gradp, L, nodes)
    % this function computes the axial pressure field by solving the linear system associated to the gradp, with a II-order accurate finite difference centred scheme 
    
    gradp_T = transpose(gradp);
    gradp_T = gradp_T*1e-3; % [kPa/mm] scaling of matrix coefficients to avoid singularities
    pTPL_T = gradp_T*0;
    
    h = (L/(nodes-1));  % [mm]
    A = zeros(nodes,nodes);
    
    for i=2:nodes-1
        A(i,i-1) = -1/(2*h);
        A(i,i+1) = 1/(2*h);
    end
    
    %inlet foward diff. scheme II order accurate
    A(1,1) = -3/(2*h);
    A(1,2) = 2/h;
    A(1,3) = -1/(2*h);
    
    %outlet bc null pressure
    A(nodes,nodes) = 1;
    gradp_T(nodes) = 0;
   
    pTPL_T = linsolve(A,gradp_T);
    pTPL = transpose(pTPL_T);
    pTPL = pTPL*1e3; % [Pa]
    
end