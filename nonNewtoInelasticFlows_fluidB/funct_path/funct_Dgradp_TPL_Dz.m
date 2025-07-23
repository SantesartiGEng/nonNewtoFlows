function Dgradp_TPL_Dz = funct_Dgradp_TPL_Dz(gradpTPL, L, nodes_axis)
    % this function computes the derivative of the gradp field with a II-order accurate finite difference centred scheme 
    
    
    Dgradp_TPL_Dz = 0*gradpTPL;
    h = L/(nodes_axis-1);
    
    %%% axis - finite centred diff. scheme - II order accuracy
    for i=2:nodes_axis-1
        Dgradp_TPL_Dz(1,i) = gradpTPL(1,i+1)-gradpTPL(1,i-1);
        Dgradp_TPL_Dz(1,i) = Dgradp_TPL_Dz(1,i)/(2*h);
    end
    
    %%% inlet - finite forward diff. scheme - II order accuracy
    Dgradp_TPL_Dz(1,1) = -3/2*gradpTPL(1,1) + 2*gradpTPL(1,2) - 1/2*gradpTPL(1,3);
    Dgradp_TPL_Dz(1,1) = Dgradp_TPL_Dz(1,1)/h;
    
    %%% outlet - finite backward diff. scheme - II order accuracy
    Dgradp_TPL_Dz(1,nodes_axis) = 3/2*gradpTPL(1,nodes_axis) - 2*gradpTPL(1,nodes_axis-1) + 1/2*gradpTPL(1,nodes_axis-2);
    Dgradp_TPL_Dz(1,nodes_axis) = Dgradp_TPL_Dz(1,nodes_axis)/h;
end