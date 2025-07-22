function visc2 = visc_modPL(p,Gamma)
%%%- Power law model 
    visc2 = p(1).*(Gamma.^(p(2)-1));
% %%%- Carreau-Yasuda model 
%     visc = p(4)+ (p(1)-p(4))*(1+(p(2).*Gamma).^p(5)).^((p(3)-1)/p(5));

end