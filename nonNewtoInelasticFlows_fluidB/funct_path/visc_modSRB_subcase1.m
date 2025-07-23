function visc = visc_modSRB_subcase1(p,Gamma)

% % %%%- SRB model
%     visc = p(1).*(  ( 1+ (Gamma.*p(3)).^p(5)  )./( 1+ (Gamma.*p(2)).^p(5) )  ).^( ( 1-p(4) )./p(5) );
% %%%- SRB model - subcase 1
    visc = p(1)./( ( 1+ (Gamma.*p(2)).^p(4) ).^( ( 1-p(3) )./p(4) ) );
end