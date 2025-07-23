function R0TPL = funct_R0_axis(R, gradpTPL, tau0)

    R0TPL(:) = -1./(gradpTPL(:)); 
    R0TPL = R0TPL*2*tau0; % R0 in [mm]
    R0TPL(:) = 100.*R0TPL(:)./R(:); % percentage [%]
    for i=1:length(R0TPL)
        if R0TPL(i) > 100
            R0TPL(i) = 100;
        end
    end

end