function RInfTPL = funct_RInf_axis(R, gradpTPL, tauInf)

    RInfTPL(:) = -1./(gradpTPL(:)); 
    RInfTPL = RInfTPL*2*tauInf; % RInf in [mm]
    RInfTPL(:) = 100.*RInfTPL(:)./R(:); % percentage [%]
    for i=1:length(RInfTPL)
        if RInfTPL(i) > 100
            RInfTPL(i) = 100;
        end
    end

end