function gradpTPL = ffunct_Qgradp_MSR_fzero(mingradp, Q, R, mu0, tau0, K, n)
   alpha = (n+1)/n;
   beta = (3*n+1)/n;
   %%%--- MSR FLOW RATE
   
   %%--- for r <= R_0
   funA0_MSR = @(x) (1/alpha).*((x./(2*K))^(1/n)).*( R.^alpha - ( (2*tau0)./x)^alpha ) + (tau0^2)./(mu0.*x);

   funQLSR_MSR = @(x) pi*( ((4*tau0^2)./(x.^2)).*funA0_MSR(x) - ((2*tau0^4)./( mu0.*(x.^3) ) ) );

   %%--- for r > R_0
   funQMSR_MSR_1 = @(x) 0.5*(R.^alpha).*( R.^2 - (4*tau0^2)./(x.^2) ) ;
   funQMSR_MSR_2 = @(x) (1/beta)*( R.^beta - (2*tau0./x).^beta ) ;
   funQMSR_MSR = @(x) (2*pi/alpha)*( (x./(2*K)).^(1/n) )*( funQMSR_MSR_1(x) - funQMSR_MSR_2(x)  );

   %%%--- TOTAL FLOW RATE
   funQ = @(x) funQLSR_MSR(x) + funQMSR_MSR(x) - Q;

   x0 = mingradp;
   options = optimset('PlotFcns',{@optimplotx,@optimplotfval});
   gradpTPL = fzero(funQ,x0);
   gradpTPL = -gradpTPL; %gradp is always negative along the nozzle

end