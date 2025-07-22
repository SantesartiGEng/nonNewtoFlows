function gradpTPL = ffunct_Qgradp_HSR_fzero(Abs_gradpNewto_mu0, Q, R, mu0, tau0, K, n, muInf, tauInf)
   alpha = (n+1)/n;
   beta = (3*n+1)/n;

   %%%--- HSR FLOW RATE
   
   %%--- for r <= R_0  
   funA0_HSR_1 = @(x) (1/alpha).*((x./(2*K))^(1/n))*( (2*tauInf./x).^alpha - (2*tau0./x).^alpha );
   funA0_HSR = @(x) funA0_HSR_1(x) + (x./4)*( ( 4*tau0^2./(mu0.*(x.^2)) ) + R.^2/muInf - ( 4*tauInf^2./(muInf.*(x.^2)) )) ;

   funQLSR_HSR = @(x) pi*( ((4*tau0^2)./(x.^2)).*funA0_HSR(x) - ((2*tau0^4)./( mu0.*(x.^3) ) ) );

   %%--- for R_0 < r <= R_Inf
   funA1_HSR_1 = @(x) (1/alpha).*((x./(2*K))^(1/n))*( (2*tauInf./x).^alpha );
   funA1_HSR = @(x) (x./(4*muInf))*(R.^2 - (2*tauInf./x).^2) + funA1_HSR_1(x);
   
   funQMSR_HSR_1 = @(x) 2/(alpha*beta)*((x./(2*K))^(1/n)).*((2*tauInf./x).^beta - (2*tau0./x).^beta);
   funQMSR_HSR = @(x) pi*( funA1_HSR(x)*( (2*tauInf./x).^2 - (2*tau0./x).^2 )  -   funQMSR_HSR_1(x) );
   
   %%--- for r > R_Inf
   funQHSR_HSR = @(x) (pi.*x./(8*muInf)).*( ( R.^2 - (2*tauInf./x).^2  ).^2 );

   %%%--- TOTAL FLOW RATE
   funQ = @(x) funQLSR_HSR(x) + funQMSR_HSR(x) + funQHSR_HSR(x) - Q;

   x0 = Abs_gradpNewto_mu0;
   options = optimset('PlotFcns',{@optimplotx,@optimplotfval});
   gradpTPL = fzero(funQ,x0);
   gradpTPL = -gradpTPL; %gradp is always negative along the nozzle

end