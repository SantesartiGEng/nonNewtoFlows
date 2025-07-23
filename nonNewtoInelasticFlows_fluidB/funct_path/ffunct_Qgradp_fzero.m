function gradpTPL = ffunct_Qgradp_fzero(mingradp, Q, R, mu0, K, n, tau0)
%%%--- LSR FLOW RATE
betaLSR = (1/(2*K))^(1/n);
betaLSR = betaLSR*(4*(tau0^2)*n)/(n+1);

funQLSR1 = @(x) 2.*(tau0.^4)./( mu0.*(x.^3) );

funQLSR2 = @(x) betaLSR.*( x.^( (1./n)-2 ) ).*( R.^((n+1)./n) - (2.*tau0./x).^((n+1)./n) );

funQLSR = @(x) pi*( funQLSR1(x) + funQLSR2(x) );

%%%--- HSR FLOW RATE
betaHSRa = (1/(2*K))^(1/n);
betaHSRa = (pi*n)*betaHSRa/( (n+1)*(3*n+1) );
    
funQHSR1a = @(x) (3.*n+1).*( R.^( (n+1)./n ) ) - 2.*n.*( (  (2.*tau0)./x  ).^( (n+1)./n ) );
funQHSR1b = @(x) -4.*(tau0.^2)./(x.^2);
funQHSR1c = @(x) funQHSR1a(x)*funQHSR1b(x);

betaHSRb = (n+1).*R.^( (3.*n+1)./n );
funQHSR1 = @(x) funQHSR1c(x) + betaHSRb;

funQHSR2 = @(x) betaHSRa.*( x.^(1./n) );

funQHSR = @(x) funQHSR1(x)*funQHSR2(x);

%%%--- TOTAL FLOW RATE
funQ = @(x) funQLSR(x) + funQHSR(x) - Q;

x0 = mingradp;
options = optimset('PlotFcns',{@optimplotx,@optimplotfval});
gradpTPL = fzero(funQ,x0);
gradpTPL = -gradpTPL; %gradp is always negative along the nozzle
 

end