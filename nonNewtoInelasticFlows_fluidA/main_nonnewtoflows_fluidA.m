%% INFO
% This is an open-source MATLAB code to analyse non-Newtonian inelastic fluids 
% flowing in slightly tapered axisymmetric pipes.
%
% The flow analysis implemented regards the reference shear-thinning "FLUID A" characterised by two constant viscosity plateau connected by a shear-thinning branch described in: 
% Santesarti, G., Marino, M., Viola, F., Verzicco, R., & Vairo, G. (2025). A Quasi-Analytical Solution for "Carreau-Yasuda-like" Shear-thinning Fluids Flowing in Slightly Tapered Pipes. arXiv preprint 2502.14991 at https://doi.org/10.48550/arXiv.2502.14991
%
% If you use this code, please cite the reference article mentioned above.
%
% MATLAB Code written by Gianluca Santesarti, University of Rome Tor Vergata, February 2024.

%% IMPORTANT NOTE
%  all measurement units follow the SI except for the length expressed in [mm]

%% SETUP OF FLOW PROBLEM DATA
clc;
clear all;

path = strcat(pwd,'\funct_path');
addpath(path)  

%%% discretization along the axis and radius
nodes_axis = 3001;
nodes_radius = 101;

%%% Geometrical data
Rin = 1.5; %inlet radius [mm]
Rout = 0.25; %outet radius [mm]
L = 20; % axial length [mm]
theta = atand((Rin-Rout)/L); % 3.5763 deg
thetarad = theta*pi/180; % taper angle [rad]

z = linspace(0, L, nodes_axis); % axis array

R = linspace(0, L, nodes_axis);
R(:) = funct_R(z(:),theta,Rin);

%%% Rheological data - SRB model
mu0 = 200 ; %zero-shear rate viscosity [Pa*s] 
lambda0 = 1; % zero-time constant [s]
gamma0 = 1/lambda0; %zero-shear rate [1/s]
tau0 = mu0*gamma0; % stess_0 [Pa]

n = 0.5; %shear-thinning index [-]
K =  mu0*(lambda0^(n-1)); %consistency index [Pas^n]

lambdaInf = 0.05; % infinity-time constant [s]
gammaInf = 1/lambdaInf; %infinity-shear rate [1/s]
muInf = mu0*(lambdaInf/lambda0)^(1-n); % [Pas]
tauInf = muInf*gammaInf; % stess_0 [Pa]

a = 2; % Yasuda exponent

alpha = (n+1)/n;
beta = (3*n+1)/n;

%%% Rheological data - TPL model
mu0_TPL = mu0; %[Pas]
lambda0_TPL = lambda0; % zero-time constant [s]
gamma0_TPL = gamma0; %[1/s]
tau0_TPL = tau0; %[Pa]

n_TPL = n;

lambdaInf_TPL = lambdaInf; % infinity-time constant [s]
gammaInf_TPL = gammaInf; %[1/s]
muInf_TPL = muInf; %[Pas]
tauInf_TPL = tauInf; %[Pa]

%%% Extrusion-bioprinting process data
vExtrusion = 20; %extrusion velocity [mm/s]

Q = pi*Rout^2*vExtrusion; %flow rate [mm^3/s]

%% NEWTONIAN mu0 ANALYSIS SECTION
%%%--- Pressure
pNewto_mu0(:) = R(:).^(-3) - Rin^(-3);
pNewto_mu0 = -pNewto_mu0*( 8*Q*mu0 )/( 3*pi*tan(thetarad) ); 
pNewto_mu0 = pNewto_mu0 - pNewto_mu0(1,size(pNewto_mu0,2));% pressure [Pa]

%%%--- Pressure Gradient
gradpNewto_mu0(:) =  1./( R(:).^4 );
gradpNewto_mu0 = -gradpNewto_mu0*(Q*8*mu0)/pi; % pressure gradient [Pa/mm]

%% POWER LAW ANALYSIS SECTION
%%%--- Pressure 
pPowLaw(:) =  R(:).^(-3.*n) - Rin^(-3*n);
pPowLaw = -pPowLaw*( (Q*(3*n+1)/(pi*n))^n )*(2*K)/(3*n*tan(thetarad)); 
pPowLaw = pPowLaw - pPowLaw(1,size(pPowLaw,2));% pressure [Pa]

%%%--- Pressure Gradient
gradpPowlaw(:) =  R(:).^( -(3.*n+1) ) ;
beta1 = ( ( Q*(3*n+1) )/(pi*n) )^n;
gradpPowlaw = -gradpPowlaw*(2*K)*beta1; % pressure gradient [Pa/mm]

%% NEWTONIAN muInf ANALYSIS SECTION
%%%--- Pressure
pNewto_muInf(:) = R(:).^(-3) - Rin^(-3);
pNewto_muInf = -pNewto_muInf*( 8*Q*muInf )/( 3*pi*tan(thetarad) ); 

pNewto_muInf = pNewto_muInf - pNewto_muInf(1,size(pNewto_muInf,2));% pressure [Pa]

%%%--- Pressure Gradient
gradpNewto_muInf = z;
gradpNewto_muInf(:) =  1./( R(:).^4 );
gradpNewto_muInf = -gradpNewto_muInf*(Q*8*muInf)/pi; % pressure gradient [Pa/mm]

%%  TPL dp/dz COMPUTATION SECTION
gradpTPL = funct_gradpTPL_fluidA(gradpNewto_mu0, gradpPowlaw, gradpNewto_muInf, Q, R, mu0, tau0, K, n, muInf, tauInf);
pTPL = funct_pTPL(gradpTPL, L, nodes_axis);

%%
%%%--- R0 radius percentage of power law region beginning
R0TPL = funct_R0_axis(R, gradpTPL, tau0);

%%%--- RInf radius percentage of power law region ending
RInfTPL = funct_RInf_axis(R, gradpTPL, tauInf);

%% PLOT R0 RInf
xmin = 0; xmax = L; % z [mm]
ymin = 0;

figure('OuterPosition', [500 200 600 500])


P1 = plot(z(:), R0TPL(:),'k -','LineWidth',1); %z [mm] and R0 [%]
hold on
P2 = plot(z(:), RInfTPL(:),'k -.','LineWidth',1); %z [mm] and R0 [%]

xlim([xmin xmax])
ylim([ymin 105])
ax = gca; ax.FontSize = 16;
xlabel('$z$ [mm] ','Interpreter','latex','FontSize',20);
ylabel('$R_x/R\,[\%]$ ','Interpreter','latex','FontSize',20);
legend([P1 P2],{'$R_0/R$ - QA (TPL)','$R_\infty/R$ - QA (TPL)'}, 'Location','southwest','Interpreter','latex','fontsize',12);

title('$R_0/R$ and $R_\infty/R$','Interpreter','latex')
grid on
hold off

%% VISCOSITY & SHEAR RATE WINDOW AT THE OUTLET

x_optSRB = [mu0 lambda0 lambdaInf n a];
x_optTPL = [mu0_TPL lambda0_TPL lambdaInf_TPL n_TPL];

%%% Shear rate window at the outlet
gamma_outlet_1 = -gradpTPL(1,end)*0.05*Rout/(2*mu0); %evaluated at r=0.05R_out
gamma_outlet_2 = -gradpTPL(1,end)*Rout/(2*muInf); %evaluated at r=R_out

%%% arrays
xgrid1_T = logspace(-2, 5, 1001);
xgrid1 = transpose(xgrid1_T);

muSRB(:,1) = visc_modSRB(x_optSRB,xgrid1(:,1));

muTPL = 0*xgrid1;
for i=1:size(muTPL,1)
    muTPL(i,1) = visc_modTPL(x_optTPL,xgrid1(i,1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%------------------FIGURE-VISCOSITY--SRB+TPL----------------%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

marker_viscTPL = [gamma0_TPL mu0_TPL; gammaInf_TPL muInf_TPL];

figure('OuterPosition', [500 200 600 500]);

L1 = loglog(xgrid1(:,1), muSRB(:,1),'--k','LineWidth',1);
hold on
L2 = loglog(xgrid1(:,1), muTPL(:,1),'-k','LineWidth',1);
loglog(marker_viscTPL(:,1), marker_viscTPL(:,2), 'o k');
L3 = xline(gamma_outlet_1,'-.','linewidth',1.5);
L4 = xline(gamma_outlet_2,'-.','linewidth',1.5);

xpoint = gamma0_TPL - 0;
ypoint = mu0_TPL + 40;
text(xpoint,ypoint,'$X_0$','VerticalAlignment','middle','HorizontalAlignment','center','rotation',-0,'Interpreter','latex','fontweight','bold','FontSize',16)
xpoint = gammaInf_TPL - 3;
ypoint = muInf_TPL - 8;
text(xpoint,ypoint,'$X_\infty$','VerticalAlignment','middle','HorizontalAlignment','center','rotation',-0,'Interpreter','latex','fontweight','bold','FontSize',16)

xpoint = gamma_outlet_1;
ypoint = 10;
text(xpoint,ypoint,'$0.05R_{out}$','VerticalAlignment','middle','HorizontalAlignment','center','rotation',90,'Interpreter','latex','fontweight','bold','FontSize',14)

xpoint = gamma_outlet_2;
ypoint = 10;
text(xpoint,ypoint,'$R_{out}$','VerticalAlignment','middle','HorizontalAlignment','center','rotation',90,'Interpreter','latex','fontweight','bold','FontSize',14)

xlim([1e-1 1e3])
ylim([15 400]) 
ax = gca; ax.FontSize = 16; 
xlabel('$\dot{\gamma}$ [s$^{-1}$]','interpreter','latex','fontsize',20);
ylabel('$\mu$ [Pa$\cdot$s]','interpreter','latex','fontsize',20);

lgd = legend([L2 L1 L3],{'TPL','SRB','Shear rate window (outlet)'}, 'Location','northeast');
lgd.FontSize = 12;

title('Viscosity \& Shear rate window (outlet)','Interpreter','latex')
hold off
grid on
box on

%% PLOT - AXIS - AXIAL PRESSURE
figure('OuterPosition', [500 200 600 500])
P1 = plot(z(:), pTPL(:).*1e-3,'k -','LineWidth',0.8); % z [mm] p [kPa]

ax = gca; ax.FontSize = 16; 
xlabel('$z$ [mm] ','Interpreter','latex','FontSize',20);
ylabel('$p$ [kPa] ','Interpreter','latex','FontSize',20);
legend([P1],{'QA (TPL)'},'Location','northeast','Interpreter','latex','FontSize',12);

title('Axis - Pressure','Interpreter','latex')
grid on
xmin = 0; xmax = L; % z [mm]
ymin = 0; ymax = 200; % p [kPa]
xlim([xmin xmax])
ylim([ymin ymax])

%% PLOT - AXIS - AXIAL VELOCITY
AxialVelTPL = funct_AxialVel_TPL_HSR(z, gradpTPL, R, R0TPL, RInfTPL, mu0, K, n, muInf);

figure('OuterPosition', [500 200 600 500])

P1 = plot(z(:), AxialVelTPL(:),'k -','LineWidth',1); % z [mm] vz [mm/s]

ax = gca; ax.FontSize = 16; 
xlabel('$z$ [mm] ','Interpreter','latex','FontSize',20);
ylabel('$v_z$ [mm/s] ','Interpreter','latex','FontSize',20);
legend([P1],{'QA (TPL)'},'Location','northwest','Interpreter','latex','FontSize',12);
xmin = 0; xmax = L; % z [mm]
xlim([xmin xmax])

title('Axis - Axial velocity','Interpreter','latex')
grid on

%% OUTLET - AXIAL VELOCITY

r_Outlet = linspace(0, R(1,end), nodes_radius);  % r [mm]
AxVelTPL_Outlet = funct_AxVelTPL_radius_HSR(gradpTPL(1,end), r_Outlet, R0TPL(1,end), RInfTPL(1,end), mu0, K, n, muInf); % [mm/s]

figure('OuterPosition', [500 200 600 500]) 

P1 = plot(r_Outlet(:)*1e3, AxVelTPL_Outlet(:),'k -','LineWidth',0.8); % r [\mum] vz [mm/s] 

ax = gca; ax.FontSize = 16;
xlabel('$r$ [$\mu$m]','Interpreter','latex','FontSize',20);
ylabel('$v_z$ [mm/s]','Interpreter','latex','FontSize',20);
legend([P1],{'QA (TPL)'},'Location','northeast','Interpreter','latex','FontSize',12);
xlim([0 R(1,end)*1e3]) % r [\mum]
ylim([0 40])

title('Outlet - Axial velocity','Interpreter','latex')
grid on

%%  OUTLET - RADIAL VELOCITY EXPRESSIONS

DgradpTPLDz = funct_Dgradp_TPL_Dz(gradpTPL, L, nodes_axis);
index_z_axis = size(z,2); %outlet section
[radvel_outlet, radiusvect_outlet] = funct_radvel_HSR(index_z_axis, nodes_radius, R, thetarad, gradpTPL, DgradpTPLDz, mu0, tau0, K, n, muInf, tauInf, alpha, beta);

figure('OuterPosition', [500 200 600 500]) 

P1 = plot(radiusvect_outlet(:)*1e3, radvel_outlet(:),'k -','LineWidth',0.8);  % r [\mum] v in [mm/s] 

ax = gca; ax.FontSize = 16;
xlabel('$r$ [$\mu$m]','interpreter','latex','fontsize',20);
ylabel('$v_r$ [mm/s]','interpreter','latex','fontsize',20);
xlim([0 R(1,end)*1e3])  % r in [\mum]
ylim([-1.0 0])
legend([P1],{'QA (TPL)'},'Location','northeast','Interpreter','latex','FontSize',12);

title('Outlet - Radial velocity','Interpreter','latex')
grid on
