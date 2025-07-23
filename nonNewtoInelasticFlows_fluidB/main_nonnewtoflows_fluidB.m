%% INFO
% This is an open-source MATLAB code to analyse non-Newtonian inelastic fluids 
% flowing in slightly tapered axisymmetric pipes.
%
% The flow analysis implemented regards the reference shear-thinning "FLUID B" characterised by one constant viscosity plateau and a shear-thinning branch described in: 
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

R = z;
R(:) = funct_R(z(:),theta,Rin);

%%% Rheological data - SRB model
mu0 = 200 ; %zero-shear rate viscosity [Pa*s] 
lambda0 = 1; % zero-time constant [s]
gamma0 = 1/lambda0; %zero-shear rate [1/s]
tau0 = mu0*gamma0; % stess_0 [Pa]

n = 0.2; %shear-thinning index  [-]
K =  mu0*(lambda0^(n-1)); %consistency index [Pas^n]
a = 2; % Yasuda exponent

alpha = (n+1)/n;
beta = (3*n+1)/n;

%%% Rheological data - TPL model
mu0_TPL = mu0; %[Pa*s]
lambda0_TPL = lambda0; % zero-time constant [s]
gamma0_TPL = gamma0; %[1/s]
tau0_TPL = mu0_TPL*gamma0_TPL;

n_TPL = n;

%%% Process data
vExtrusion = 20; %extrusion velocity [mm/s]
Q = pi*Rout^2*vExtrusion; %flow rate [mm^3/s]

%% NEWTONIAN mu0 ANALYSIS SECTION

%%%--- Pressure
pNewto_mu0(:) = R(:).^(-3) - Rin^(-3);
pNewto_mu0 = -pNewto_mu0*( 8*Q*mu0 )/( 3*pi*tan(thetarad) ); 
pNewto_mu0 = pNewto_mu0 - pNewto_mu0(1,size(pNewto_mu0,2)); % pressure [Pa]

%%%--- Pressure Gradient
gradpNewto_mu0 = z;
gradpNewto_mu0(:) =  1./( R(:).^4 );
gradpNewto_mu0 = -gradpNewto_mu0*(Q*8*mu0)/pi; % pressure gradient in [Pa/mm]


%% POWER LAW ANALYSIS SECTION

%%%--- Pressure 
pPowLaw(:) =  R(:).^(-3.*n) - Rin^(-3*n) ;
pPowLaw = -pPowLaw*( (Q*(3*n+1)/(pi*n))^n )*(2*K)/(3*n*tan(theta*pi/180)); % pressure [Pa]
pPowLaw = pPowLaw - pPowLaw(1,size(pPowLaw,2));

%%%--- Pressure Gradient
gradpPowlaw(:) =  R(:).^( -(3.*n+1) ) ;
beta1 = ( ( Q*(3*n+1) )/(pi*n) )^n;
gradpPowlaw = -gradpPowlaw*(2*K)*beta1; % pressure gradient in [Pa/mm]

%%  TPL dp/dz COMPUTATION SECTION
gradpTPL = funct_gradpTPL_fluidB(gradpNewto_mu0, gradpPowlaw, Q, R, mu0, K, n, tau0);
pTPL = funct_pTPL(gradpTPL, L, nodes_axis);

%%
%%%--- R0 radius percentage of power law region beginning
R0TPL = funct_R0_axis(R, gradpTPL, tau0);

%% PLOT R0
figure('OuterPosition', [500 200 600 500])

xmin = 0; xmax = L; % z [mm]
ymin = 0;

P1 = plot(z(:), R0TPL(:),'k -','LineWidth',1); %z [mm] and R0 [%]

ax = gca; ax.FontSize = 16;
xlim([xmin xmax])
ylim([ymin 100])
xlabel('$z$ [mm] ','Interpreter','latex','FontSize',20);
ylabel('$R_0/R\,[\%]$ ','Interpreter','latex','FontSize',20);
legend([P1],{'$R_0/R$ - QA (TPL)'}, 'Location','southwest','Interpreter','latex','fontsize',12);

title('$R_0/R$','Interpreter','latex')
grid on

%% VISCOSITY & SHEAR RATE WINDOW AT THE OUTLET

x_optSRB = [mu0 lambda0 n a];
x_optTPL = [mu0_TPL lambda0_TPL n_TPL];

%%% Shear rate window at the outlet
gamma_outlet_1 = -gradpTPL(1,end)*0.05*Rout/(2*mu0); %evaluated at r=0.05R_out
gamma_outlet_2 = ( -gradpTPL(1,end)/(2*K)*Rout )^(1/n); %evaluated at r=R_out

%%% arrays
xgrid1_T = logspace(-2, 5, 1001);
xgrid1 = transpose(xgrid1_T);

muSRB(:,1) = visc_modSRB_subcase1(x_optSRB,xgrid1(:,1));

muTPL = 0*xgrid1;
for i=1:size(muTPL,1)
    muTPL(i,1) = visc_modTPL_subcase1(x_optTPL,xgrid1(i,1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%------------------FIGURE-VISCOSITY--SRB+TPL----------------%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
marker_viscTPL = [gamma0_TPL mu0_TPL];

figure('OuterPosition', [500 200 600 500])

L1 = loglog(xgrid1(:,1), muSRB(:,1),'--k','LineWidth',1);
hold on
L2 = loglog(xgrid1(:,1), muTPL(:,1),'-k','LineWidth',1);
loglog(marker_viscTPL(:,1), marker_viscTPL(:,2), 'o k');
L3 = xline(gamma_outlet_1,'-.','linewidth',1.2);
L4 = xline(gamma_outlet_2,'-.','linewidth',1.2);

xpoint = gamma0_TPL+2;
ypoint = mu0_TPL + 60;
text(xpoint,ypoint,'$X_0$','VerticalAlignment','middle','HorizontalAlignment','center','rotation',-0,'Interpreter','latex','fontweight','bold','FontSize',16)

xpoint = gamma_outlet_1;
ypoint = 1;
text(xpoint,ypoint,'$0.05R_{out}$','VerticalAlignment','middle','HorizontalAlignment','center','rotation',90,'Interpreter','latex','fontweight','bold','FontSize',14)

xpoint = gamma_outlet_2;
ypoint = 1;
text(xpoint,ypoint,'$R_{out}$','VerticalAlignment','middle','HorizontalAlignment','center','rotation',90,'Interpreter','latex','fontweight','bold','FontSize',14)

xlim([1e-2 1e4])
ylim([1 800])
xlabel('$\dot{\gamma}$ [s$^{-1}$]','interpreter','latex','fontsize',20);
ylabel('$\mu$ [Pa$\cdot$s]','interpreter','latex','fontsize',20);

lgd = legend([L2 L1 L3],{'TPL','SRB','Shear rate window (outlet)'}, 'Location','northeast');
lgd.FontSize = 12;

title('Viscosity \& Shear rate window (outlet)','Interpreter','latex')
ax = gca; ax.FontSize = 16; 
hold off
grid on
box on

%% PLOT - AXIS - AXIAL PRESSURE
figure('OuterPosition', [500 200 600 500])
P1 = plot(z(:), pTPL(:).*1e-3,'k -','LineWidth',0.8); % z [mm] p [kPa]

ax = gca; ax.FontSize = 16; 
xmin = 0; xmax = L; % z [mm]
xlim([xmin xmax])
xlabel('$z$ [mm] ','Interpreter','latex','FontSize',20);
ylabel('$p$ [kPa] ','Interpreter','latex','FontSize',20);
lgd = legend([P1],{'QA (TPL)'},'Location','northeast','Interpreter','latex','FontSize',12);
lgd.FontSize = 12;

title('Axis - Pressure','Interpreter','latex')
grid on
hold off

%% PLOT - AXIS - AXIAL VELOCITY
AxialTPL = funct_AxialVel_TPL(z, gradpTPL, R, R0TPL, mu0, K, n);

figure('OuterPosition', [500 200 600 500])
P1 = plot(z(:), AxialTPL(:),'k -','LineWidth',0.8); % z [mm] vz [mm/s]

ax = gca; ax.FontSize = 16; 
xmin = 0; xmax = L; % z [mm]
xlim([xmin xmax])
xlabel('$z$ [mm] ','Interpreter','latex','FontSize',20);
ylabel('$v_z$ [mm/s] ','Interpreter','latex','FontSize',20);
lgd = legend([P1],{'QA (TPL)'},'Location','northwest','Interpreter','latex','FontSize',12);
lgd.FontSize = 12;

title('Axis - Axial velocity','Interpreter','latex')
grid on

%% ANSYS FLUENT COMPARISON SECTION - AXIAL VELOCITY - OUTLET SECTION
%%%---Velocity at the OUTLET
r_Outlet = linspace(0, R(1,end), nodes_radius);  % r [mm]
AxVelTPL_Outlet = funct_AxVelTPL_radius(gradpTPL(1,end), r_Outlet, R(1,end), R0TPL(1,end), mu0, K, n);

figure('OuterPosition', [500 200 600 500])
P1 = plot(r_Outlet(:)*1e3, AxVelTPL_Outlet(:),'k -','LineWidth',0.8); % r in [mum] vz in [mm/s]

ax = gca; ax.FontSize = 16;
xmin = 0; xmax = r_Outlet(1,end)*1e3; % r [\mum]
ymin = 0; ymax = 40; %[mm/s]
xlim([xmin xmax])
ylim([ymin ymax])
xlabel('$r$ [$\mu$m]','Interpreter','latex','FontSize',20);
ylabel('$v_z$ [mm/s]','Interpreter','latex','FontSize',20);
lgd = legend([P1],{'QA (TPL)'},'Location','northeast','Interpreter','latex','FontSize',12);
lgd.FontSize = 12;

title('Outlet - Axial velocity','Interpreter','latex')
grid on

%% ANSYS FLUENT COMPARISON SECTION - RADIAL VELOCITY - OUTLET SECTION

DgradpTPLDz = funct_Dgradp_TPL_Dz(gradpTPL, L, nodes_axis);
index_z_axis = size(z,2); %outlet section
[radvel_outlet, radiusvect_outlet] = funct_radvel_MSR(index_z_axis, nodes_radius, R, thetarad, gradpTPL, DgradpTPLDz, mu0, tau0, K, n, alpha);

figure('OuterPosition', [500 200 600 500]) 
P1 = plot(radiusvect_outlet(:)*1e3, radvel_outlet(:),'k -','LineWidth',0.8); % r [\mum] vr [mm/s] 

ax = gca; ax.FontSize = 16;
xmin = 0; xmax = r_Outlet(1,end)*1e3; % r [\mum]
xlim([xmin xmax])
xlabel('$r$ [$\mu$m]','interpreter','latex','fontsize',20);
ylabel('$v_r$ [mm/s]','interpreter','latex','fontsize',20);
lgd = legend([P1],{'QA (TPL)'},'Location','northeast','Interpreter','latex','FontSize',12);
lgd.FontSize = 12;

title('Outlet - Radial velocity','Interpreter','latex')
grid on
