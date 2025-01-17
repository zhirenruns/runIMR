%                   Simulation of Laser-Induced Cavitation
%               
% Zhiren Zhu (zhiren@umich.edu)
% Dec. 2024
% =========================================================================
% Usage:
%
%   This is intended to be a clean update to the IMR simulation code
%   written by Jon Estrada leading up to his 2018 paper. (Ref [1])
%   The code structure is partially inspired by the code written by Jin
%   Yang during his post-doc at Wisconsin.
%
%   The present script is a self-contained run file for parameter sweep.
%   It is modified from 
%   Brute-force sweep will be implemented for now. 
%   Next step is Nelder-Mead search to match experiment.
%
%   Only do NHKV model for now.
%
% =========================================================================
% References:
%   [1] Estrada, Barajas, Henann, Johnsen & Franck (2018) JMPS 
%           (https://doi.org/10.1016/j.jmps.2017.12.006)
%   [2] Prosperetti, Crum & Commander (1988) JASA
%           (https://doi.org/10.1121/1.396145)
%   [3] Preston (2004) Ph.D. Thesis @ Caltech
%           (https://www.proquest.com/dissertations-theses/modeling-heat-mass-transfer-bubbly-cavitating/docview/305200526/se-2)
%   [4] Ascher & Greif, A First Course in Numerical Methods (SIAM)
%
% =========================================================================

clc; close all; clearvars;
addpath('graphics');

% Parameter Range
%   Adjust loop structure if additional parameters are studied in the future.

% Loop parameters
Rmax_Array = (1:0.2:4)*(1E-4);      % Max. bubble radius (m)
Lmax_Array = (4:0.2:10);           % Max. stretch of bubble (=Rmax/Req)
Tinf_Array = linspace(273.15,323.15,3); % Far-field temperature

G_Array = 0;        % Elastic modulus (Pa)
mu_Array = 0;       % Viscous modulus (Pa*s)

% Get size of arrays
nR = length(Rmax_Array);
nL = length(Lmax_Array);
nT = length(Tinf_Array);
nG = length(G_Array);
nmu = length(mu_Array);

% Constant parameters:
p_inf = 101325;             % Far-field pressure (Pa)   
rho = 998.2;                % Mass density (kg/m^3)
c_wave = 1484;              % Longitudinal wave speed (m/s)
gam = 7.0E-2;               % Surface tension (N/m)

% Store constant:
in_para = struct;

in_para.p_inf = p_inf;  
in_para.rho = rho;
in_para.c_wave = c_wave;
in_para.gam = gam;

% Make data storage structure
imr_soln(nR,nL,nT,nG,nmu) = struct;

% Initialize figure:
% figure(100);
% cmap = viridis(nT);

for aa = 1:nR
    Rmax_here = Rmax_Array(aa);
    in_para.Rmax = Rmax_here;

    % Adjust time span according to Rmax:
    tspan_here = Rmax_here/3;
    in_para.tspan = tspan_here;

    for bb = 1:nL
        Lmax_here = Lmax_Array(bb);
        in_para.Lmax = Lmax_here;

        %plot_indx = nL*(aa-1) + bb;
        %subplot(nR,nL,plot_indx);

        disp([aa,bb]); % Just to show progress
        tic;

        for cc = 1:nT
            Tinf_here = Tinf_Array(cc);
            in_para.T_inf = Tinf_here;
            
            for ii = 1:nG
                G_here = G_Array(ii);
                in_para.G_elast = G_here;

                for jj = 1:nmu
                    mu_here = mu_Array(jj);
                    in_para.mu_visc = mu_here;

                    % Run simulation for this case:
                    sol_here = IMR_Driver(in_para);

                    t_here = sol_here(:,1);
                    R_here = sol_here(:,2);
                    pb_here = sol_here(:,3);
                    
                    % Plot:
                    % cplot = cmap(cc,:);
                    % plot_history(t_here,R_here,cplot)

                    % Find collapse time here:
                    nmin_list = find(islocalmin(R_here) == 1);

                    if isempty(nmin_list)
                        tc_here = nan;
                    else
                        tc_here = t_here(nmin_list(1));
                    end

                    % Store:
                    imr_soln(aa,bb,cc,ii,jj).t = t_here;
                    imr_soln(aa,bb,cc,ii,jj).R = R_here;
                    imr_soln(aa,bb,cc,ii,jj).tc = tc_here;

                    imr_soln(aa,bb,cc,ii,jj).Rmax = Rmax_here;
                    imr_soln(aa,bb,cc,ii,jj).Lmax = Lmax_here;
                    imr_soln(aa,bb,cc,ii,jj).Tinf = Tinf_here;
                    imr_soln(aa,bb,cc,ii,jj).G = G_here;
                    imr_soln(aa,bb,cc,ii,jj).mu = mu_here;

                    imr_soln(aa,bb,cc,ii,jj).pb = pb_here;
                end
            end
        end
        toc;
    end
end

%%
% Save Result
savename = ['check-pressure_incompressible_no-surf-tension_',datestr(now,'yyyymmdd')];
save([savename,'.mat'],'imr_soln');

%%
% =========================================================================
%                              SUBROUTINES
% =========================================================================

%% (I) Driver

function sol_out = IMR_Driver(para_in)

% Unpack parameters
% (A) LIC Experiment 
Rmax = para_in.Rmax;     % Max. bubble radius (m)
Lmax = para_in.Lmax;     % Max. stretch of bubble (=Rmax/Req)
tspan = para_in.tspan;   % Temporal duration studied (s)

% (B) Far-field/atmosphere
T_inf = para_in.T_inf;             % Far-field temperature (K)
p_inf = para_in.p_inf;             % Far-field pressure (Pa) 

% (C) Surrounding material
G_elast = para_in.G_elast;          % Elastic modulus (Pa)
mu_visc = para_in.mu_visc;          % Viscous modulus (Pa*s)
rho = para_in.rho;                  % Mass density (kg/m^3)
c_wave = para_in.c_wave;            % Longitudinal wave speed (m/s)
gam = para_in.gam;                  % Surface tension (N/m)

% (D) Bubble content        (Define here for now, rather than input)
kap = 1.4;                  % Heat capacity ratio
D_bin = 2.42E-5;            % Binary diffusion coefficient (m^2/s)
Ru = 8.3144598;             % Universal gas constant (J/mol*K)
mwv = 18.01528E-3;          % Molecular weight of vapor (kg/mol)
mwg = 28.966E-3;            % Molecular weight of non-vapor gas (kg/mol) [~ Air]
Rv = Ru/mwv;                % Gas constant for vapor
Rg = Ru/mwg;                % Gas constant for non-vapor gas

% (D-i) Empirical constant for thermal conductivity 
%   Per Ref [2], K = A*T + B works well for air at 200 to 3000 Kelvin
A_thm = 5.28E-5;            % (W/m*K^2)
B_thm = 1.17E-2;            % (W/m*K)

% (D-ii) Empirical constant for saturated temperature of vapor
%   Per Appendix A3 of Ref [3], pvsat(T) = p_ref*exp(-T_ref/T) 
p_ref = 1.17E11;            % Reference pressure (Pa)
T_ref = 5.2E3;              % Reference temperature (K)

%% Parameters for computation
NY = 500;                   % # of mesh points for bubble content calculation
solver_RelTolX = 1E-7;      % Relative tolerance for ODE solver

%% Intermediate & non-dimensionalization parameters
% Req = Rmax/Lmax;                      % Equilibrium radius (m) [Not used for calculation]
vc = sqrt(p_inf/rho);                   % Characteristic velocity (m/s)
tc = Rmax/vc;                           % Characteristic time scale (s)
K_inf = A_thm*T_inf + B_thm;            % Thermal conductivity (W/m*K)
pv_inf = p_ref*exp(-T_ref/T_inf);       % Vapor pressure at equilibrium/far-field (Pa)

c_star = c_wave/vc;                     % Non-dim of wave speed
We = p_inf*Rmax/(2*gam);                % Weber number
Ca = p_inf/G_elast;                     % Cauchy number
Re = sqrt(p_inf*rho)*Rmax/(mu_visc);    % Reynolds number
Fom = D_bin/(vc*Rmax);                  % Mass Fourier number
chi = K_inf*T_inf/(p_inf*vc*Rmax);      % Lockhart-Martinelli number
A_star = A_thm*T_inf/K_inf;             % Non-dim of A (= 1 - B-star)
B_star = B_thm/K_inf;                   % Non-dim of B
pv_star = pv_inf/p_inf;                 % Non-dim of far-field vapor pressure
tspan_star = tspan/tc;                  % Non-dim of time span studied

%% Package all parameters
% For easier update down the road, let's use struct rather than array

params = struct;

params.NY = NY;
params.Ca = Ca;
params.Re = Re;
params.We = We;
params.Fom = Fom;
params.chi = chi;
params.kap = kap;
params.A_star = A_star;
params.B_star = B_star;
params.c_star = c_star;
params.pv_star = pv_star;
params.Rv = Rv;
params.Rg = Rg;
params.Lmax = Lmax;

%% Initial conditions

Rb0 = 1;                                            % Bubble radius
Vb0 = 0;                                            % Bubble wall velocity
Theta0 = zeros(1,NY);                               % Thermal field - See Eqn.(52) of Ref [1] 
pb0 = pv_star + (1 + Lmax/We - pv_star)/(Lmax^3);   % Isobaric pressure
kv0_mag = 1/(1 + (Rv/Rg)*(pb0/pv_star - 1));        % Magnitude of vapor mass fraction - See Eqn.(56) of Ref [1] 
kv0 = kv0_mag*ones(1,NY);                           % Vapor mass fraction field
SI0 = -(5 - Lmax^(-4) - 4/Lmax)/(2*Ca);             % Stress integral

X0 = [Rb0,Vb0,pb0,SI0,Theta0,kv0]';

%% March Keller-Miksis forward

options = odeset('RelTol', solver_RelTolX);
[t_sol,X_sol] = ode23tb(@(t,X) KM_bubble(t,X,params), [0,tspan_star], X0, options);

%% Post-process

% Extract (dimensionless) solution:
R_sol = X_sol(:,1);
V_sol = X_sol(:,2);
pb_sol = X_sol(:,3);
SI_sol = X_sol(:,4);

indx_T = 4;
indx_k = indx_T + NY;

Theta_sol = X_sol(:,(indx_T + 1):(indx_T + NY));
kv_sol = X_sol(:,(indx_k + 1):(indx_k + NY));

% Convert to dimensional values:
t_dim = t_sol * tc;
R_dim = R_sol * Rmax;
pb_dim = pb_sol * p_inf;
T_field = T_inf*(-B_star + sqrt(1 + (2*A_star)*Theta_sol) )./A_star;  % Dimensional temperature field

% Output:
sol_out = [t_dim,R_dim,pb_dim];

end

%% (II) Keller-Miksis steps
% Call this function to march governing equation forward

function dxdt = KM_bubble(t,X,params)

% READ INFO ================================================
% Unpack params:
NY = params.NY;
Ca = params.Ca;
Re = params.Re;
We = params.We;
Fom = params.Fom;
chi = params.chi;
kap = params.kap;
A_star = params.A_star;
B_star = params.B_star;
c_star = params.c_star;
pv_star = params.pv_star;
Rv = params.Rv;
Rg = params.Rg;
Lmax = params.Lmax;

% Current values:
Rb = X(1);      % Bubble radius
Vb = X(2);      % Bubble wall velocity
pb = X(3);      % Isobaric pressure
% SI = X(4);    % Stress integral. Evaluate directly from kinematics.

indx_T = 4;
indx_k = indx_T + NY;
Theta = X((indx_T + 1):(indx_T + NY));  % Thermal field ('theta')
kv = X((indx_k + 1):(indx_k + NY));     % Vapor mass fraction ('k')

% BUBBLE CONTENT ================================================

% Mesh for bubble content
dy = 1/(NY - 1);                % Uniform spacing
Ymesh = (dy*((1:NY)-1))';       % Dimensionless mesh from center to wall of bubble

% Dirichlet BC for vapor mass fraction
kv(end) = 1/(1 + (Rv/Rg)*(pb/pv_star - 1));

% Calculate mixture field
Tb = (-B_star + sqrt(1 + (2*A_star)*Theta) )./A_star;        % Temperature field (normalized by T_inf) - Integrate Eqn. (52) of Ref [1] to get quadratic equation
K_star = A_star.*Tb + B_star;      % Mixture thermal conductivity
Rmix = kv.*Rv + (1 - kv).*Rg;      % Mixture gas-constant field

% Find spatial derivatives
%       Neumann BC at center + Central difference for interior points + Backward difference at bubble wall
%       See: Ref [4] - Section 16.4 for backward difference formula. The formula here is BDF2.
%       For Laplacian in axisymmetric problem, nabla^2(f) = (1/r^2)*d(r^2*df/dr)/dr.
%           Expand: nabla^2(f) = d^2(f)/dr^2 + (2/r)*df/dr. Hence, the DD shown below.
%       For center of bubble, use L'Hopital to obtain Laplacian: 
%           nabla^2(f[0]) = 3 d^2(f[0])/dr^2 ~= (df[h/2]/dr - 0)/(h/2)  
%       Laplacian at bubble wall is not actually used, since we impose BC.
%           Equation shown uses (f'(1-h/2) - f'(1-h))/(h/2) for d^2(f[1])/dr^2.

BC1 = 0;
CD1= (Theta(3:end) - Theta(1:(end-2)))/(2*dy);             
BD1 = (3*Theta(end) - 4*Theta(end - 1) + Theta(end-2))/(2*dy);

DTheta = [BC1; CD1; BD1];

BC2 = (6*(Theta(2) - Theta(1)))/(dy^2);
CD2 = diff(diff(Theta))/(dy^2) + (2./Ymesh(2:(end-1))).*DTheta(2:(end-1));
BD2 = 0;        % (2*Theta(end) - 5*Theta(end-1) + 4*Theta(end-2) - Theta(end - 3))/(dy^2) + (2/Ymesh(end))*DTheta(end); % Not used
DDTheta = [BC2; CD2; BD2];

BC3 = 0;
CD3= (kv(3:end) - kv(1:(end-2)))/(2*dy);             
BD3 = (3*kv(end) - 4*kv(end - 1) + kv(end-2))/(2*dy);
Dkv = [BC3; CD3; BD3];

BC4 = (6*(kv(2) - kv(1)))/(dy^2);
CD4 = diff(diff(kv))/(dy^2) + (2./Ymesh(2:(end-1))).*Dkv(2:(end-1));
BD4 = 0;        % (2*kv(end) - 5*kv(end-1) + 4*kv(end-2) - kv(end - 3))/(dy^2) + (2/Ymesh(end))*Dkv(end); % Not used
DDkv = [BC4; CD4; BD4];

% Time rate of isobaric bubble pressure
%       See Eqn. (55) of Ref [1]
pb_part1 = -kap*pb*Vb + (kap-1)*(chi/Rb)*DTheta(end);
pb_part2 = kap*pb*Rv*Fom*Dkv(end)/(Rmix(end)*(1-kv(end))*Rb); 
pbdot = (3/Rb)*(pb_part1 + pb_part2);

% Mixture velocity field
%       See Eqn. (54) of Ref [1]
Vm_part1 = ((kap - 1)*(chi/Rb)*DTheta - Rb*pbdot*Ymesh/3)/(kap*pb);
Vm_part2 = Fom*(Rv - Rg)*Dkv./(Rb*Rmix);
Vmix = Vm_part1 + Vm_part2;

% Evolution of thermal and vapor mass fraction
%       See Eqn. (53) of Ref [1]
rhs_theta = pbdot + chi*DDTheta/(Rb^2) + (kap/(kap-1))*(pb./(K_star.*Tb))*(Fom/(Rb^2)).*((Rv-Rg)./Rmix).*Dkv.*DTheta;
denom_theta = (kap/(kap-1))*pb./(K_star.*Tb);
subtr_theta = (Vmix - Vb*Ymesh).*DTheta/Rb;
Theta_prime = rhs_theta./denom_theta - subtr_theta;
Theta_prime(end) = 0;   % Equilibrium imposed at boundary
   
rhomy = (1./Tb).*DTheta./sqrt(1 + 2*A_star*Theta) + (Rv-Rg)./Rmix.*Dkv; % = -(1/rho_m)*d(rho_m)/dy. Note: from Eqn. (20) of Ref [1], take spatial derivative to get this result.
rhs_kv = Fom/(Rb^2) * (DDkv - Dkv .* rhomy);
subtr_kv = (Vmix - Vb*Ymesh).*(Dkv)/Rb;
kv_prime = rhs_kv - subtr_kv;
kv_prime(end) = -(Rv/Rg)*(pbdot/pv_star)*(kv(end))^2;   % Equilibrium imposed at boundary (!= 0)


% SURROUNDING MATERIAL ================================================
Lb = Lmax*Rb;      % Current stretch of bubble
S_Ca = -(5 - Lb^(-4) - 4/Lb)/(2*Ca);
S_Re = -4*Vb/(Re*Rb);
SI = S_Ca + S_Re;

Sdot_Ca = (2/Ca)*(Vb/Rb)*(Lb^(-4) + 1/Lb); 
Sdot_Re1 = (4/Re)*(Vb/Rb)^2; % Acceleration-independent part
Sdot_Re2 = -(4/Re)/Rb; % Acceleration-dependent part, divided by acceleration
SIdot = Sdot_Ca + Sdot_Re1;

% INCREMENTS ================================================
Rdot = Vb;

KM_RHS = (1 + Vb/c_star) * (pb - 1/(We*Rb) + SI - 1) + (Rb/c_star)*(pbdot + SIdot + Vb/(We*Rb^2));
KM_LHS = (3/2)*(1 - Vb/(3*c_star))*(Vb^2);
KM_denom = (1 - Vb/c_star)*Rb - Sdot_Re2*(Rb/c_star);
Vdot = (KM_RHS - KM_LHS)/KM_denom;

SIdot = SIdot + Sdot_Re2*Vdot;  % Not actually used, but update properly

dxdt = [Rdot; Vdot; pbdot; SIdot; Theta_prime; kv_prime];

end

%% (III) Plot history of bubble dynamics
function [] = plot_history(t,R,rgb)

hold on; box on;

% Parameters:
to_micro = 1E6;
lw = 1;         % Line width
fs = 12;        % General font size
fs_ax = 16;     % Axis label font size

% Assign default color
if nargin<3
    rgb = [8 81 156]/255; % Blue color
end

plot(t*to_micro,R*to_micro,'LineWidth',lw,'Color',rgb)

set(gca,'TickLabelInterpreter','Latex','FontSize',fs)
xl = xlabel('$t ({\rm \mu s})$');
yl = ylabel('$R ({\rm \mu m})$');
set(xl,'Interpreter','Latex','FontSize',fs_ax);
set(yl,'Interpreter','Latex','FontSize',fs_ax)

pbaspect([ 2 1 1 ])

end

%% (IV) Plot contour of bubble content field
function [] = plot_field(t,X,NY,cmap)

hold on; box on;

% Parameters:
to_micro = 1E6;
lw = 1;         % Line width
fs = 12;        % General font size
fs_ax = 16;     % Axis label font size

% Assign default color map
if nargin<4
    cmap = flip(brewermap(256,'Spectral'),1); 
end

dY = 1/NY;
Y = dY:dY:1;

contourf(t,Y,X', 0:0.1:10, 'EdgeColor','none')

set(gca,'TickLabelInterpreter','Latex','FontSize',fs)
xl = xlabel('$t ({\rm \mu s})$');
yl = ylabel('$y = r/R$');
set(xl,'Interpreter','Latex','FontSize',fs_ax);
set(yl,'Interpreter','Latex','FontSize',fs_ax)
pbaspect([ 2 1 1 ])

colormap(cmap);
clim([0.5,1.5])
cbar = colorbar;
set(cbar,'TickLabelInterpreter','Latex','FontSize',fs,'Location','EastOutside')
cbar_ttl = get(cbar,'Title');
set(cbar_ttl ,'String',"$T/T_{\infty}$",'Interpreter',"Latex",'FontSize',fs);

end


