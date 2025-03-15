%              Initial stress calculation in linear Maxwell material
%
% Zhiren, Dec. 2024
% =========================================================================
%   This is my implementation of Sawyer's strategy to estimate initial
%   stress integral in a linear Maxwell material during bubble growth stage
%   of LIC experiment.
% =========================================================================

function S_peak = get_S0_SLS(Lam,Re,De,Ca)

%% Parameter definition
% Lam = Amplification factor of bubble
% Ca = Cauchy number
% Re = Reynolds number
% De = Deborah number

params = struct;
params.Lam = Lam;
params.Ca = Ca;
params.Re = Re;
params.De = De;

%% Iterate for initial velocity
R0 = 1/Lam;
V0_guess = sqrt( (2/3) * (R0^(-3) - 1) ); % This uses Rayleigh's solution for void collapse

Rmax_err = @(X) abs(findRmax(X,params) - 1);

V0_opt = fminsearch(Rmax_err,V0_guess);

%fm_options = optimset('Display','iter','PlotFcns',@optimplotfval);
%[V0_opt, fval, exitflag, output] = fminsearch(Rmax_err, V0_guess, fm_options); % Get access to iteration info

%% Run simulation again to get stress integral at end:

% Set up other parameters
Rb0 = 1/(params.Lam); % Initial radius
pb0 = 1;
SI0 = 0; % Stress integral

V0 = 1.1*V0_opt; % Use this structure so that we can modify and debug

tspan_star = 2.0; % Run long enough to reach end of growth
options = odeset('RelTol', 1E-7);

% Run simulation
X0 = [Rb0,V0,pb0,SI0]';
[t_sol,X_sol] = ode23tb(@(t,X) KM_bubble(t,X,params), [0,tspan_star], X0, options);

R_sol = X_sol(:,1);
[Rmax,indx] = max(R_sol); 

% DEBUG: Plot results
debug_plot(t_sol,X_sol);

tol = 1E-3;

if abs(Rmax - 1) > tol
    S_peak = nan;
else
    S_peak = X_sol(indx,4);

    % Deduct the NH part
    S_NH = -(5 - Lam^(-4) -4/Lam)/(2*Ca);
    S_peak = S_peak - S_NH;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                SUBROUTINES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% (A) Find Rmax for given initial velocity
function Rmax = findRmax(Vb0,params)

% Set up other parameters
Rb0 = 1/(params.Lam); % Initial radius
pb0 = 1;
SI0 = 0; % Stress integral

tspan_star = 2.0; % Run long enough to reach end of growth
options = odeset('RelTol', 1E-7);

% Run simulation
X0 = [Rb0,Vb0,pb0,SI0]';
[t_sol,X_sol] = ode23tb(@(t,X) KM_bubble(t,X,params), [0,tspan_star], X0, options);

% Get output
R_sol = X_sol(:,1);
Rmax = max(R_sol);

end

%% (B) Watered-down Keller-Miksis simulation 

function dxdt = KM_bubble(t,X,params)

Ca = params.Ca;
Re = params.Re;
De = params.De;
Lmax = params.Lam;

% Parameters we won't touch for now. Define here directly:
kap = 1.4;
c_star = 1484/10.0751; 

% Current values:
Rb = X(1);      % Bubble radius
Vb = X(2);      % Bubble wall velocity
pb = X(3);      % Isobaric pressure
SI = X(4);      % Stress integral

Lb = Lmax*Rb; % Current stretch

pbdot = -3*kap*pb*Vb/Rb;
SIdot_NH = -(2/Ca)*(Vb/Rb)*(Lb^(-4) + 1/Lb);
SI_NH = -(5 - Lb^(-4) - 4/Lb)/(2*Ca);
SIdot = -((SI - SI_NH) + (4/Re)*Vb/Rb)/De + SIdot_NH;

Rdot = Vb;

KM_RHS = (1 + Vb/c_star) * (pb + SI - 1) + (Rb/c_star)*(pbdot + SIdot);
KM_LHS = (3/2)*(1 - Vb/(3*c_star))*(Vb^2);
KM_denom = (1 - Vb/c_star)*Rb;
Vdot = (KM_RHS - KM_LHS)/KM_denom;

dxdt = [Rdot; Vdot; pbdot; SIdot];

end

%% (C) Plot for debug
function [] = debug_plot(t_sol,X_sol)

    figure(999);
    
    nsub = 4;
    fs = 20;

    str = ["$R^*$","$|\dot{R}^*|$","$p_{\rm b}^*$","$|S^*|$"];
    rgb = [8 81 156; 165,15,21; 35,139,69; 84,39,143]/255;

    for ii = 1:nsub

        subplot(1,nsub,ii)
        hold on; box on;
        set(gca,'TickLabelInterpreter','Latex','FontSize',fs)

        plot(abs(X_sol(:,ii)),t_sol,'Color',rgb(ii,:));

        xl = xlabel(str(ii));
        set(xl,'Interpreter','Latex','FontSize',fs)

        if ii == 1
            yl = ylabel("$t^*$");
            set(yl,'Interpreter','Latex','FontSize',fs)

            xlim([0,1])
        else
            set(gca,'xscale','log')

            % if ii == 4 || ii == 2
            %     xlim([1E-2 1E2])
            % end
        end
    end

end