%% Abstract 
% The analysis and control of powertrain systems of electric vehicle, which 
% is an important type of new energy vehicle, have been the focus of extensive 
% research, but determining the motor modeling parameters remains a problem. 
% A method of parameter determination for brushless DC motor modeling based 
% on vehicle power performance was developed in this study. The power and torque 
% of the driving motor of an electric vehicle were obtained by using the 
% dynamic equation of the electric vehicle to satisfy the requirements of 
% power performance. The ranges of the back electromotive force coefficient and 
% the winding inductance were derived from the voltage and dynamic equations of
% brushless DC motor, which were deduced from the expected power and torque 
% of the motor. The modeling parameters were then determined on the basis 
% of the influence of power source voltage, back electromotive force coefficient,
% winding inductance, and winding resistance on vehicle power performance. 
% A hardware-in-loop simulation of vehicle power performance was performed 
% to verify the effectiveness of the proposed method. Results indicate that 
% the maximum vehicle velocity is 172 km/h, and the acceleration time of 100 km/h 
% is 13 s, which reveal that the motor modeling parameters obtained with the 
% method satisfy relevant requirements.

%% Reference Paper: Determination of modeling parameters for a brushless DC motor that satisfies the power performance of an electric vehicle
%   Author:          Zhiyong Zhang 
%   Link to paper:   https://journals.sagepub.com/doi/full/10.1177/0020294019842607

%% Load Mareks Model 
% AA_Mareks_Model;
% clear;
%% Vehicle Performance Parameters
global EV_a_front EV_mass EV_rmcc EV_Crr EV_Cd EV_GR_max EV_GR EV_Rw_m;
global g rho ;

% GR_bldc = 33;

g                   = 9.81;
rho                 = 1.205;
EV_a_front          = 0.63;
EV_mass             = 100;
EV_rmcc             = 1;
EV_Crr              = 0.012;
EV_Cd               = 0.6;
EV_GR_max           = 31.2;      
EV_Rw_m             = inch2mtr(26)/2;

% Speeds
EV_v_max_mps        = kph2mps(25);              % m/s
EV_v_cruise_mps     = kph2mps(25);              % m/s
EV_v_acc_mps        = kph2mps(25);              % m/s
EV_v_nom_mps        = kph2mps(20);              % m/s
EV_v_max_rpm        = mps2rpm(EV_v_max_mps);    % RPM
EV_v_nom_rpm        = mps2rpm(EV_v_nom_mps);    % RPM
EV_v_max_rps        = rpm2rps(EV_v_max_rpm);    % rad/s

% acceleration times
EV_acc_time         = 0.1; % in m.s^-2
%% Other Vehicle Performance
EV_slope_rad        = deg2rad(8);

%is the base speed ratio, which is set to 4–5 for passenger cars . Gao Y and Ehsani M. A torque and speed coupling hybrid drivetrain-architecture, control, and simulation. 
%IEEE T Power Electr 2006; 21(3): 741–748.
EV_beta             = 2; 

%% Battery Parameters
BATTERY_V_max = 40.7; % [V] unit test

%% Motor Predetermined Parameters
% Fixed Parameters
MECH_itta = 1;
MOTOR_itta = 1;       % efficiency 
% Sweeped Parameters
% MOTOR_nPoles = 10;
% MOTOR_R = 0.025;
% MOTOR_I_max = 35;
% MOTOR_Jm = 1e-6;
%Link with Determine_optimal_BLDC_Params_full previous script
BATTERY_V_max = V_max_bldc; % [V]
MOTOR_nPoles = nPoles;
MOTOR_R = R_bldc;
MOTOR_I_max = I_max_bldc;
MOTOR_Jm = Jm_bldc;
EV_GR               = GR_bldc;   

%% Vehicle Power calculatios
EV_P_v = (calc_F_rr(0) + calc_F_drag(EV_v_max_mps))*EV_v_max_mps/(MECH_itta*MOTOR_itta);
EV_P_i = (calc_F_rr(EV_slope_rad) + calc_F_drag(EV_v_cruise_mps) + calc_F_ramp(EV_slope_rad))*EV_v_cruise_mps/(MECH_itta*MOTOR_itta);
EV_P_j = (calc_F_rr(0) + calc_F_drag(EV_v_acc_mps) + calc_F_acc(EV_acc_time))*EV_v_acc_mps/(MECH_itta*MOTOR_itta);

%% Determination of characteriatic speed
MOTOR_n_nom = (9.55*EV_GR*EV_v_nom_mps)/EV_Rw_m;
MOTOR_n_max = (9.55*EV_GR*EV_v_max_mps)/EV_Rw_m;
MOTOR_n_base = MOTOR_n_max/EV_beta;

%% Power and Torque Matching
MOTOR_Te = EV_Rw_m*(g*EV_mass*EV_Crr + EV_rmcc*EV_mass*EV_acc_time + 0.5*rho*EV_a_front*EV_Cd*EV_v_max_mps^2)/EV_GR;
MOTOR_P_max = max([EV_P_v, EV_P_i, EV_P_j]);                % OUPUT

MOTOR_T_j = EV_P_j*60/(EV_GR*2*pi*EV_v_max_rpm);        % OUPUT
MOTOR_T_i = EV_P_i*60/(EV_GR*2*pi*EV_v_max_rpm);        % OUPUT
MOTOR_Te_max = max([MOTOR_T_j,MOTOR_T_i]);                   % OUPUT

% MOTOR_F_drive  = MOTOR_Te*EV_GR/EV_Rw_m
%% Determination fo back-EMF coefficient
MOTOR_ke_min = MOTOR_Te_max/(2*MOTOR_I_max);
MOTOR_ke_max = 9.55*BATTERY_V_max/(2*MOTOR_n_max);
MOTOR_Kv = calc_Kv(MOTOR_ke_max);
MOTOR_Ke = (MOTOR_ke_min+MOTOR_ke_max)/2;

%% Determination of Max Pm flux linkage
MOTOR_Psi = calc_PM_max_flux_link(MOTOR_Ke, MOTOR_nPoles/2);

MOTOR_P_nom = (calc_F_rr(0) + calc_F_drag(EV_v_nom_mps))*EV_v_nom_mps/MOTOR_itta;
MOTOR_Te_nom = (3/2)*(MOTOR_nPoles/2)*MOTOR_Psi*(MOTOR_I_max*0.8);

%% Inductance Determination
MOTOR_n0 = 15*BATTERY_V_max/(pi*MOTOR_ke_min);
MOTOR_L_minus_M_max = calc_L_M(MOTOR_ke_min,MOTOR_n0, MOTOR_n_max, MOTOR_R, MOTOR_Te, MOTOR_nPoles);
MOTOR_L_minus_M = MOTOR_L_minus_M_max;

%% Friction Drag Coefficient
MOTOR_Bv = calc_Bv(MOTOR_Te_nom, EV_v_max_rpm*EV_GR);
% MOTOR_Bv = Bv_bldc;

%% calculate time constnts 
[MOTOR_t_mech, MOTOR_t_em] = calc_BLDC_time_constants(MOTOR_R, MOTOR_L_minus_M, MOTOR_Ke, MOTOR_Jm, MOTOR_Bv);
%% Compile MOTOR_params;
MOTOR_params = struct(               ...
    'GR',               EV_GR,              ...
    'P_max',            MOTOR_P_max,        ...    
    'P_nom',            MOTOR_P_nom,        ...
    'V_max',            BATTERY_V_max,      ...
    'I_max',            MOTOR_I_max,        ...
    'nPolePairs',       MOTOR_nPoles/2,     ...
    'R_bldc_Ohm',       MOTOR_R,            ...
    'L_minus_M',        MOTOR_L_minus_M_max,...
    'K_emf_min_rpVs',   MOTOR_ke_min,       ...
    'K_emf_max_rpVs',   MOTOR_ke_max,       ...
    'Psi',              MOTOR_Psi,          ...
    'Ke',               MOTOR_Ke,           ...
    'Kv',               MOTOR_Kv,           ...
    'B_v',              MOTOR_Bv,           ...
    'N0_load_speed',    MOTOR_n0,           ...
    'N_max',            MOTOR_n_max,        ...
    'N_nom',            MOTOR_n_nom,        ...
    'N_base',           MOTOR_n_base,       ...
    'Te_max',           MOTOR_Te_max,        ...
    'Te_max_speed',     MOTOR_Te,           ...
    'Te_nom',           MOTOR_Te_nom,        ...
    'Jm',               MOTOR_Jm,            ...
    't_mech',           MOTOR_t_mech,        ...
    't_em',             MOTOR_t_em            ...
);
MOTOR_params;

%% clear workspace
clearvars -except MOTOR_params  idx Bv_bldc V_max_bldc I_max_bldc Psi_bldc R_bldc nPoles GR_bldc Jm_bldc ... 
    RUN_NEW nDesignPoints Bv_bldc_lst GR_bldc_lst Psi_bldc_lst Jm_bldc_lst R_bldc_lst I_bldc_lst V_dc_lst nPoles_lst theoretical_results clean_theoretical_resultls

%% ------------------- Functions -------------------------------------
function Rw_m  = inch2mtr(Rw_in)
    Rw_m = Rw_in*0.0254;
end

function v_rps = rpm2rps(v_rpm)
    v_rps = v_rpm*0.10471975513824;
end

function v_rpm = mps2rpm(v_mps)
    global EV_Rw_m
    v_rpm  = 60*v_mps/(2*pi*EV_Rw_m);
end

function v_mps = kph2mps(v_kph)
    v_mps = v_kph*0.277778;
end

function F_rr = calc_F_rr(slope)
    global EV_mass EV_Crr g
    F_rr=EV_mass*g*EV_Crr*cos(slope);
end

function F_drag = calc_F_drag(v_mps)
    global rho EV_a_front EV_Cd
    F_drag=0.5*rho*EV_a_front*EV_Cd*(v_mps^2);
end

function F_a = calc_F_acc(acc_time)
    global EV_mass EV_rmcc
    F_a=EV_rmcc*EV_mass*acc_time;
end

function F_ramp = calc_F_ramp(slope)
    global EV_mass g
    F_ramp = EV_mass*g*sin(slope);
end

function Kv = calc_Kv(Ke)
    Kv = 60/(2*pi*Ke);
end

function L_minus_M_max = calc_L_M(MOTOR_ke_min,MOTOR_n0, MOTOR_n_max, MOTOR_R, MOTOR_Te, MOTOR_nPoles)
    L_minus_M_max = (4*pi*(MOTOR_ke_min^2)*(MOTOR_n0 - MOTOR_n_max)-60*MOTOR_R*MOTOR_Te)/(3*MOTOR_nPoles*MOTOR_n_max*MOTOR_Te);
end

function Bv = calc_Bv(Te_Nm, v_rpm)
    Bv = Te_Nm / rpm2rps(v_rpm) ;
end

function [t_mech, t_em] = calc_BLDC_time_constants(R, L, Ke, Jm, Bv)
    R = R*2;
    L=L*2;
    t_mech = (R*Jm + L*Bv) / (R*Bv + Ke*Ke);
    t_em = (L*Jm) / (R*Jm + L*Bv);
end

function Psi_m = calc_PM_max_flux_link(Ke, nPolePairs)
    Psi_m = Ke/(sqrt(3)*2*pi*1000*nPolePairs/60);
end