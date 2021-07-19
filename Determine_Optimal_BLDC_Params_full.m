%% Run AA_Mareks_Model.m 
%% Src to Determine_Optimal_BLDC_Params_Script.m
clear;

%% -------------------- ITERATION SETUP ------------------------------------------
% V_BATT
V_bldc_min = 4.2*10; V_bldc_max = 4.2*10; % 6S -11S 
V_dc_lst = V_bldc_min:4.2:V_bldc_max;                                 % From 6S - 10S

%%  Reduction ratio list
GR_min = 6.85; GR_max = 70; 
GR_bldc_lst = GR_min : 10 :GR_max;
%% nPole
nPoles_min =    14; nPoles_max = 18;
nPoles_lst      =   nPoles_min : 2 : nPoles_max;   % to start Parfor Loop nPoles' must evaluate to a row vector of consecutive increasing integers. For more information,
%% R_BLDC 
R_bldc_min =  40;                         % Number of points per
R_bldc_max   = 300;                          % miliOhm
a = R_bldc_max/20;
R_bldc_lst   = linspace(1, R_bldc_max, 20)/1000;
%% I_BLDC
I_bldc_min      =   1;
I_bldc_max      =   100;
I_bldc_lst      =   linspace(I_bldc_min,I_bldc_max,15);
%% Jm 
Jm_bldc_lst = [ 1e2, 1e3, 1e4, 1e5]*5e-9; % 1e-7*kg*m^2
%% Bv
Bv_bldc_lst = [1e2, 1e3, 1e4 ,1e5]*5e-7; % 1e-7*kg*m^2

nDesignPoints = length(Bv_bldc_lst)*length(GR_bldc_lst)*length(Jm_bldc_lst)*length(V_dc_lst)*(length(nPoles_lst))*length(R_bldc_lst)*length(I_bldc_lst);

idx = 0;
theoretical_results = cell(nDesignPoints,1);
%% Start Recording time
tic
%% ------------------------ EXHAUSTIVE SEARCH ---------------------------
        
for V_max_bldc = V_dc_lst   %%%%%%%%%% EDIT!
    for I_max_bldc = I_bldc_lst
        for nPoles = nPoles_lst
            for R_bldc = R_bldc_lst
                for Jm_bldc = Jm_bldc_lst
                    for GR_bldc = GR_bldc_lst
%                         for Bv_bldc = Bv_bldc_lst
                            idx = idx+1;
%                             fprintf("Iteration #: %d\nV_dc: \t%d\t\tI_bldc: %d\t\tnPoles: %d\t\tR_bldc: %d\n",idx, V_max_bldc, I_max_bldc, nPoles, R_bldc);
                            fprintf("Completed:\t%d", int16(100*idx/nDesignPoints));
                            disp("%");
                            Optimal_BLDC_Params_
                            bool = sparse_parse_matrix(MOTOR_params);
                            if bool
                                theoretical_results{end+1} = MOTOR_params;
                            end
%                         end
                    end
                end
            end
        end
   end
end

%% Clear Vas ro free memory
clearvars -except... 
    Bv_bldc_lst GR_bldc_lst nDesignPoints R_bldc_lst I_bldc_lst V_dc_lst nPoles_lst Jm_bldc_lst... 
    theoretical_results clean_theoretical_resultls 
%% Clean up results
clean_theoretical_results = theoretical_results(~cellfun('isempty', theoretical_results));
%% Save Clean Theoretical BLDc results
save('clean_theoretical_results.mat');

%% call plot_BLDC_Design_space.mat
% plot_BLDC_Design_space

%% Stop recording time
toc

%% ----------------------------- FUNCTIONS ----------------------------------------------------------------------------
%% 1. Sparse Parse Materix Creation Function
function bool = sparse_parse_matrix(params)
    % This function creates return a boolean based on if the current design point being investigated  

    cond_1 = (params.L_minus_M > 0);                                            % Physical feasibility: ( L - M ) > 0                                             % Only Even idx stored
%     cond_4 = (params.R_bldc_Ohm >= 10e-3)  || (params.R_bldc_Ohm <= 250e-3);    % 10mΩ  <  R_bldc   <  150mΩ
%      cond_4 = (params.R_bldc_Ohm >= 250e-3)  || (params.R_bldc_Ohm <= 500e-3);    % 10mΩ  <  R_bldc   <  150mΩ
    cond_7 = (params.K_emf_min_rpVs <= params.K_emf_max_rpVs);
    cond_8 = (params.N0_load_speed > params.N_max*1.1);
    cond_9 = (params.t_mech >0.00001) && (params.t_em >0.00001) ...
        && (params.t_mech > params.t_em) && (params.t_mech < 7*params.t_em);
    
    bool =  cond_1 ... 
            && cond_7 ...
            && cond_8 ...
            && cond_9 ...
            ;
end