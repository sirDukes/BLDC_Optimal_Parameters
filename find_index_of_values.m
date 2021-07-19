clear; clc; load('clean_sim_results_full.mat');

clearvars -except clean_sim_results clean_theoretical_results H_list ...
    I_list nDesignPoints P_list S_list T_list id_list ... 
    R_list Jm_list nP_list Bv_list L_list

global R_list Jm_list nP_list Bv_list L_list Kv_list Ke_list;
R_list = [];
Jm_list = [];
nP_list = [];
Bv_list = [];
L_list = [];
Kv_list = [];
Ke_list = [];

%% Gather BLDC Data from motor reuslts
for i=1:length(clean_sim_results)
    gather_bldc_sim_results(clean_sim_results{i})
end

%% Plot SimulationDesign Space
% plot_simulated_space();

%% Iteration through best results
search_list  = id_list(H_list > 0.82);
for i=search_list
    [Vdc,Rs,Ls,Ms,Psi,Ke,Kv,...
        N_max,No,Bv, Jm,Jsys, nPolePairs] ...
        = set_simulation_parameters(clean_theoretical_results{i});
    
    disp(clean_theoretical_results{i});
%     disp(clean_theoretical_results{i}.t_mech/clean_theoretical_results{i}.t_em);
    
     ee_bldc_ideal_input_performance_curves_plot;

      plot_design_point_performance_curve(idx,...
        torques,powers,efficiencys,speeds,currents,...
        mPower,mSpeed,mCurrent,mEfficiency)  % Plot valid design points
end 

%% Clear Vars
clearvars -except clean_sim_results clean_theoretical_results H_list ...
    I_list nDesignPoints P_list S_list T_list id_list ...
    R_list Jm_list nP_list Bv_list L_list


%% ---------------~FUNCTIONS-------------------
function [Vdc,Rs,Ls,Ms,Psi,Ke,Kv,N_max,No,Bv,Jm,Jsys,nPolePairs]= set_simulation_parameters(results)
    Vdc         = results.V_max;
    Rs          = results.R_bldc_Ohm;
    Ls          = results.L_minus_M;        % (L-M)_max*0.8
    Ms          = 0;
    Jm          = results.Jm;
    Jsys        = 0;
    Bv          = results.B_v;                % Bv value is set to half 
    Psi         = results.Psi;
    Ke          = results.Ke;
    Kv          = results.Kv;
    N_max       = results.N_max;
    nPolePairs  = results.nPolePairs;
    No          = results.N0_load_speed;
end 

function plot_design_point_performance_curve(idx,...
    torques,powers,efficiencys,speeds,currents,...
    mPower,mSpeed,mCurrent,mEfficiency)

   % Reuse figure if it exists, else create new figure
    if ~exist('h1_ee_bldc_ideal_input_performance_curves', 'var') || ...
            ~isgraphics(h1_ee_bldc_ideal_input_performance_curves, 'figure')
        h1_ee_bldc_ideal_input_performance_curves = figure('Name', 'ee_bldc_ideal_input_performance_curves');
    end
    figure(h1_ee_bldc_ideal_input_performance_curves)
    clf(h1_ee_bldc_ideal_input_performance_curves)
    hf = figure(1);

    set(hf,'Color',[1 1 1])
    
    plot(torques.signals.values(idx,2),powers.signals.values(idx,2)/max(powers.signals.values(:,2)),torques.signals.values(idx,2),currents.signals.values(idx,2)/max(currents.signals.values(:,2)),torques.signals.values(idx,2),speeds.signals.values(idx,2)/max(speeds.signals.values(:,2)),torques.signals.values(idx,2),efficiencys.signals.values(idx,2))
    grid
    axis([0 max(torques.signals.values(idx)) 0 1])
    xlabel('Torque (Nm)')
    ylabel('Magnitude (pu)')
        legend(['Power : ',num2str(ceil(10*mPower/1000)/10),'kW Max'],['Current : ',num2str(ceil(10*mCurrent)/10),'A Max'],['Speed :',num2str(ceil(10*mSpeed)/10),'RPM Max'],...
            ['Efficiency : ',num2str(mEfficiency*100), '%']);
    % Remove temporary variables
    clear hf hp
end

function gather_bldc_sim_results(sim_results)
    global R_list Jm_list nP_list Bv_list Kv_list Ke_list L_list;
    R_list(end+1) = sim_results.R_bldc_Ohm;
    Jm_list(end+1) = sim_results.Jm;
    nP_list(end+1) = sim_results.nPolePairs;
    Bv_list(end+1) = sim_results.B_v;
    Kv_list(end+1) = sim_results.Kv;
    Ke_list(end+1) = sim_results.Ke;
    L_list(end+1) = sim_results.L_minus_M;
end

function plot_simulated_space()
    global H_list R_list Jm_list nP_list Bv_list Kv_list Ke_list L_list;
    
    tiledlayout(2,2);
    % Tile 1
    nexttile
    scatter3(nP_list, Ke_list, H_list)
    % Tile 2
    nexttile
    scatter3(R_list, L_list, H_list)
    % TIle 3
    nexttile
    scatter3(R_list, Ke_list, H_list)
    
end