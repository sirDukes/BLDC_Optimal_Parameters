load('clean_sim_results_xx.mat')
global GR_list Jm_list Bv_list...
    R_min R_max i_max i_min GR_max GR_min nPoles_max nPoles_min Jm_max Jm_min Bv_max Bv_min
R_max = 0; R_min = 0;i_max = 0; i_min = 0;
GR_max=0; GR_min=0; nPoles_max=0; nPoles_min =0;
Jm_max =0;Jm_min =0;Bv_max =0;Bv_min=0;
GR_list = [];
Jm_list = [];
Bv_list = [];

for i =1:length(clean_sim_results)
    calc_param_ranges(clean_sim_results{i});
    
    if (clean_sim_results{i}.GR > 0) && (clean_sim_results{i}.GR < 10)  
                    fprintf(...
"\nGR:\t%d \nVdc:\t%d \nR:\t %d \nL:\t %d \nΨm:\t %d\nKe:\t %d\nKv:\t %d\nNmax:\t %d\nNo:\t %d\nBv:\t %d\nJm:\t %d\nnPolePairs:\t %d\n"...
,clean_sim_results{i}.GR,Vdc,Rs,Ls,Psi,Ke,Kv,N_max,No,Bv, Jm, nPolePairs...
            );
        
        [Vdc,Rs,Ls,Ms,Psi,Ke,Kv,...
        N_max,No,Bv, Jm,Jsys, nPolePairs] ...
            = set_simulation_parameters(clean_sim_results{i});
    
        ee_bldc_ideal_input_performance_curves_plot
        
        plot_design_point_performance_curve(idx,...
        torques,powers,efficiencys,speeds,currents,...
        mPower,mSpeed,mCurrent,mEfficiency)  % Plot valid design points
    end 
end

disp(R_min)
disp(R_max)

function calc_param_ranges(clean_sim_results)
    global Bv_list Jm_list GR_list ...
        R_min R_max i_max i_min GR_max GR_min nPoles_max nPoles_min Jm_max Jm_min Bv_max Bv_min
    %% Resistance
    if (clean_sim_results.R_bldc_Ohm > R_max)
        R_min = R_max;
        R_max = clean_sim_results.R_bldc_Ohm;
    end
    %% Current
    if (clean_sim_results.I_max > i_max)
        i_min = i_max;
        i_max = clean_sim_results.I_max;
    end 
    %% GR
    GR_list(end+1) = calc_GR(clean_sim_results.N_nom);
    GR_min = min(GR_list); GR_max = max(GR_list);
    %% Jm
    Jm_list(end+1) = clean_sim_results.Jm;
    Jm_min = min(Jm_list); Jm_max = max(Jm_list);
    %% Bv
    Bv_list(end+1) = clean_sim_results.B_v;
    Bv_min = min(Bv_list); Bv_max = max(Bv_list);
end 

function GR = calc_GR(MOTOR_n_nom)
    GR = (0.66*MOTOR_n_nom)/ (12.43*9.55);
end

function [Vdc,Rs,Ls,Ms,Psi,Ke,Kv,N_max,No,Bv,Jm,Jsys,nPolePairs]= set_simulation_parameters(results)
    Vdc     = results.V_max;
    Rs      = results.R_bldc_Ohm;
    Ls      = results.L_minus_M;        % (L-M)_max*0.8
    Ms      = 0;
    Jm      = results.Jm;
    Jsys    = 0;
    Bv      = results.B_v;                % Bv value is set to half 
    Psi     = results.Psi;
    Ke      = results.Ke;
    Kv      = results.Kv;
    N_max    = results.N_max;
    nPolePairs = results.nPolePairs;
    No       = results.N0_load_speed;
end 

function plot_design_point_performance_curve(idx,...
    torques,powers,efficiencys,speeds,currents,...
    mPower,mSpeed,mCurrent,mEfficiency)

   % Reuse figure if it exists, else create new figure
    if ~exist('h1_ee_bldc_ideal_input_performance_curves', 'var') || ...
            ~isgraphics(h1_ee_bldc_ideal_input_performance_curves, 'figure')
        h1_ee_bldc_ideal_input_performance_curves = figure('Name', 'ee_bldc_ideal_input_performance_curves');
    end
    figure(h1_ee_bldc_ideal_input_performance_curves);
    clf(h1_ee_bldc_ideal_input_performance_curves);
    hf = figure(1);

    set(hf,'Color',[1 1 1]);
    
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