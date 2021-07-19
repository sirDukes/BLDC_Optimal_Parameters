tic
clear;

%% Setup
global clean_theoretical_results;
load('clean_theoretical_results.mat')
global clean_sim_results;
clean_sim_results = cell(length(clean_theoretical_results),1);
global P_list I_list H_list S_list T_list id_list;
P_list = [];
I_list = [];
H_list = [];
S_list = [];
T_list = [];
id_list = [];

%% Gather data from thoeretical model - Needed for Parallel Computing

%% Start of for-loop
for i=1:length(clean_theoretical_results) 
%     disp('Checking BLDC:');
%     disp(clean_theoretical_results{i});
    [Vdc,Rs,Ls,Ms,Psi,Ke,Kv,...
        N_max,No,Bv, Jm,Jsys, nPolePairs] ...
        = set_simulation_parameters(clean_theoretical_results{i});

       ee_bldc_ideal_input_performance_curves_plot;
        
       bool = filter_simulation_results(...
           mPower, mEfficiency,mTorque...
           );% Use results generated from bldc_performance_curve.m
% % % 
%        plot_design_point_performance_curve(idx,...
%            torques,powers,efficiencys,speeds,currents,...
%            mPower,mSpeed,mCurrent,mEfficiency);  % Plot valid design points
%        
       if bool
            fprintf(...
"\nVdc:\t%d \nR:\t %d \nL:\t %d \nΨm:\t %d\nKe:\t %d\nKv:\t %d\nNmax:\t %d\nNo:\t %d\nBv:\t %d\nJm:\t %d\nnPolePairs:\t %d\n"...
,Vdc,Rs,Ls,Psi,Ke,Kv,N_max,No,Bv, Jm, nPolePairs...
            );

            clean_sim_results{i} = clean_theoretical_results{i};
            
            gather_data_from_simulation(i, ...
                 mPower, mSpeed, mTorque,mCurrent, mEfficiency...
                 ); % Gather lists of KPI from valid design points
             
            plot_design_point_performance_curve(idx,...
                torques,powers,efficiencys,speeds,currents,...
                mPower,mSpeed,mCurrent,mEfficiency)  % Plot valid design points

       end
       fprintf("Completed:\t %d %c\n",int16(100*i/length(clean_theoretical_results)),"%");

end %% End of for-loop

%% Remove empty entries
clean_sim_results = clean_sim_results(~cellfun('isempty', clean_sim_results));

%% Record time 
toc

%% Save Clean sim results 
save('clean_sim_results_yy.mat');
%% Clear vars
% clearvars -except params clean_theoretical_results R_list L_list Ke_list Jm_list nPoles_list id_list tm_list tem_list

%% Functions
%% Gather Sim data - incomplete
function gather_data_from_simulation(i, mPower, mSpeed,mTorque,mCurrent, mEfficiency)
global P_list I_list H_list S_list T_list id_list;    
id_list(end+1) = i; 
P_list(end+1) = mPower;
I_list(end+1) = mCurrent;
H_list(end+1) = mEfficiency;
S_list(end+1) = mSpeed;
T_list(end+1) = mTorque;
end
%% Scatter Pllot - Redundant 
function scatter_plots()
    global R_list L_list R_list L_list Ke_list Jm_list nPoles_list id_list tem_list tm_list;
    
    t = tiledlayout(7,1);
    nexttile
    scatter(id_list, L_list);
    nexttile
    scatter(id_list, nPoles_list);
    nexttile
    scatter(id_list, 9.55*Ke_list);
    nexttile
    scatter(id_list, R_list);
    nexttile
    scatter(id_list, Jm_list);
    nexttile
    scatter(id_list, 1e3*tem_list);
    nexttile
    scatter(id_list, 1e3*tm_list);
end
%% Set Simulation Parameters 
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
%% Filter Simulation results
function bool = filter_simulation_results(mPower, mEfficiency,mTorque)
%     mSpeed_cond         = mSpeed > 5000;
%     minTorque_cond      = minTorque < 0.3;
    mTorque_cond        = mTorque > 2 && mTorque < 10;
%     mCurrent_cond       = mCurrent < 80;
    mEfficiency_cond    = mEfficiency >0.75;
%     mTemp_cond          = any(temperatures.max < 35); 
% %     mPower_cond         = mPower > 600;
     
    if mEfficiency_cond && mTorque_cond
        bool =1;
    else
        bool =0;
    end
    
%     disp(mTemp_cond); disp(mSpeed_cond),disp(minTorque_cond),disp(mEfficiency_cond);
    
end
%% Plot Performance Curves
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