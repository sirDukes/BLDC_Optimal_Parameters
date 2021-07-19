function simulate_design_space(results) 
    global Vdc Rs Ls Ms Psi Ke Kv N_max No Bv Jm Jsys nPolePairs;
    [Vdc,Rs,Ls,Ms,Psi,Ke,Kv,...
            N_max,No,Bv, Jm,Jsys, nPolePairs] ...
            = set_simulation_parameters(results);

       ee_bldc_ideal_input_performance_curves_plot;
       
       plot_design_point_performance_curve(idx,...
           torques,powers,efficiencys,speeds,currents,...
           mPower,mSpeed,mCurrent,mEfficiency)  % Plot valid design points
 
     %% Plot Performance Curves
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
end