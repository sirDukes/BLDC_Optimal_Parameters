% Code to plot simulation results from ee_bldc_ideal_input_performance_curves
%% Plot Description:
% The plot below shows the performance curves of BLDC using ideal input.

% Generate new simulation results if they don't exist or if they need to be updated
if ~exist('simlog_ee_bldc_ideal_input_performance_curves', 'var') || ...
        simlogNeedsUpdate(simlog_ee_bldc_ideal_input_performance_curves, 'ee_bldc_ideal_input_performance_curves') 
    sim('ee_bldc_ideal_input_performance_curves')
    % Model StopFcn callback adds a timestamp to the Simscape simulation data log
end

idx = torques.signals.values(:,2) > 0;

mPower = max(powers.signals.values(:,2));
mCurrent = max(currents.signals.values(:,2));
mSpeed = max(speeds.signals.values(:,2));
mEfficiency = max(efficiencys.signals.values(:,2));
mTorque = max(torques.signals.values(:,2));
minTorque = min(torques.signals.values(idx,2));
minCurrent = min(currents.signals.values(idx,2));