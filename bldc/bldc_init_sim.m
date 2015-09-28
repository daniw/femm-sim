%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Magnetic simulation of a BLDC motor using FEMM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants for simulation setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Resolution for simulation iterations
res_cogg = 1;
res_rot = 1;
res_el  = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate needed variables for simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation ranges
cogg_range = res_cogg:res_cogg:(360 * 3 / stat_nof_slots);
rot_range = res_rot:res_rot:(360 * 3 / stat_nof_slots);
i_range =  res_el:res_el:360;
disp(strcat(num2str(length(cogg_range) + length(rot_range) * length(i_range)), ' Simulations'));
disp(strcat('  ', num2str(length(cogg_range)),                                 ' for Cogging Torque'));
disp(strcat('  ', num2str(length(rot_range) * length(i_range)),                ' for Torque'));
