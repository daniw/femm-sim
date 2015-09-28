%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Magnetic simulation of a BLDC motor using FEMM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants for defining motor dimensions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rot_nof_poles = 10;
rot_nof_mag_per_pole = 2;
stat_nof_slots = 15;
mot_thickness = 10;
stat_nof_wdg = 50;

rot_mag_thick = 2;
rot_mag_width = 9;
rot_inner_rad = 30.5;
rot_outer_rad = 35;

stat_out_rad = 30;
stat_base_rad = 20;
stat_inner_rad = 16;

stat_slot_width = 4;
stat_head_width = 10;
%                   Width   Radius
stat_slot_shape =  [0.      0.      ;
                    0.      0.88    ;
                    0.1     0.9     ;
                    0.95    0.9     ;
                    1.      1.      ];

mat_air     = 'Air';
mat_stat    = 'Vanadium Permedur';
mat_rot     = 'Vanadium Permedur';
mat_magnet  = 'NdFeB 52 MGOe';
%mat_magnet  = 'NdFeB 40 MGOe';
mat_wind    = '1mm';

current = 5;

% Groups for stator and rotor
group_stator = 1;
group_rotor = 2;
group_air = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate needed variables for motor dimensions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stator edge
stat_head_edge = sqrt((stat_out_rad)^2 - (stat_head_width/2)^2);
% stator base edge
stat_base_edge = stat_base_rad * cos(asin(stat_slot_width / 2 / stat_base_rad));
% scale stator slot
stat_slot_shape_r(:,1) = (stat_slot_shape(:,1) .* (stat_head_width - stat_slot_width) .+ stat_slot_width) ./ 2;
stat_slot_shape_r(:,2) = stat_slot_shape(:,2) .* (stat_head_edge - stat_base_edge) .+ stat_base_edge;
stat_slot_shape_l(:,1) = -stat_slot_shape_r(:,1);
stat_slot_shape_l(:,2) = stat_slot_shape_r(:,2);
% stator head angle
stat_head_angle = asin(stat_head_width/2 / stat_out_rad) * 180 / pi;
% stator base angle
stat_slot_edge_angle = atan(stat_slot_shape_r(1,1) / stat_slot_shape_r(1,2)) *180 / pi;
stat_slot_base_angle = 360 / stat_nof_slots / 2 - stat_slot_edge_angle;
% Stator winding radius
stat_wind_rad = sqrt(stat_slot_shape_r(length(stat_slot_shape) - 1, 1) ^ 2 + stat_slot_shape_r(length(stat_slot_shape) - 1, 2) ^ 2);
% Stator winding angle
stat_wind_angle = atan(stat_slot_shape_r(length(stat_slot_shape) - 1, 1) / stat_slot_shape_r(length(stat_slot_shape) - 1, 2)) * 180 / pi;
% Stator winding contour
stat_wind_cont = [stat_wind_rad * sin(stat_wind_angle * pi / 180) stat_wind_rad * cos(stat_wind_angle * pi / 180)
                  stat_wind_rad * sin(pi / stat_nof_slots) stat_wind_rad * cos(pi / stat_nof_slots)
                  stat_base_rad * sin(pi / stat_nof_slots) stat_base_rad * cos(pi / stat_nof_slots)
                  stat_wind_rad * sin(pi / stat_nof_slots) stat_wind_rad * cos(pi / stat_nof_slots)
                  stat_wind_rad * sin(2 * pi / stat_nof_slots - stat_wind_angle * pi / 180) stat_wind_rad * cos(2 * pi / stat_nof_slots - stat_wind_angle * pi / 180)];
% magnet shape
rot_mag_shape = [ rot_mag_width / 2     rot_inner_rad;
                  rot_mag_width / 2     rot_inner_rad + rot_mag_thick;
                 -rot_mag_width / 2     rot_inner_rad + rot_mag_thick;
                 -rot_mag_width / 2     rot_inner_rad];
% inner radius of rotor ring
rot_ring_inner_rad = sqrt((rot_inner_rad + rot_mag_thick)^2 + (rot_mag_width / 2)^2);
% position of next slot
stat_base_n_x = stat_base_rad * sin((stat_slot_edge_angle + 2 * stat_slot_base_angle) / 180 * pi);
stat_base_n_y = stat_base_rad * cos((stat_slot_edge_angle + 2 * stat_slot_base_angle) / 180 * pi);
