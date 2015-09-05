%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Magnetic simulation of a BLDC motor using FEMM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear variables and windows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants for defining motor dimensions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rot_nof_poles = 10;
rot_nof_mag_per_pole = 2;
stat_nof_slots = 15;

rot_mag_thick = 3;
rot_mag_width = 9;
rot_inner_rad = 30.5;
rot_outer_rad = 35;

stat_out_rad = 30;
stat_base_rad = 20;
stat_inner_rad = 16;

stat_slot_width = 4;
stat_head_width = 8;
%                   Width   Radius
stat_slot_shape =  [0.      0.      ;
                    0.      0.88    ;
                    0.1     0.9     ;
                    0.95    0.9     ;
                    1.      1.      ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate needed variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stator edge
stat_head_edge = sqrt((stat_out_rad)^2 - (stat_head_width/2)^2);
% scale stator slot
stat_slot_shape_r(:,1) = (stat_slot_shape(:,1) .* (stat_head_width - stat_slot_width) .+ stat_slot_width) ./ 2;
stat_slot_shape_r(:,2) = stat_slot_shape(:,2) .* (stat_head_edge - stat_base_rad) .+ stat_base_rad;
stat_slot_shape_l(:,1) = -stat_slot_shape_r(:,1);
stat_slot_shape_l(:,2) = stat_slot_shape_r(:,2);
% stator head angle
stat_head_angle = asin(stat_head_width/2 / stat_out_rad) * 180 / pi;
% stator base angle
stat_slot_edge_angle = atan(stat_slot_shape_r(1,1) / stat_slot_shape_r(1,2)) *180 / pi
stat_slot_base_angle = 360 / stat_nof_slots / 2 - stat_slot_edge_angle
% magnet shape
rot_mag_shape = [ rot_mag_width / 2     rot_inner_rad;
                  rot_mag_width / 2     rot_inner_rad + rot_mag_thick;
                 -rot_mag_width / 2     rot_inner_rad + rot_mag_thick;
                 -rot_mag_width / 2     rot_inner_rad];
% inner radius of rotor ring
rot_ring_inner_rad = sqrt((rot_inner_rad + rot_mag_thick)^2 + (rot_mag_width / 2)^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open FEMM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
openfemm;

% Create a new magnetics problem
newdocument(0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Draw motor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rotor
% -----
% rotor ring
mi_drawarc(0, rot_outer_rad, 0, -rot_outer_rad, 180, 1);
mi_drawarc(0, -rot_outer_rad, 0, rot_outer_rad, 180, 1);
mi_drawarc(0, rot_ring_inner_rad, 0, -rot_ring_inner_rad, 180, 1);
mi_drawarc(0, -rot_ring_inner_rad, 0, rot_ring_inner_rad, 180, 1);
% magnets
mi_drawpolygon(rot_mag_shape);
for i = 1:length(rot_mag_shape)
    mi_selectnode(rot_mag_shape(i,:));
endfor
mi_copyrotate2(0, 0, 360 / rot_nof_poles / rot_nof_mag_per_pole, rot_nof_poles * rot_nof_mag_per_pole - 1, 0);
for i = 1:length(rot_mag_shape) - 1
    mi_selectsegment((rot_mag_shape(i,:) + rot_mag_shape(i+1,:)) / 2);
endfor
mi_selectsegment((rot_mag_shape(1,:) + rot_mag_shape(length(rot_mag_shape),:)) / 2);
mi_copyrotate2(0, 0, 360 / rot_nof_poles / rot_nof_mag_per_pole, rot_nof_poles * rot_nof_mag_per_pole - 1, 1);
mi_clearselected();
mi_selectnode(0, rot_ring_inner_rad);
if (mod(rot_nof_poles * rot_nof_mag_per_pole, 2) == 0)
    mi_selectnode(0, -rot_ring_inner_rad);
endif
for i = 1:(rot_nof_poles * rot_nof_mag_per_pole / 2) - 1
    mi_selectarcsegment(sin(2*pi/(rot_nof_poles * rot_nof_mag_per_pole) * i), cos(2*pi/(rot_nof_poles * rot_nof_mag_per_pole) * i));
    mi_selectarcsegment(-sin(2*pi/(rot_nof_poles * rot_nof_mag_per_pole) * i), cos(2*pi/(rot_nof_poles * rot_nof_mag_per_pole) * i));
endfor
mi_deleteselected();

% Stator
% ------
% inner circle
mi_drawarc(0, stat_inner_rad, 0, -stat_inner_rad, 180, 1);
mi_drawarc(0, -stat_inner_rad, 0, stat_inner_rad, 180, 1);
% slots
mi_drawpolyline(stat_slot_shape_r);
mi_drawpolyline(stat_slot_shape_l);
mi_addarc(stat_head_width/2, stat_head_edge, -stat_head_width/2, stat_head_edge, stat_head_angle*2, 1);
for i = 1:length(stat_slot_shape)
    mi_selectnode(stat_slot_shape_r(i,:));
    mi_selectnode(stat_slot_shape_l(i,:));
endfor
mi_copyrotate2(0, 0, 360 / stat_nof_slots, stat_nof_slots - 1, 0);
for i = 1:length(stat_slot_shape) - 1
    mi_selectsegment((stat_slot_shape_r(i,:) + stat_slot_shape_r(i+1,:)) / 2);
    mi_selectsegment((stat_slot_shape_l(i,:) + stat_slot_shape_l(i+1,:)) / 2);
endfor
mi_copyrotate2(0, 0, 360 / stat_nof_slots, stat_nof_slots - 1, 1);
mi_selectarcsegment(0, stat_out_rad);
mi_copyrotate2(0, 0, 360 / stat_nof_slots, stat_nof_slots - 1, 3);
x2 = 6.3
y2 = 19.07
mi_addarc(x2, y2, stat_slot_shape_r(1,1), stat_slot_shape_r(1,2), stat_slot_base_angle * 2, 1);
mi_selectarcsegment((x2 + stat_slot_shape_r(1,1)) / 2, (y2 + stat_slot_shape_r(1,2)) / 2);
mi_copyrotate2(0, 0, 360 / stat_nof_slots, stat_nof_slots - 1, 3);

% Boundary around motor
% ---------------------
mi_zoomnatural();
mi_drawarc(0, rot_outer_rad * 3, 0, -rot_outer_rad * 3, 180, 1);
mi_drawarc(0, -rot_outer_rad * 3, 0, rot_outer_rad * 3, 180, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zoom view to current problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mi_zoomnatural();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close FEMM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sleep(10);
%closefemm;
