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
% Setting up plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
graphics_toolkit('gnuplot');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Switch to either control simulation
% 0 -> Test sequence only
% 1 -> Perform actual simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simulate = 1;
debug = 0;

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
mat_wind    = '1mm';

current = 5;

% Resolution for simulation iterations
res_clogg = 1;
res_rot = 5;
res_el  = 20;

% Groups for stator and rotor
group_stator = 1;
group_rotor = 2;
group_air = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate needed variables
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open FEMM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
openfemm;

% Create a new magnetics problem
newdocument(0);

% Set up problem
mi_probdef(0, 'millimeters', 'planar', 1e-8, mot_thickness, '30', 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Materials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Materials from Library
% ----------------------
mi_getmaterial(mat_air);
mi_getmaterial(mat_magnet);
mi_getmaterial(mat_stat);
if mat_rot != mat_stat
    mi_getmaterial(mat_rot);
endif
mi_getmaterial(mat_wind);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Circuits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mi_addcircprop('L1+', 0, 1);
mi_addcircprop('L1-', 0, 1);
mi_addcircprop('L2+', 0, 1);
mi_addcircprop('L2-', 0, 1);
mi_addcircprop('L3+', 0, 1);
mi_addcircprop('L3-', 0, 1);

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
mi_selectnode(0, rot_outer_rad);
mi_selectnode(0, -rot_outer_rad);
mi_selectnode(0, rot_ring_inner_rad);
mi_selectnode(0, -rot_ring_inner_rad);
mi_selectarcsegment(rot_outer_rad, 0);
mi_selectarcsegment(-rot_outer_rad, 0);
mi_selectarcsegment(rot_ring_inner_rad, 0);
mi_selectarcsegment(-rot_ring_inner_rad, 0);
mi_setgroup(group_rotor);
% magnets
mi_drawpolygon(rot_mag_shape);
for i = 1:length(rot_mag_shape)
    mi_selectnode(rot_mag_shape(i,:));
endfor
mi_setgroup(group_rotor);
for i = 1:length(rot_mag_shape)
    mi_selectnode(rot_mag_shape(i,:));
endfor
mi_copyrotate2(0, 0, 360 / rot_nof_poles / rot_nof_mag_per_pole, rot_nof_poles * rot_nof_mag_per_pole - 1, 0);
for i = 1:length(rot_mag_shape) - 1
    mi_selectsegment((rot_mag_shape(i,:) + rot_mag_shape(i+1,:)) / 2);
endfor
mi_selectsegment((rot_mag_shape(1,:) + rot_mag_shape(length(rot_mag_shape),:)) / 2);
mi_setgroup(group_rotor);
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
mi_selectnode(0, stat_inner_rad);
mi_selectnode(0, -stat_inner_rad);
mi_selectarcsegment(stat_inner_rad, 0);
mi_selectarcsegment(-stat_inner_rad, 0);
mi_setgroup(group_stator);
% slots
mi_drawpolyline(stat_slot_shape_r);
mi_drawpolyline(stat_slot_shape_l);
mi_addarc(stat_head_width/2, stat_head_edge, -stat_head_width/2, stat_head_edge, stat_head_angle*2, 1);
for i = 1:length(stat_slot_shape)
    mi_selectnode(stat_slot_shape_r(i,:));
    mi_selectnode(stat_slot_shape_l(i,:));
endfor
mi_setgroup(group_stator);
for i = 1:length(stat_slot_shape)
    mi_selectnode(stat_slot_shape_r(i,:));
    mi_selectnode(stat_slot_shape_l(i,:));
endfor
mi_copyrotate2(0, 0, 360 / stat_nof_slots, stat_nof_slots - 1, 0);
for i = 1:length(stat_slot_shape) - 1
    mi_selectsegment((stat_slot_shape_r(i,:) + stat_slot_shape_r(i+1,:)) / 2);
    mi_selectsegment((stat_slot_shape_l(i,:) + stat_slot_shape_l(i+1,:)) / 2);
endfor
mi_setgroup(group_stator);
for i = 1:length(stat_slot_shape) - 1
    mi_selectsegment((stat_slot_shape_r(i,:) + stat_slot_shape_r(i+1,:)) / 2);
    mi_selectsegment((stat_slot_shape_l(i,:) + stat_slot_shape_l(i+1,:)) / 2);
endfor
mi_copyrotate2(0, 0, 360 / stat_nof_slots, stat_nof_slots - 1, 1);
mi_selectarcsegment(0, stat_out_rad);
mi_setgroup(group_stator);
mi_selectarcsegment(0, stat_out_rad);
mi_copyrotate2(0, 0, 360 / stat_nof_slots, stat_nof_slots - 1, 3);
mi_addarc(stat_base_n_x, stat_base_n_y, stat_slot_shape_r(1,1), stat_slot_shape_r(1,2), stat_slot_base_angle * 2, 1);
mi_selectarcsegment((stat_base_n_x + stat_slot_shape_r(1,1)) / 2, (stat_base_n_y + stat_slot_shape_r(1,2)) / 2);
mi_setgroup(group_stator);
mi_selectarcsegment((stat_base_n_x + stat_slot_shape_r(1,1)) / 2, (stat_base_n_y + stat_slot_shape_r(1,2)) / 2);
mi_copyrotate2(0, 0, 360 / stat_nof_slots, stat_nof_slots - 1, 3);
% Windings
mi_drawpolyline(stat_wind_cont);
mi_selectnode(stat_wind_cont(2,:));
mi_selectnode(stat_wind_cont(3,:));
mi_setgroup(group_stator);
mi_selectnode(stat_wind_cont(2,:));
mi_selectnode(stat_wind_cont(3,:));
mi_copyrotate2(0, 0, 360 / stat_nof_slots, stat_nof_slots - 1, 0);
for i = [1 2 4]
    mi_selectsegment((stat_wind_cont(i,:) + stat_wind_cont(i+1,:)) / 2);
endfor
mi_setgroup(group_stator);
for i = [1 2 4]
    mi_selectsegment((stat_wind_cont(i,:) + stat_wind_cont(i+1,:)) / 2);
endfor
mi_copyrotate2(0, 0, 360 / stat_nof_slots, stat_nof_slots - 1, 1);

% Boundary around motor
% ---------------------
mi_zoomnatural();
mi_drawarc(0, rot_outer_rad * 3, 0, -rot_outer_rad * 3, 180, 1);
mi_drawarc(0, -rot_outer_rad * 3, 0, rot_outer_rad * 3, 180, 1);

% Block labels
% ------------
% Air
mi_addblocklabel(0, 0);
mi_selectlabel(0, 0);
mi_setblockprop(mat_air, 1, 0, '<none>', 0, group_air, 0)
mi_addblocklabel(0, 2 * rot_outer_rad);
mi_selectlabel(0, 2 * rot_outer_rad);
mi_setblockprop(mat_air, 1, 0, '<none>', 0, group_air, 0)
mi_addblocklabel(0, -stat_head_edge);
mi_selectlabel(0, -stat_head_edge);
mi_setblockprop(mat_air, 1, 0, '<none>', 0, group_air, 0)
mi_clearselected();
% Stator iron
mi_addblocklabel(0, (stat_inner_rad + stat_base_rad) / 2);
mi_selectlabel(0, (stat_inner_rad + stat_base_rad) / 2);
mi_setgroup(group_stator);
mi_selectlabel(0, (stat_inner_rad + stat_base_rad) / 2);
mi_setblockprop(mat_stat, 1, 0, '<none>', 0, group_stator, 0)
mi_clearselected();
% Rotor iron
mi_addblocklabel(0, (rot_outer_rad + rot_ring_inner_rad) / 2);
mi_selectlabel(0, (rot_outer_rad + rot_ring_inner_rad) / 2);
mi_setgroup(group_rotor);
mi_selectlabel(0, (rot_outer_rad + rot_ring_inner_rad) / 2);
mi_setblockprop(mat_rot, 1, 0, '<none>', 0, group_rotor, 0)
mi_clearselected();
% Magnets
mi_addblocklabel(0, rot_inner_rad + (rot_mag_thick / 2));
mi_selectlabel(0, rot_inner_rad + (rot_mag_thick / 2));
mi_setgroup(group_rotor);
mi_selectlabel(0, rot_inner_rad + (rot_mag_thick / 2));
mi_copyrotate2(0, 0, 360 / rot_nof_poles / rot_nof_mag_per_pole, (rot_nof_poles * rot_nof_mag_per_pole) - 1, 2)
for i = 0:(rot_nof_poles * rot_nof_mag_per_pole - 1)
    mi_selectlabel((rot_inner_rad + (rot_mag_thick / 2)) * sin(i / (rot_nof_poles * rot_nof_mag_per_pole) * 2 * pi), -(rot_inner_rad + (rot_mag_thick / 2)) * cos(i / (rot_nof_poles * rot_nof_mag_per_pole) * 2 * pi));
    if (mod(floor(i / rot_nof_mag_per_pole), 2))
        dir = 90 + (360 / rot_nof_poles / rot_nof_mag_per_pole * i);
    else
        dir = 270 + (360 / rot_nof_poles / rot_nof_mag_per_pole * i);
    endif
    mi_setblockprop(mat_magnet, 1, 0, '<none>', dir, group_rotor, 0);
    mi_clearselected();
endfor
% Windings
for i = 0:stat_nof_slots-1
    wind_angle_p = 2 * pi / stat_nof_slots * i + (180 / stat_nof_slots + stat_wind_angle) / 2 * pi / 180;
    wind_angle_n = 2 * pi / stat_nof_slots * i - (180 / stat_nof_slots + stat_wind_angle) / 2 * pi / 180;
    mi_addblocklabel([(stat_base_rad + stat_wind_rad) / 2 * sin(wind_angle_p)
                      (stat_base_rad + stat_wind_rad) / 2 * cos(wind_angle_p)]);
    mi_addblocklabel([(stat_base_rad + stat_wind_rad) / 2 * sin(wind_angle_n)
                      (stat_base_rad + stat_wind_rad) / 2 * cos(wind_angle_n)]);
    mi_selectlabel([(stat_base_rad + stat_wind_rad) / 2 * sin(wind_angle_p)
                    (stat_base_rad + stat_wind_rad) / 2 * cos(wind_angle_p)]);
    mi_selectlabel([(stat_base_rad + stat_wind_rad) / 2 * sin(wind_angle_n)
                    (stat_base_rad + stat_wind_rad) / 2 * cos(wind_angle_n)]);
    mi_setgroup(group_stator);
    if mod(i, 3) == 0
        mi_selectlabel([(stat_base_rad + stat_wind_rad) / 2 * sin(wind_angle_p)
                        (stat_base_rad + stat_wind_rad) / 2 * cos(wind_angle_p)]);
        mi_setblockprop(mat_wind, 1, 0, 'L1+', 0, group_stator, stat_nof_wdg)
        mi_clearselected();
        mi_selectlabel([(stat_base_rad + stat_wind_rad) / 2 * sin(wind_angle_n)
                        (stat_base_rad + stat_wind_rad) / 2 * cos(wind_angle_n)]);
        mi_setblockprop(mat_wind, 1, 0, 'L1-', 0, group_stator, stat_nof_wdg)
        mi_clearselected();
    elseif mod(i, 3) == 1
        mi_selectlabel([(stat_base_rad + stat_wind_rad) / 2 * sin(wind_angle_p)
                        (stat_base_rad + stat_wind_rad) / 2 * cos(wind_angle_p)]);
        mi_setblockprop(mat_wind, 1, 0, 'L2+', 0, group_stator, stat_nof_wdg)
        mi_clearselected();
        mi_selectlabel([(stat_base_rad + stat_wind_rad) / 2 * sin(wind_angle_n)
                        (stat_base_rad + stat_wind_rad) / 2 * cos(wind_angle_n)]);
        mi_setblockprop(mat_wind, 1, 0, 'L2-', 0, group_stator, stat_nof_wdg)
        mi_clearselected();
    elseif mod(i, 3) == 2
        mi_selectlabel([(stat_base_rad + stat_wind_rad) / 2 * sin(wind_angle_p)
                        (stat_base_rad + stat_wind_rad) / 2 * cos(wind_angle_p)]);
        mi_setblockprop(mat_wind, 1, 0, 'L3+', 0, group_stator, stat_nof_wdg)
        mi_clearselected();
        mi_selectlabel([(stat_base_rad + stat_wind_rad) / 2 * sin(wind_angle_n)
                        (stat_base_rad + stat_wind_rad) / 2 * cos(wind_angle_n)]);
        mi_setblockprop(mat_wind, 1, 0, 'L3-', 0, group_stator, stat_nof_wdg)
        mi_clearselected();
    endif
endfor
mi_clearselected();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zoom view to current problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mi_zoomnatural();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rotate rotor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mi_clearselected();
loops = 0;
for i = res_clogg:res_clogg:(360 * 3 / stat_nof_slots)
%for i = 1
    name = strcat('bldc_r', num2str(i - res_clogg), '_clogg');
    mi_saveas(strcat(name, '.fem'));
    if simulate == 1
        mi_createmesh();
        mi_analyze(1);
        mi_loadsolution();
        mo_zoom(-rot_outer_rad, -rot_outer_rad, rot_outer_rad, rot_outer_rad);
        mo_showdensityplot(1, 0, 3, 0, 'mag');
        mo_savebitmap(strcat(name, '.bmp'));
        mo_seteditmode('area');
        mo_clearblock();
        mo_groupselectblock(group_stator);
        torque(i / res_clogg) = mo_blockintegral(22);
        mo_close();
    endif
    loops += 1;
    mi_selectgroup(group_rotor);
    mi_moverotate(0, 0, res_clogg);
endfor
mi_selectgroup(group_rotor);
mi_moverotate(0, 0, -(360 * 3 / stat_nof_slots));
if simulate == 1
    save('clogging_torque.mat', 'torque');
    figure(1);
    plot(torque);
    title('Cogging torque');
    xlabel('rotation [^\circ]');
    ylabel('cogging torque [Nm]');
    print -dpdf 'clogging_torque';
    close all;
endif

mi_clearselected();
for i = res_rot:res_rot:(360 * 3 / stat_nof_slots)
%for i = 1
    for phi = res_el:res_el:(360)
        il1 = current * sin((phi      ) * pi / 180);
        il2 = current * sin((phi + 120) * pi / 180);
        il3 = current * sin((phi - 120) * pi / 180);
        is(phi / res_el, 1) = il1;
        is(phi / res_el, 2) = il2;
        is(phi / res_el, 3) = il3;
        is(phi / res_el, 4) = phi;
        mi_setcurrent('L1+',  il1);
        mi_setcurrent('L1-', -il1);
        mi_setcurrent('L2+',  il2);
        mi_setcurrent('L2-', -il2);
        mi_setcurrent('L3+',  il3);
        mi_setcurrent('L3-', -il3);

        name = strcat('bldc_r', num2str(i - res_rot), '_i', num2str(phi - res_el), '_torque');
        mi_saveas(strcat(name, '.fem'));
        if simulate == 1
            mi_createmesh();
            mi_analyze(1);
            mi_loadsolution();
            mo_zoom(-rot_outer_rad, -rot_outer_rad, rot_outer_rad, rot_outer_rad);
            mo_showdensityplot(1, 0, 3, 0, 'mag');
            mo_savebitmap(strcat(name, '.bmp'));
            mo_seteditmode('area');
            mo_clearblock();
            mo_groupselectblock(group_stator);
            torque(i / res_rot, phi / res_el) = mo_blockintegral(22);
            mo_close();
        endif
        loops += 1;
    endfor
    if debug == 1
        figure(2);
        plot(is(:,1:3));
        title('Sine current in windings L_1, L_2, L_3');
        xlabel('\phi [^\circ]');
        ylabel('I_{L_x} [A]');
        legend('I_{L_1}', 'I_{L_2}', 'I_{L_3}');
        print -dpdf 'phase_current';
        close all;
    endif

    mi_selectgroup(group_rotor);
    mi_moverotate(0, 0, res_rot);
endfor
mi_selectgroup(group_rotor);
mi_moverotate(0, 0, -(360 * 3 / stat_nof_slots));
if simulate == 1
    save('torque.mat', 'torque');
    figure(1);
    plot3(torque);
    title('Torque');
    xlabel('rotation [^\circ]');
    ylabel('phase angle [^\circ]');
    zlabel('torque [Nm]');
    print -dpdf 'clogging_torque'
    close all;
endif
loops


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close FEMM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sleep(10);
%closefemm;
