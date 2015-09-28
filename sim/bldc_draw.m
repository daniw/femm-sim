%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Magnetic simulation of a BLDC motor using FEMM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
end
mi_setgroup(group_rotor);
for i = 1:length(rot_mag_shape)
    mi_selectnode(rot_mag_shape(i,:));
end
mi_copyrotate2(0, 0, 360 / rot_nof_poles / rot_nof_mag_per_pole, rot_nof_poles * rot_nof_mag_per_pole - 1, 0);
for i = 1:length(rot_mag_shape) - 1
    mi_selectsegment((rot_mag_shape(i,:) + rot_mag_shape(i+1,:)) / 2);
end
mi_selectsegment((rot_mag_shape(1,:) + rot_mag_shape(length(rot_mag_shape),:)) / 2);
mi_setgroup(group_rotor);
for i = 1:length(rot_mag_shape) - 1
    mi_selectsegment((rot_mag_shape(i,:) + rot_mag_shape(i+1,:)) / 2);
end
mi_selectsegment((rot_mag_shape(1,:) + rot_mag_shape(length(rot_mag_shape),:)) / 2);
mi_copyrotate2(0, 0, 360 / rot_nof_poles / rot_nof_mag_per_pole, rot_nof_poles * rot_nof_mag_per_pole - 1, 1);
mi_clearselected();
mi_selectnode(0, rot_ring_inner_rad);
if (mod(rot_nof_poles * rot_nof_mag_per_pole, 2) == 0)
    mi_selectnode(0, -rot_ring_inner_rad);
end
for i = 1:(rot_nof_poles * rot_nof_mag_per_pole / 2) - 1
    mi_selectarcsegment(sin(2*pi/(rot_nof_poles * rot_nof_mag_per_pole) * i), cos(2*pi/(rot_nof_poles * rot_nof_mag_per_pole) * i));
    mi_selectarcsegment(-sin(2*pi/(rot_nof_poles * rot_nof_mag_per_pole) * i), cos(2*pi/(rot_nof_poles * rot_nof_mag_per_pole) * i));
end
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
end
mi_setgroup(group_stator);
for i = 1:length(stat_slot_shape)
    mi_selectnode(stat_slot_shape_r(i,:));
    mi_selectnode(stat_slot_shape_l(i,:));
end
mi_copyrotate2(0, 0, 360 / stat_nof_slots, stat_nof_slots - 1, 0);
for i = 1:length(stat_slot_shape) - 1
    mi_selectsegment((stat_slot_shape_r(i,:) + stat_slot_shape_r(i+1,:)) / 2);
    mi_selectsegment((stat_slot_shape_l(i,:) + stat_slot_shape_l(i+1,:)) / 2);
end
mi_setgroup(group_stator);
for i = 1:length(stat_slot_shape) - 1
    mi_selectsegment((stat_slot_shape_r(i,:) + stat_slot_shape_r(i+1,:)) / 2);
    mi_selectsegment((stat_slot_shape_l(i,:) + stat_slot_shape_l(i+1,:)) / 2);
end
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
end
mi_setgroup(group_stator);
for i = [1 2 4]
    mi_selectsegment((stat_wind_cont(i,:) + stat_wind_cont(i+1,:)) / 2);
end
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
    end
    mi_setblockprop(mat_magnet, 1, 0, '<none>', dir, group_rotor, 0);
    mi_clearselected();
end
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
    end
end
mi_clearselected();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zoom view to current problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mi_zoomnatural();
