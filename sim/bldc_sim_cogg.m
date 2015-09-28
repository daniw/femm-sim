%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Magnetic simulation of a BLDC motor using FEMM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform cogging torque simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear previous selection
mi_clearselected();

% Iterate through rotor rotation
for i = cogg_range
    name = strcat('bldc_r', num2str(i - res_cogg), '_cogg');
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
        cogg_torque(i / res_cogg) = mo_blockintegral(22);
        mo_close();
    end
    mi_selectgroup(group_rotor);
    mi_moverotate(0, 0, res_cogg);
end

% rotate rotor back to original position
mi_selectgroup(group_rotor);
mi_moverotate(0, 0, -(360 * 3 / stat_nof_slots));

% plot simulation results
if simulate == 1
    save('bldc_cogging_torque.mat', 'cogg_torque');
    figure(1);
    plot(cogg_range, cogg_torque);
    title('Cogging torque');
    xlabel('rotation [^\circ]');
    ylabel('cogging torque [Nm]');
    print -dpdf 'bldc_cogging_torque';
    close all;
end

% Clear previous selection
mi_clearselected();
