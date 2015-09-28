%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Magnetic simulation of a BLDC motor using FEMM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform torque simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear previous selection
mi_clearselected();

% Iterate through rotor rotation and current phase
for i = rot_range
    for phi = i_range
        % Set current for actual phase phi
        il1 = current * sin((phi      ) * pi / 180);
        il2 = current * sin((phi + 120) * pi / 180);
        il3 = current * sin((phi - 120) * pi / 180);
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
        end
    end
    mi_selectgroup(group_rotor);
    mi_moverotate(0, 0, res_rot);
end

% rotate rotor back to original position
mi_selectgroup(group_rotor);
mi_moverotate(0, 0, -(360 * 3 / stat_nof_slots));

% plot simulation results
if simulate == 1
    save('bldc_torque.mat', 'torque');
    figure(1);
    mesh(i_range, rot_range, torque);
    title('Torque');
    xlabel('phase angle [^\circ]');
    ylabel('rotation [^\circ]');
    zlabel('torque [Nm]');
    print -dpdf 'bldc_torque'
    close all;
end

% Clear previous selection
mi_clearselected();
