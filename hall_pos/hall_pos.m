%% Magnetic simulation of a BLDC motor using FEMM

%% Clear variables and windows
clear;
clc;
close all;

%% Setting up plots
graphics_toolkit('gnuplot');

% Switch to either control simulation
% 0 -> Test sequence only
% 1 -> Perform actual simulations
simulate = 1;
debug = 0;

%% Constants for defining dimensions
mag_x = 2;
mag_y = 2;
mag_thick = 4;
nof_mag = 50;

iron_thick = 2;

mat_air     = 'Air';
mat_iron    = 'Vanadium Permedur';
mat_magnet  = 'NdFeB 52 MGOe';

% Groups for stator and rotor
group_magnet = 1;
group_iron = 2;
group_air = 3;

% Simulation setup
res = 0.1;
dist_max = 2;

%% Calculate needed variables
mag_shape = [0 0; mag_x 0; mag_x mag_y; 0 mag_y];
mag_length = mag_x * nof_mag;
iron_shape = [0 0; mag_length 0; mag_length -iron_thick; 0 -iron_thick];

%% Open FEMM
openfemm;

% Create a new magnetics problem
newdocument(0);

% Set up problem
mi_probdef(0, 'millimeters', 'planar', 1e-8, mag_thick, '30', 0);

%% Materials
% Materials from Library
% ----------------------
mi_getmaterial(mat_air);
mi_getmaterial(mat_magnet);
mi_getmaterial(mat_iron);

%% Draw magnets

% Iron
% ----
mi_drawpolygon(iron_shape);
for i = 1:length(iron_shape)
    mi_selectnode(iron_shape(i,:));
endfor
mi_setgroup(group_iron);
for i = 1:length(iron_shape) - 1
    mi_selectsegment((iron_shape(i,:) + iron_shape(i+1,:)) / 2);
endfor
mi_selectsegment((iron_shape(1,:) + iron_shape(length(iron_shape),:)) / 2);
mi_setgroup(group_iron);
mi_clearselected();

% Magnets
% -------
mi_drawpolygon(mag_shape);
for i = 1:length(mag_shape)
    mi_selectnode(mag_shape(i,:));
endfor
mi_setgroup(group_magnet);
for i = 1:length(mag_shape)
    mi_selectnode(mag_shape(i,:));
endfor
mi_copytranslate2(mag_x, 0, nof_mag - 1, 0);
for i = 1:length(mag_shape) - 1
    mi_selectsegment((mag_shape(i,:) + mag_shape(i+1,:)) / 2);
endfor
mi_selectsegment((mag_shape(1,:) + mag_shape(length(mag_shape),:)) / 2);
mi_setgroup(group_magnet);
for i = 1:length(mag_shape) - 1
    mi_selectsegment((mag_shape(i,:) + mag_shape(i+1,:)) / 2);
endfor
mi_selectsegment((mag_shape(1,:) + mag_shape(length(mag_shape),:)) / 2);
mi_copytranslate2(mag_x, 0, nof_mag - 1, 1);
mi_clearselected();

% Boundary around magnets
% -----------------------
mi_zoomnatural();
mi_drawarc(mag_length / 2,  mag_length * 10, mag_length / 2, -mag_length * 10, 180, 1);
mi_drawarc(mag_length / 2, -mag_length * 10, mag_length / 2,  mag_length * 10, 180, 1);

% Block labels
% ------------
% Air
mi_addblocklabel(-1, 0);
mi_selectlabel(-1, 0);
mi_setblockprop(mat_air, 1, 0, '<none>', 0, group_air, 0)
mi_clearselected();
% Iron
mi_addblocklabel(mag_length / 2, -iron_thick / 2);
mi_selectlabel(mag_length / 2, -iron_thick / 2);
mi_setblockprop(mat_iron, 1, 0, '<none>', 0, group_iron, 0)
mi_clearselected();
% Magnets
mi_addblocklabel(mag_x / 2, mag_y / 2);
mi_selectlabel(mag_x / 2, mag_y / 2);
mi_setgroup(group_magnet);
mi_selectlabel(mag_x / 2, mag_y / 2);
mi_copytranslate2(mag_x, 0, nof_mag - 1, 2);
for i = 0:(nof_mag - 1)
    mi_selectlabel(mag_x * i + mag_x / 2, mag_y / 2);
    if (mod(i, 2))
        dir = 90;
    else
        dir = 270;
    endif
    mi_setblockprop(mat_magnet, 1, 0, '<none>', dir, group_magnet, 0);
    mi_clearselected();
endfor

%% Zoom view to current problem
%mi_zoomnatural();

%% Perform simulation
mi_clearselected();
name = strcat('hall_pos');
mi_saveas(strcat(name, '.fem'));
if simulate == 1
    mi_createmesh();
    mi_analyze(1);
    mi_loadsolution();
    mo_zoom(0, 0, mag_length, mag_y);
    mo_showdensityplot(1, 0, 3, 0, 'mag');
    %mo_close();
endif
if simulate == 1
    for dist_int = 1:1:(dist_max/res) + 1
        dist = (dist_int - 1) * res;
        name_out = strcat(name, '_d', num2str(dist));
        mo_addcontour(-10, mag_y + dist);
        mo_addcontour(mag_length + 10, mag_y + dist);
        % B
        mo_makeplot(2, 1000, strcat(name_out, '_Bn.txt'), 1);
        mo_makeplot(3, 1000, strcat(name_out, '_Bt.txt'), 1);
        % H
        mo_makeplot(5, 1000, strcat(name_out, '_Hn.txt'), 1);
        mo_makeplot(6, 1000, strcat(name_out, '_Ht.txt'), 1);
        mo_clearcontour();
        % import results
        tmp = load(strcat(name_out, '_Bn.txt'));
        bn(:,dist_int) = tmp(:,2);
        tmp = load(strcat(name_out, '_Bt.txt'));
        bt(:,dist_int) = tmp(:,2);
        tmp = load(strcat(name_out, '_Hn.txt'));
        hn(:,dist_int) = tmp(:,2);
        tmp = load(strcat(name_out, '_Ht.txt'));
        ht(:,dist_int) = tmp(:,2);
    endfor
endif
if simulate == 1
    save('bn.mat', 'bn');
    save('bt.mat', 'bt');
    save('hn.mat', 'hn');
    save('ht.mat', 'ht');

    figure(1);
    mesh(bn);
    title('B_n');
    xlabel('distance');
    ylabel('y');
    zlabel('B_n');
    print -dpdf 'hall_pos_Bn';

    figure(2);
    mesh(bt);
    title('B_t');
    xlabel('distance');
    ylabel('y');
    zlabel('B_t');
    print -dpdf 'hall_pos_Bt';

    figure(3);
    mesh(hn);
    title('H_n');
    xlabel('distance');
    ylabel('y');
    zlabel('H_n');
    print -dpdf 'hall_pos_Hn';

    figure(4);
    mesh(ht);
    title('H_t');
    xlabel('distance');
    ylabel('y');
    zlabel('H_t');
    print -dpdf 'hall_pos_Ht';

    close all;
endif

%% close FEMM
closefemm;
