%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Magnetic simulation of a BLDC motor using FEMM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read data from FEMM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read number of nodes from FEMM
nof_nodes = mo_numnodes();
%nof_nodes = 1000;
% read coordinates for every node from FEMM
for i = 1:nof_nodes
    xy(i,:) = mo_getnode(i);
end
% Calculate radius of every node
rad = sqrt(xy(:,1).^2 + xy(:,2).^2);
% remove nodes with coordinates outside of border (problems caused by rounding errors)
xyfilt = xy(find(rad < (3 * rot_outer_rad)),:);
% read flux density for every node from FEMM
b = mo_getb(xyfilt);
