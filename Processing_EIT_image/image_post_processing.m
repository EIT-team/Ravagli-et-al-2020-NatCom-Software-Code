
clc;
clear all;
close all;

%% Load image and re-compute centre of each voxel in chopped mesh
filename='IMG_RECON_22_EIT_RL_T_1mA_lite_BW2K.mat'; % Change image file here
load(filename);

%Compute centre of each voxel
Mesh_hex_cut.centre = ( Mesh_hex_cut.Nodes(Mesh_hex_cut.Hex(:,1),:)+Mesh_hex_cut.Nodes(Mesh_hex_cut.Hex(:,2),:)+Mesh_hex_cut.Nodes(Mesh_hex_cut.Hex(:,3),:)+Mesh_hex_cut.Nodes(Mesh_hex_cut.Hex(:,4),:)+...
Mesh_hex_cut.Nodes(Mesh_hex_cut.Hex(:,5),:)+Mesh_hex_cut.Nodes(Mesh_hex_cut.Hex(:,6),:)+Mesh_hex_cut.Nodes(Mesh_hex_cut.Hex(:,7),:)+Mesh_hex_cut.Nodes(Mesh_hex_cut.Hex(:,8),:))/8;

%% Choose frame of peak dZ/dSigma
T_peak=58;
% T_peak=T_peak+[-50:1:50]; % Averaging over frames around peak
sigma_pp=mean(sigma(:,T_peak),2); % Take closest slice to reference ring, best image here (see Ravagli et al. 2019)

%% Find all unique slices along the x-axis and choose one
[x_pos, ia, ib] = unique(Mesh_hex_cut.centre(:,1),'rows');
slices = accumarray(ib, find(ib), [], @(rows){rows});
nSlices=length(slices); 
ind_slice=slices{1}; % Closest slice to reference ring

sigma_pp=sigma_pp(ind_slice,1);
Mesh_hex_cut.Hex=Mesh_hex_cut.Hex(ind_slice,:);
Mesh_hex_cut.centre=Mesh_hex_cut.centre(ind_slice,:);

%% Median filtering for removal of "spike" voxels

dist_th_median=41e-6;   % Radius in voxels=1

sigma_temp=zeros(size(sigma_pp));
for iVoxel=1:size(sigma_pp,1)    
    c=Mesh_hex_cut.centre(iVoxel,:);
    ind= abs(Mesh_hex_cut.centre(:,1)-c(1))<dist_th_median & abs(Mesh_hex_cut.centre(:,2)-c(2))<dist_th_median & abs(Mesh_hex_cut.centre(:,3)-c(3))<dist_th_median;
    sigma_temp(iVoxel,1)=median(sigma_pp(ind,1));    
end

sigma_pp=sigma_temp;

%% Mean filtering for smoothing

dist_th_mean=121e-6; % Radius in voxels=4

sigma_temp=zeros(size(sigma_pp));
for iVoxel=1:size(sigma_pp,1)    
    c=Mesh_hex_cut.centre(iVoxel,:);
    ind= abs(Mesh_hex_cut.centre(:,1)-c(1))<dist_th_mean & abs(Mesh_hex_cut.centre(:,2)-c(2))<dist_th_mean & abs(Mesh_hex_cut.centre(:,3)-c(3))<dist_th_mean;
    sigma_temp(iVoxel,1)=mean(sigma_pp(ind,1));    
end

sigma_pp=sigma_temp;

%% Normalize in range [0-1]
sigma_pp=(sigma_pp-min(sigma_pp))/(max(sigma_pp)-min(sigma_pp));

%% Write the resulting dSigma distribution to VTK
writeVTKcell_hex([ filename(1:end-4)],Mesh_hex_cut.Hex,Mesh_hex_cut.Nodes,sigma_pp,'sigma');































