close all; clc;
disp('Running script that drives pwld');

global verbose
verbose = true;
verbose = false;

% data.geofile = '';
% geo = 'shestakov_quad_L1_nc8_emb4_a0.1.txt';
% geo = 'random_poly_mesh_L1_n64_a0.965.txt' ;
% geo = 'random_poly_mesh_L1_n64_a0.9.txt'
% geo = 'z_mesh_quad_L1_n9_a0.05.txt'
% geo = 'z_mesh_quad_L1_n6_a0.05.txt'
% geo = 'z_mesh_quad_L1_n5_a0.05.txt'
% geo = 'shestakov_poly_mesh_L1_nc4_a0.02.txt';
% geo = 'z_mesh_poly_L1_n5_a0.05.txt';
geo = 'z_mesh_quad_L1_n6_a0.35.txt';

data.geofile = sprintf('%s%s','..\geom_codes\figs\',geo);

data.logi_mms = true;
data.mms_type = 1;
if ~data.logi_mms
    data_pbtype = 'linear'
end

data.max_ref_cycles = 1;
data.ref_threshold  = 0; % 0=uniform refinement
data.logi_amr_output  = false; % 0=uniform refinement

data.logi_plot           = true;
data.logi_plot_err_i     = false;
data.generate_vtk_output = true;

% result_basename = 'results\z_mesh_poly_n5_a0.05_';
result_basename = 'results\z_mesh_quad_n6_a0.35_';
% result_basename = 'results\toto_';

if(data.logi_mms)
    result_basename = sprintf('%s%s%c',result_basename,'mms_',int2str(data.mms_type));
else
    result_basename = sprintf('%s%s',result_basename,data_pbtype);
end
result_basename
data.vtk_basename        = result_basename;

data.save_workspace      = true;
data.workspace_name      = strcat(result_basename,'.mat');

pwld_solve_problem(data);


