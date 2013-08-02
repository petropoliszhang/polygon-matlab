close all; clc;
disp('Running script that drives pwld');

global verbose
verbose = true;

% data.geofile = '';
geo = 'shestakov_quad_L1_nc8_emb4_a0.1.txt';
data.geofile = sprintf('%s%s','..\geom_codes\figs\',geo);

data.logi_mms = true;
data.mms_type = 1;
if ~data.logi_mms
    data_pbtype = 'linear'
end

data.max_ref_cycles = 1;
data.ref_threshold  = 0; % 0=uniform refinement

data.logi_plot           = true;
data.logi_plot_err_i     = false;
data.generate_vtk_output = false;

result_basename = 'results\reg_quad_';
if(data.logi_mms)
    result_basename = sprintf('%s%s%c',result_basename,'mms_',int2str(data.mms_type));
else
    result_basename = sprintf('%s%s',result_basename,data_pbtype);
end
result_basename
data.vtk_basename        = result_basename;

data.save_workspace      = false;
data.workspace_name      = result_basename;

pwld_solve_problem(data);


