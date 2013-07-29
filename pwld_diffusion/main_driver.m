close all; clc;
disp('Running script that drives pwld');

global verbose
verbose = false;

data.geofile = '';

data.logi_mms = true;
data.mms_type = 1;
if ~data.logi_mms
    data_pbtype = 'linear'
end

data.max_ref_cycles = 7;
data.ref_threshold  = 0; % 0=uniform refinement

data.logi_plot           = true;
data.logi_plot_err_i     = false;
data.generate_vtk_output = true;

result_basename = 'results\reg_quad_';
if(data.logi_mms)
    result_basename = sprintf('%s%s%c',result_basename,'mms_',int2str(data.mms_type));
else
    result_basename = sprintf('%s%s',result_basename,data_pbtype);
end
result_basename
data.vtk_basename        = result_basename;

data.save_workspace      = true;
data.workspace_name      = result_basename;

pwld_solve_problem(data);


