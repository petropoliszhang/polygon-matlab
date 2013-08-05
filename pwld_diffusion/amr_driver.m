close all; clc;
disp('Running script that drives pwld');

global verbose
verbose = false;

data.geofile = '';

data.logi_mms = true;
data.mms_type = 2;
if ~data.logi_mms
    data_pbtype = 'linear';
end

data.max_ref_cycles = 40;
data.ref_threshold  = 0.9; % 0=uniform refinement

data.logi_plot           = true;
data.logi_plot_err_i     = false;
data.generate_vtk_output = true;

if(data.logi_mms)
    result_basename = strcat('results\amr_mms_',int2str(data.mms_type),...
        '_threshold',num2str(data.ref_threshold,2));
else
    error('only mms for amr')
    result_basename = strcat('results\amr_',data_pbtype);
end
result_basename
data.vtk_basename        = result_basename;

data.save_workspace      = true;
data.workspace_name      = strcat(result_basename,'.mat');

pwld_solve_problem(data);


