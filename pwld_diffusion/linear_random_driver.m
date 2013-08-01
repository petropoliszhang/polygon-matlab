close all; clc;
disp('Running script that drives the pwld code');

global verbose
verbose = false;

geo = 'random_quad_mesh_L100_n10_a0.2.txt'   ;
geo = 'random_quad_mesh_L100_n10_a0.33.txt'  ;
geo = 'random_quad_mesh_L100_n2_a0.1.txt'    ;
geo = 'random_quad_mesh_L100_n30_a0.1.txt'   ;
geo = 'random_quad_mesh_L100_n30_a0.33.txt'  ;
geo = 'random_quad_mesh_L100_n3_a0.25.txt'   ;
% geo = 'random_quad_mesh_L100_n50_a0.33.txt'  ;
   

if(strcmp(geo,''))
    data.geofile = '';
else
    data.geofile = sprintf('%s%s','..\geom_codes\figs\',geo);
end

data.logi_mms = false;
data.mms_type = 2;
if ~data.logi_mms
    data.pbtype = 'linear';
end

data.max_ref_cycles = 1; % must be >=1
data.ref_threshold  = 0; % 0=uniform refinement

data.logi_plot           = true;
data.logi_plot_err_i     = false;
data.generate_vtk_output = true;

result_basename = 'results\reg_quad_';
if(data.logi_mms)
    result_basename = sprintf('%s%s%c',result_basename,'mms_',int2str(data.mms_type));
else
    k1=strfind(geo,'_L');
    k2=strfind(geo,'.txt');
    gg = geo(k1+1:k2-1);
    result_basename = strcat(result_basename,data.pbtype,'_',gg);
end
result_basename
data.vtk_basename        = result_basename;

data.save_workspace      = true;
data.workspace_name      = strcat(result_basename,'.mat');

pwld_solve_problem(data);



% useless, a0 --> rectangular grids
% % % random_quad_mesh_L100_n10_a0.txt     
% % % random_quad_mesh_L100_n50_a0.txt     
% % % random_quad_mesh_L100_n30_a0.txt     
% % % random_quad_mesh_L1_n128_a0.txt      
% % % random_quad_mesh_L1_n16_a0.txt       
% % % random_quad_mesh_L1_n2_a0.txt        
% % % random_quad_mesh_L1_n32_a0.txt       
% % % random_quad_mesh_L1_n40_a0.txt       
% % % random_quad_mesh_L1_n4_a0.txt        
% % % random_quad_mesh_L1_n64_a0.txt       
% % % random_quad_mesh_L1_n8_a0.txt  
