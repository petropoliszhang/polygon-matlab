close all; clc;
disp('Running script that drives the pwld code');

global verbose
verbose = false;

geo = 'shestakov_quad_L1_nc8_emb1_a0.1.txt';

% choose geo file 
if(strcmp(geo,''))
    data.geofile = '';
else
    data.geofile = sprintf('%s%s','..\geom_codes\figs\',geo);
end

% data
data.logi_mms = true;
data.mms_type = 1;
if ~data.logi_mms
    data.pbtype = 'linear';
end

data.max_ref_cycles = 1; % must be >=1
data.ref_threshold  = 0; % 0=uniform refinement

data.logi_plot           = true;
data.logi_plot_err_i     = false;
data.generate_vtk_output = true;

% build result name
result_dir = 'results\';
% is quad or poly mesh
found_poly = strfind(geo,'poly');
found_quad = strfind(geo,'quad');
if( ( isempty(found_poly) &&  isempty(found_quad) )|| ...
    (~isempty(found_poly) && ~isempty(found_quad) ) )
    error('impossible mesh type');
elseif ( ~isempty(found_poly) && isempty(found_quad) )
    mesh_type='poly';
elseif ( isempty(found_poly) && ~isempty(found_quad) )
    mesh_type='quad';
else
    error('unknown mesh type');
end
% is rand, shes, or z-mesh ?
found_rand = strfind(geo,'random');
if(isempty(found_rand))
    found_shes = strfind(geo,'shestakov');
    if(isempty(found_shes))
        found_z = strfind(geo,'z_mesh');
        if(isempty(found_z))
            found_misha = strfind(geo,'misha');
            if(isempty(found_misha))
                error('unknown mesh');
            else
                result_basename = strcat(result_dir,'mish_',mesh_type,'_');
            end
        else
            result_basename = strcat(result_dir,'zzzz_',mesh_type,'_');
        end
    else
        result_basename = strcat(result_dir,'shes_',mesh_type,'_');
    end
else
    result_basename = strcat(result_dir,'rand_',mesh_type,'_');
end

if(data.logi_mms)
    result_basename = sprintf('%s%s%c',result_basename,'mms_',int2str(data.mms_type));
else
    result_basename = strcat(result_basename,data.pbtype);
end


% save common portion of the name
result_basename_emb = result_basename;

data.save_workspace      = true;

%%%%% portion to repeat
k1=strfind(geo,'_L');
k2=strfind(geo,'.txt');
gg = geo(k1+1:k2-1);
result_basename = strcat(result_basename_emb,'_',gg);
result_basename
data.vtk_basename        = result_basename;
data.workspace_name      = strcat(result_basename,'.mat');

pwld_solve_problem(data);

%--------------------------------------------
geo = 'shestakov_quad_L1_nc8_emb2_a0.1.txt';
data.geofile = sprintf('%s%s','..\geom_codes\figs\',geo);
%%%%% portion to repeat
k1=strfind(geo,'_L'); k2=strfind(geo,'.txt');
gg = geo(k1+1:k2-1);
result_basename = strcat(result_basename_emb,'_',gg); result_basename
data.result_basename        = result_basename;
data.vtk_basename        = result_basename;
data.workspace_name      = strcat(result_basename,'.mat');
pwld_solve_problem(data);
%--------------------------------------------
geo = 'shestakov_quad_L1_nc8_emb3_a0.1.txt';
data.geofile = sprintf('%s%s','..\geom_codes\figs\',geo);
%%%%% portion to repeat
k1=strfind(geo,'_L'); k2=strfind(geo,'.txt');
gg = geo(k1+1:k2-1);
result_basename = strcat(result_basename_emb,'_',gg); result_basename
data.result_basename        = result_basename;
data.vtk_basename        = result_basename;
data.workspace_name      = strcat(result_basename,'.mat');
pwld_solve_problem(data);
%--------------------------------------------
geo = 'shestakov_quad_L1_nc8_emb4_a0.1.txt';
data.geofile = sprintf('%s%s','..\geom_codes\figs\',geo);
%%%%% portion to repeat
k1=strfind(geo,'_L'); k2=strfind(geo,'.txt');
gg = geo(k1+1:k2-1);
result_basename = strcat(result_basename_emb,'_',gg); result_basename
data.result_basename        = result_basename;
data.vtk_basename        = result_basename;
data.workspace_name      = strcat(result_basename,'.mat');
pwld_solve_problem(data);
%--------------------------------------------
geo = 'shestakov_quad_L1_nc8_emb5_a0.1.txt';
data.geofile = sprintf('%s%s','..\geom_codes\figs\',geo);
%%%%% portion to repeat
k1=strfind(geo,'_L'); k2=strfind(geo,'.txt');
gg = geo(k1+1:k2-1);
result_basename = strcat(result_basename_emb,'_',gg); result_basename
data.result_basename        = result_basename;
data.vtk_basename        = result_basename;
data.workspace_name      = strcat(result_basename,'.mat');
pwld_solve_problem(data);
%--------------------------------------------
geo = 'shestakov_quad_L1_nc8_emb6_a0.1.txt';
data.geofile = sprintf('%s%s','..\geom_codes\figs\',geo);
%%%%% portion to repeat
k1=strfind(geo,'_L'); k2=strfind(geo,'.txt');
gg = geo(k1+1:k2-1);
result_basename = strcat(result_basename_emb,'_',gg); result_basename
data.result_basename        = result_basename;
data.vtk_basename        = result_basename;
data.workspace_name      = strcat(result_basename,'.mat');
pwld_solve_problem(data);
%--------------------------------------------
geo = 'shestakov_quad_L1_nc8_emb7_a0.1.txt';
data.geofile = sprintf('%s%s','..\geom_codes\figs\',geo);
%%%%% portion to repeat
k1=strfind(geo,'_L'); k2=strfind(geo,'.txt');
gg = geo(k1+1:k2-1);
result_basename = strcat(result_basename_emb,'_',gg); result_basename
data.result_basename        = result_basename;
data.vtk_basename        = result_basename;
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
