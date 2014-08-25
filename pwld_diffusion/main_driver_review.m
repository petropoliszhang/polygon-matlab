close all; clc; clear all;
disp('Running script that drives pwld');

global verbose
verbose = false;

% ----- geo file -----
% geo = 'z_mesh_quad_L1_n20_a0.05.txt';
geo = 'z_mesh_quad_L1_n20_a0.25.txt';
geo = 'z_mesh_poly_L1_n20_a0.05.txt';
geo ='z_mesh_quad_L1_n20_a0.05.txt';

% geo = 'z_mesh_quad_L1_n20_a0.2.txt';
% % geo = 'z_mesh_quad_L1_n20_a0.4.txt';
% % geo = 'z_mesh_quad_L1_n20_a0.25.txt';
% geo = 'z_mesh_quad_L1_n20_a0.35.txt';
% geo = 'random_quad_mesh_L1_n2_a0.txt'
geo = 'shestakov_poly_mesh_L1_nc5_a0.1.txt'
geo ='shestakov_quad_nc5_a0.1.txt';
% geo=''
if(strcmp(geo,''))
    data.geofile = '';
else
    data.geofile = sprintf('%s%s','..\geom_codes\figs\',geo);
end

% ---- mms/linear choice
data.logi_mms = true;
data.mms_type = 3;
if ~data.logi_mms
    data.pbtype = 'linear'
end

% ----- refinement choices
data.max_ref_cycles = 1;
data.ref_threshold  = 0; % 0=uniform refinement
data.logi_amr_output = false;

% ---- plotting choices
data.logi_plot           = true;
data.logi_plot_err_i     = false;
data.generate_vtk_output = false;

data.save_workspace      = false;

if data.generate_vtk_output
    % ---- create result filename
    result_dir = 'results\';
    % ----  is quad or poly mesh ?
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
    % ---- is rand, shes, z-mesh, or misha ?
    found_rand = strfind(geo,'random');
    if(isempty(found_rand))
        found_shes = strfind(geo,'shestakov');
        if(isempty(found_shes))
            found_z = strfind(geo,'z_mesh');
            if(isempty(found_z))
                found_smooth = strfind(geo,'smooth');
                if(isempty(found_smooth))
                    error('unknown mesh');
                else
                    result_basename = strcat(result_dir,'smoo_',mesh_type,'_');
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
    % --- put pieces together for filename
    if(data.logi_mms)
        result_basename = sprintf('%s%s%c',result_basename,'mms_',int2str(data.mms_type));
    else
        result_basename = strcat(result_basename,data.pbtype);
        
    end
    %%%%% portion to repeat
    k1=strfind(geo,'_L');
    k2=strfind(geo,'.txt');
    gg = geo(k1+1:k2-1);
    result_basename = strcat(result_basename,'_',gg);
    
else
    result_basename='';
end
result_basename
data.vtk_basename        = result_basename;
data.workspace_name      = strcat(result_basename,'.mat');


% --- PWLD SOLVE finally
pwld_solve_problem(data);


