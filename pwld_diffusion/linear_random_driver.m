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

% august 1 below:
% geo = 'random_quad_mesh_L1_n16_a0.65.txt' ; % one small oscillation, goes away with use of <0 det
% geo = 'random_quad_mesh_L1_n16_a0.5.txt' ;  % good 
% geo = 'random_quad_mesh_L1_n16_a0.66.txt' ; % good with use of <0 det
% geo = 'z_mesh_quad_L1_n20_a0.05.txt' ;  % good
% geo = 'shestakov_quad_L1_nc4_a0.15.txt'; % good, some wiggles, does not improve with <0 det
% geo = 'shestakov_quad_L1_nc4_a0.25.txt'; % good
% geo = 'misha_quad_L1_n4.txt'; % not working
% geo = 'random_quad_mesh_L1_n3_a0.txt' ; % good, must use <0 det

% geo = 'random_poly_mesh_L1_n15_a0.95.txt'; % good
% geo = 'random_poly_mesh_L1_n15_a0.99.txt'; % good
% geo = 'shestakov_poly_mesh_L1_nc4_a0.15.txt'; % good, not very crazy for a mesh ...
geo = 'shestakov_poly_mesh_L1_nc4_a0.02.txt'; % good, not very crazy for a mesh ...
geo = 'z_mesh_poly_L1_n20_a0.05.txt';
% 
% ----- geo file -----
if(strcmp(geo,''))
    data.geofile = '';
else
    data.geofile = sprintf('%s%s','..\geom_codes\figs\',geo);
end

% ---- mms/linear choice
data.logi_mms = false;
data.mms_type = 2;
if ~data.logi_mms
    data.pbtype = 'linear';
end

% ----- refinement choices
data.max_ref_cycles = 1; % must be >=1
data.ref_threshold  = 0; % 0=uniform refinement

% ---- plotting choices
data.logi_plot           = true;
data.logi_plot_err_i     = false;
data.generate_vtk_output = true;

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
% --- put pieces together for filename    
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

% --- PWLD SOLVE finally
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
