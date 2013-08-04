close all; clc;
disp('Running script that drives the pwld code');

global verbose
verbose = false;

% geo = 'shestakov_quad_L1_nc8_emb1_a0.1.txt';
% % geo = 'smooth_quad_mesh_L100_nc6_emb1_a0.15.txt';
% geo = 'random_poly_mesh_L1_n4_a0.965.txt';
% geo = 'random_poly_mesh_L1_n2_a0.9.txt';
% geo = 'smooth_poly_mesh_L1_n2_a0.15.txt';
% geo = 'shestakov_poly_mesh_L1_nc1_a0.02.txt';
% geo = 'shestakov_poly_mesh_L1_nc1_a0.1.txt';
% geo = 'shestakov_poly_mesh_L1_nc1_a0.15.txt';
% geo = 'shestakov_poly_mesh_L1_nc1_a0.25.txt';

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

norm_data{1}=pwld_solve_problem(data);

for k=2:6
%     k1=strfind(geo,'_n')+1; 
    k1=strfind(geo,'_nc')+2; 
    k2=strfind(geo,'_a0'); 
    % new geo
% %     geo = strcat(geo(1:k1),int2str(2^(k)),geo(k2:end))
%    geo = strcat(geo(1:k1),int2str(k+1),geo(k2:end))
    geo = strcat(geo(1:k1),int2str(k),geo(k2:end))
    data.geofile = strcat('..\geom_codes\figs\',geo);
    k1=strfind(geo,'_L'); k2=strfind(geo,'.txt');
    gg = geo(k1+1:k2-1);
    result_basename = strcat(result_basename_emb,'_',gg); result_basename
    data.result_basename        = result_basename;
    data.vtk_basename        = result_basename;
    data.workspace_name      = strcat(result_basename,'.mat');
    norm_data{k}=pwld_solve_problem(data);
end

workspace_name = strcat(result_basename,'_final','.mat');
expression = sprintf('save %s',workspace_name');
eval(expression);

for k=1:length(norm_data)
    ndof(k)=norm_data{k}(1);
    erro(k)=norm_data{k}(2);
end

nel_ = ndof/4;
h_ = 1./ sqrt(nel_);
figure(30)
plot(log10(nel_),log10(erro),'+-')
figure(31)
plot(log10(1./h_),log10(erro),'+-')

