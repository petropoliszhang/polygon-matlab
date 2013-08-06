clear all; close all; clc;

basedir = '.\results\convergence\';

% ------------
% ------------
% --- POLY ---
% ------------
% ------------
figure(1)
% uniform
% % str = strcat(basedir,'quads\uniform\mms-1\reg_quad_mms_1.mat');
% % expression=sprintf('%s %s','load',str);eval(expression);
% % plot(log10(norm_data(1:end-1,1)),log10(norm_data(1:end-1,2)),'+-');
% % hold all
% % disp('quad uniform')
% % norm_data(:,1)
% % clear norm_data;
% % leg=['Uniform       '];

dir =strcat(basedir,'poly\');

hold all

% random
% geo = 'rand_quad_mms_1_L1_nc8_emb1_a0.66.mat';
geo = 'rand_poly_mms_1_L1_n2_a0.9.mat';
for k=1:6
    k1=strfind(geo,'_n');
    str = strcat(dir,'rand\',geo(1:k1+1),int2str(2^k),geo(k1+3:end));
    expression=sprintf('%s %s','load',str);eval(expression)
    ndof_(k)=norm_data(1);
    erro_(k)=norm_data(2);
end
plot(log10(ndof_),log10(erro_),'d-');
disp('poly rand')
ndof_'
clear norm_data ndof_ erro_;
% leg=[leg ; 'Random        '];
leg=['Random        '];

% smooth
geo = 'smoo_poly_mms_1_L1_n2_a0.15.mat';
for k=1:6
    k1=strfind(geo,'_n');
    str = strcat(dir,'smooth\',geo(1:k1+1),int2str(2^k),geo(k1+3:end));
    expression=sprintf('%s %s','load',str);eval(expression)
    ndof_(k)=norm_data(1);
    erro_(k)=norm_data(2);
end
plot(log10(ndof_),log10(erro_),'o-');
disp('poly smooth')
ndof_'
clear norm_data ndof_ erro_;
leg=[leg ; 'Sinusoidal    '];

% shestakov
% % geo = 'shes_poly_mms_1_L1_nc1_a0.02.mat';
% % for k=1:6
% %     k1=strfind(geo,'_nc');
% %     str = strcat(dir,'shes\',geo(1:k1+2),int2str(k),geo(k1+4:end));
% %     expression=sprintf('%s %s','load',str);eval(expression)
% %     ndof_(k)=norm_data(1);
% %     erro_(k)=norm_data(2);
% % end
% % plot(log10(ndof_),log10(erro_),'x-');
% % disp('poly shestakov')
% % ndof_'
% % clear norm_data ndof_ erro_;
% % leg=[leg ; 'Shestakov 0.02'];
% % 
% % geo = 'shes_poly_mms_1_L1_nc1_a0.1.mat';
% % for k=1:6
% %     k1=strfind(geo,'_nc');
% %     str = strcat(dir,'shes\',geo(1:k1+2),int2str(k),geo(k1+4:end));
% %     expression=sprintf('%s %s','load',str);eval(expression)
% %     ndof_(k)=norm_data(1);
% %     erro_(k)=norm_data(2);
% % end
% % plot(log10(ndof_),log10(erro_),'x-');
% % disp('poly shestakov')
% % ndof_'
% % clear norm_data ndof_ erro_;
% % leg=[leg ; 'Shestakov 0.10'];

geo = 'shes_poly_mms_1_L1_nc1_a0.15.mat';
for k=1:6
    k1=strfind(geo,'_nc');
    str = strcat(dir,'shes\',geo(1:k1+2),int2str(k),geo(k1+4:end));
    expression=sprintf('%s %s','load',str);eval(expression)
    ndof_(k)=norm_data(1);
    erro_(k)=norm_data(2);
end
plot(log10(ndof_),log10(erro_),'x-');
disp('poly shestakov')
ndof_'
clear norm_data ndof_ erro_;
leg=[leg ; 'Shestakov 0.15'];

geo = 'shes_poly_mms_1_L1_nc1_a0.25.mat';
for k=1:6
    k1=strfind(geo,'_nc');
    str = strcat(dir,'shes\',geo(1:k1+2),int2str(k),geo(k1+4:end));
    expression=sprintf('%s %s','load',str);eval(expression)
    ndof_(k)=norm_data(1);
    erro_(k)=norm_data(2);
end
plot(log10(ndof_),log10(erro_),'v-');
disp('poly shestakov')
ndof_'
clear norm_data ndof_ erro_;
leg=[leg ; 'Shestakov 0.25'];



% z-mesh 0.05
str = strcat(dir,'z-mesh\z_mesh_poly_n6_a0.05_mms_1.mat'); 
expression=sprintf('%s %s','load',str);eval(expression);
ndof_(1)=norm_data(1);
erro_(1)=norm_data(2);
str = strcat(dir,'z-mesh\z_mesh_poly_n9_a0.05_mms_1.mat');
expression=sprintf('%s %s','load',str);eval(expression);
ndof_(2,1)=norm_data(1);
erro_(2,1)=norm_data(2);
str = strcat(dir,'z-mesh\z_mesh_poly_n20_a0.05_mms_1.mat');
expression=sprintf('%s %s','load',str);eval(expression);
ndof_(3,1)=norm_data(1);
erro_(3,1)=norm_data(2);
str = strcat(dir,'z-mesh\z_mesh_poly_n40_a0.05_mms_1.mat');
expression=sprintf('%s %s','load',str);eval(expression);
ndof_(4,1)=norm_data(1);
erro_(4,1)=norm_data(2);
str = strcat(dir,'z-mesh\z_mesh_poly_n80_a0.05_mms_1.mat');
expression=sprintf('%s %s','load',str);eval(expression);
ndof_(5,1)=norm_data(1);
erro_(5,1)=norm_data(2);
plot(log10(ndof_),log10(erro_),'s-');
clear norm_data ndof_ erro_;
leg=[leg ; 'Z             '];

% ref
slop=-1;
x0=1.52;y0=-.5;
x0=1.8;y0=-.5;
yref=@(x) slop*(x-x0)+y0;
x1=4.4;
plot([x0 x1],[yref(x0) yref(x1)],'k--','LineWidth',2);
leg=[leg ; 'Slope = 1     '];


%add legend
legend(leg,'Location','SouthWest')
xlabel('log(number of unknowns)','Fontsize',12)
ylabel('log(error)','Fontsize',12)
axis tight
axis([1.5 4.65 -3.25 -0.25])
print('-dpdf',strcat('results\convergence\cv_poly','.pdf'));
print('-dpng',strcat('results\convergence\cv_poly','.png'));

