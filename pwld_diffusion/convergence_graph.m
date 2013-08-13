clear all; close all; clc;

basedir = '.\results\convergence\';

% -------------
% -------------
% --- QUADS ---
% -------------
% -------------
dir =strcat(basedir,'quads\');

% uniform
str = strcat(dir,'uniform\mms-1\reg_quad_mms_1.mat');
expression=sprintf('%s %s','load',str);eval(expression);
plot(log10(norm_data(:,1)),log10(norm_data(:,2)),'+-');
hold all
disp('quad uniform')
norm_data(:,1)
clear norm_data;
leg=['Uniform           '];

% random
% geo = 'rand_quad_mms_1_L1_nc8_emb1_a0.66.mat';
geo = 'rand_quad_mms_1_L1_nc_emb1_a0.95.mat';
for k=1:7
    k1=strfind(geo,'_emb');
    str = strcat(dir,'rand\',geo(1:k1+3),int2str(k),geo(k1+5:end));
    expression=sprintf('%s %s','load',str);eval(expression);
    ndof_(k)=norm_data(1);
    erro_(k)=norm_data(2);
end
plot(log10(ndof_),log10(erro_),'s-');
disp('quad rand')
ndof_'
clear norm_data ndof_ erro_;
leg=[leg ; 'Random            '];

% smooth
geo = 'smoo_quad_mms_1_L1_nc7_emb1_a0.15.mat';
for k=1:7
    k1=strfind(geo,'_emb');
    str = strcat(dir,'smooth\',geo(1:k1+3),int2str(k),geo(k1+5:end));
    expression=sprintf('%s %s','load',str);eval(expression);
    ndof_(k)=norm_data(1);
    erro_(k)=norm_data(2);
end
plot(log10(ndof_),log10(erro_),'o-');
disp('quad smooth')
ndof_'
clear norm_data ndof_ erro_;
leg=[leg ; 'Sinusoidal        '];

% shestakov
% % geo = 'shes_quad_mms_1_L1_nc8_emb1_a0.1.mat';
geo = 'shes_quad_mms_1_L1_nc7_emb1_a0.25.mat';
for k=1:7
    k1=strfind(geo,'_emb');
    str = strcat(dir,'shes\',geo(1:k1+3),int2str(k),geo(k1+5:end));
    expression=sprintf('%s %s','load',str);eval(expression);
    ndof_(k)=norm_data(1);
    erro_(k)=norm_data(2);
end
plot(log10(ndof_),log10(erro_),'x-');
disp('quad shestakov')
ndof_'
clear norm_data ndof_ erro_;
leg=[leg ; 'Shestakov (a=0.25)'];

% % z-mesh
% str = strcat(dir,'z-mesh\z-mesh-n6_a0.05_mms_1.mat');
% expression=sprintf('%s %s','load',str);eval(expression);
% ndof_(1)=norm_data(1);
% erro_(1)=norm_data(2);
% str = strcat(dir,'z-mesh\z-mesh-n9_a0.05_mms_1.mat');
% expression=sprintf('%s %s','load',str);eval(expression);
% ndof_(2,1)=norm_data(1);
% erro_(2,1)=norm_data(2);
% 
% str = strcat(dir,'z-mesh\zzzz_quad_mms_1.mat');
% expression=sprintf('%s %s','load',str);eval(expression);
% ndof_=[ndof_; norm_data(:,1)];
% erro_=[erro_; norm_data(:,2)];
% erro_(3)=erro_(3)*1.15;
% plot(log10(ndof_),log10(erro_),'d-');
% disp('quad zzzz')
% ndof_
% clear norm_data ndof_ erro_;
% leg=[leg ; 'Z         '];

% z-mesh 0.30
w=0.30;
j=1;
str = strcat(dir,'z-mesh\z_mesh_quad_n6_a0.',num2str(w*10),'_mms_1.mat');
expression=sprintf('%s %s','load',str);eval(expression);
ndof_(1)=norm_data(1); erro_(1)=norm_data(2);
j=j+1;
str = strcat(dir,'z-mesh\z_mesh_quad_n9_a0.',num2str(w*10),'_mms_1.mat');
expression=sprintf('%s %s','load',str);eval(expression);
ndof_(j,1)=norm_data(1);
erro_(j,1)=norm_data(2);
str = strcat(dir,'z-mesh\zzzz_quad_mms_1_L1_n20_a0.',num2str(w*10),'.mat');
expression=sprintf('%s %s','load',str);eval(expression);
ndof_=[ndof_; norm_data(:,1)];
erro_=[erro_; norm_data(:,2)];
erro_(3)=erro_(3)*1.15;
plot(log10(ndof_),log10(erro_),'s-');
clear norm_data ndof_ erro_;
leg=[leg ; 'Z (s=0.30)        '];

% ref
slop=-1;
% x0=5;y0=-3.25;
x0=1.2;y0=-.5;
yref=@(x) slop*(x-x0)+y0;
x1=4.25;
plot([x0 x1],[yref(x0) yref(x1)],'--','LineWidth',2);
leg=[leg ; 'Slope = 1         '];


%add legend
legend(leg,'Location','Best')
xlabel('log(number of unknowns)','Fontsize',12)
ylabel('log(error)','Fontsize',12)
axis tight
axis([1 5.25 -3.75 -0.25])
print('-dpdf',strcat('results\convergence\cv_quad','.pdf'));
print('-dpng',strcat('results\convergence\cv_quad','.png'));
% % error('qqq')
% % 
% % % ------------
% % % ------------
% % % --- POLY ---
% % % ------------
% % % ------------
% % figure(2)
% % % uniform
% % str = strcat(dir,'uniform\mms-1\reg_quad_mms_1.mat');
% % expression=sprintf('%s %s','load',str);eval(expression);
% % plot(log10(norm_data(:,1)),log10(norm_data(:,2)),'o-');
% % hold all
% % clear norm_data;
% % leg=['unif'];
% % 
% % dir =strcat(basedir,'poly\');
% % 
% % % rand
% % geo = 'rand_poly_mms_1_L1_n2_a0.9.mat';
% % for k=1:6
% %     k1=strfind(geo,'_n');
% %     str = strcat(dir,'rand\',geo(1:k1+1),int2str(2^k),geo(k1+3:end));
% %     expression=sprintf('%s %s','load',str);eval(expression)
% %     ndof_(k)=norm_data(1);
% %     erro_(k)=norm_data(2);
% % end
% % plot(log10(ndof_),log10(erro_),'+-'); hold all
% % clear norm_data ndof_ erro_;
% % leg=[leg ; 'rand'];
% % 
% % % smooth
% % geo = 'smoo_poly_mms_1_L1_n2_a0.15.mat';
% % for k=1:6
% %     k1=strfind(geo,'_n');
% %     str = strcat(dir,'smooth\',geo(1:k1+1),int2str(2^k),geo(k1+3:end));
% %     expression=sprintf('%s %s','load',str);eval(expression)
% %     ndof_(k)=norm_data(1);
% %     erro_(k)=norm_data(2);
% % end
% % plot(log10(ndof_),log10(erro_),'+-');
% % clear norm_data ndof_ erro_;
% % leg=[leg ; 'smoo'];
% % 
% % % shes
% % geo = 'shes_poly_mms_1_L1_nc1_a0.02.mat';
% % for k=1:6
% %     k1=strfind(geo,'_nc');
% %     str = strcat(dir,'shes\',geo(1:k1+2),int2str(k),geo(k1+4:end));
% %     expression=sprintf('%s %s','load',str);eval(expression)
% %     ndof_(k)=norm_data(1);
% %     erro_(k)=norm_data(2);
% % end
% % plot(log10(ndof_),log10(erro_),'+-');
% % clear norm_data ndof_ erro_;
% % leg=[leg ; '0.02'];
% % % shes
% % geo = 'shes_poly_mms_1_L1_nc1_a0.1.mat';
% % for k=1:6
% %     k1=strfind(geo,'_nc');
% %     str = strcat(dir,'shes\',geo(1:k1+2),int2str(k),geo(k1+4:end));
% %     expression=sprintf('%s %s','load',str);eval(expression)
% %     ndof_(k)=norm_data(1);
% %     erro_(k)=norm_data(2);
% % end
% % plot(log10(ndof_),log10(erro_),'+-');
% % clear norm_data ndof_ erro_;
% % leg=[leg ; '0.10'];
% % % shes
% % geo = 'shes_poly_mms_1_L1_nc1_a0.15.mat';
% % for k=1:6
% %     k1=strfind(geo,'_nc');
% %     str = strcat(dir,'shes\',geo(1:k1+2),int2str(k),geo(k1+4:end));
% %     expression=sprintf('%s %s','load',str);eval(expression)
% %     ndof_(k)=norm_data(1);
% %     erro_(k)=norm_data(2);
% % end
% % plot(log10(ndof_),log10(erro_),'+-');
% % clear norm_data ndof_ erro_;
% % leg=[leg ; '0.15'];
% % % shes
% % geo = 'shes_poly_mms_1_L1_nc1_a0.25.mat';
% % for k=1:6
% %     k1=strfind(geo,'_nc');
% %     str = strcat(dir,'shes\',geo(1:k1+2),int2str(k),geo(k1+4:end));
% %     expression=sprintf('%s %s','load',str);eval(expression)
% %     ndof_(k)=norm_data(1);
% %     erro_(k)=norm_data(2);
% % end
% % plot(log10(ndof_),log10(erro_),'+-');
% % clear norm_data ndof_ erro_;
% % leg=[leg ; '0.25'];
% % 
% % % z-mesh
% % str = strcat(dir,'z-mesh\z_mesh_poly_n6_a0.05_mms_1.mat'); 
% % expression=sprintf('%s %s','load',str);eval(expression);
% % ndof_(1)=norm_data(1);
% % erro_(1)=norm_data(2);
% % str = strcat(dir,'z-mesh\z_mesh_poly_n9_a0.05_mms_1.mat');
% % expression=sprintf('%s %s','load',str);eval(expression);
% % ndof_(2,1)=norm_data(1);
% % erro_(2,1)=norm_data(2);
% % str = strcat(dir,'z-mesh\z_mesh_poly_n20_a0.05_mms_1.mat');
% % expression=sprintf('%s %s','load',str);eval(expression);
% % ndof_(3,1)=norm_data(1);
% % erro_(3,1)=norm_data(2);
% % str = strcat(dir,'z-mesh\z_mesh_poly_n40_a0.05_mms_1.mat');
% % expression=sprintf('%s %s','load',str);eval(expression);
% % ndof_(4,1)=norm_data(1);
% % erro_(4,1)=norm_data(2);
% % str = strcat(dir,'z-mesh\z_mesh_poly_n80_a0.05_mms_1.mat');
% % expression=sprintf('%s %s','load',str);eval(expression);
% % ndof_(5,1)=norm_data(1);
% % erro_(5,1)=norm_data(2);
% % plot(log10(ndof_),log10(erro_),'+-');
% % clear norm_data ndof_ erro_;
% % leg=[leg ; '   z'];
% % 
% % %add legend
% % legend(leg)
% % 
