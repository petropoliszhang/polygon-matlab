clear all; close all; clc;

basedir = '.\results\convergence\';

% --- QUADS ---
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
leg=['Uniform'];


% z-mesh 0.05
w=0.05;
j=1;
str = strcat(dir,'z-mesh\z-mesh-n6_a0.',num2str(w*100,'%2.2d'),'_mms_1.mat');
expression=sprintf('%s %s','load',str);eval(expression);
ndof_(1)=norm_data(1); erro_(1)=norm_data(2);
j=j+1;
str = strcat(dir,'z-mesh\z-mesh-n9_a0.',num2str(w*100,'%2.2d'),'_mms_1.mat');
expression=sprintf('%s %s','load',str);eval(expression);
ndof_(j,1)=norm_data(1);
erro_(j,1)=norm_data(2);
str = strcat(dir,'z-mesh\zzzz_quad_mms_1_a0.',num2str(w*100,'%2.2d'),'.mat');
expression=sprintf('%s %s','load',str);eval(expression);
ndof_=[ndof_; norm_data(:,1)];
erro_=[erro_; norm_data(:,2)];
erro_(3)=erro_(3)*1.15;
plot(log10(ndof_),log10(erro_),'d-');
clear norm_data ndof_ erro_;
leg=[leg ; strcat('Z, 0.',num2str(w*100,'%2.2d'))];

% z-mesh 0.10
w=0.10;
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
plot(log10(ndof_),log10(erro_),'o-');
clear norm_data ndof_ erro_;
leg=[leg ; strcat('Z, 0.',num2str(w*100,'%2.2d'),'   ')];

% z-mesh 0.15
w=0.15;
j=1;
str = strcat(dir,'z-mesh\z_mesh_quad_n6_a0.',num2str(w*100),'_mms_1.mat');
expression=sprintf('%s %s','load',str);eval(expression);
ndof_(1)=norm_data(1); erro_(1)=norm_data(2);
j=j+1;
str = strcat(dir,'z-mesh\z_mesh_quad_n9_a0.',num2str(w*100),'_mms_1.mat');
expression=sprintf('%s %s','load',str);eval(expression);
ndof_(j,1)=norm_data(1);
erro_(j,1)=norm_data(2);
str = strcat(dir,'z-mesh\zzzz_quad_mms_1_L1_n20_a0.',num2str(w*100,'%2.2d'),'.mat');
expression=sprintf('%s %s','load',str);eval(expression);
ndof_=[ndof_; norm_data(:,1)];
erro_=[erro_; norm_data(:,2)];
erro_(3)=erro_(3)*1.15;
plot(log10(ndof_),log10(erro_),'s-');
clear norm_data ndof_ erro_;
leg=[leg ; strcat('Z, 0.',num2str(w*100,'%2.2d'),'   ')];

% z-mesh 0.20
w=0.20;
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
leg=[leg ; strcat('Z, 0.',num2str(w*100,'%2.2d'),'   ')];

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
leg=[leg ; strcat('Z, 0.',num2str(w*100,'%2.2d'),'   ')];

% z-mesh 0.40
w=0.40;
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
leg=[leg ; strcat('Z, 0.',num2str(w*100,'%2.2d'),'   ')];


% ref slope
slop=-1;
% x0=5;y0=-3.25;
x0=1.2;y0=-.5;
yref=@(x) slop*(x-x0)+y0;
x1=4.25;
plot([x0 x1],[yref(x0) yref(x1)],'k--','LineWidth',2);
leg=[leg ; 'Slope=1'];

slop=-.80;
x0=3.5;y0=-0.5;
yref=@(x) slop*(x-x0)+y0;
x1=5.1;
plot([x0 x1],[yref(x0) yref(x1)],'k--','LineWidth',2);

% slop=-1;
% x0=3;y0=-0.5;
% yref=@(x) slop*(x-x0)+y0;
% x1=5.1;
% plot([x0 x1],[yref(x0) yref(x1)],'k--','LineWidth',2);

%add legend
legend(leg,'Location','Best')
xlabel('log(number of unknowns)','Fontsize',12)
ylabel('log(error)','Fontsize',12)
axis([1 5.25 -3.75 -0.25])
print('-dpdf',strcat('results\convergence\cv_Z_quad','.pdf'));
print('-dpng',strcat('results\convergence\cv_Z_quad','.png'));

