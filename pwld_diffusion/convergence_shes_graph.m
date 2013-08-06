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
leg=['Uniform '];

% shestakov
geo = 'shes_quad_mms_1_L1_nc8_emb1_a0.1.mat';
for k=1:7
    k1=strfind(geo,'_emb');
    str = strcat(dir,'shes\',geo(1:k1+3),int2str(k),geo(k1+5:end))
    expression=sprintf('%s %s','load',str);eval(expression);
    ndof_(k)=norm_data(1);
    erro_(k)=norm_data(2);
end
plot(log10(ndof_),log10(erro_),'s-');
leg=[leg ; 'a = 0.10'];

geo = 'shes_quad_mms_1_L1_nc7_emb1_a0.15.mat';
for k=1:7
    k1=strfind(geo,'_emb');
    str = strcat(dir,'shes\',geo(1:k1+3),int2str(k),geo(k1+5:end))
    expression=sprintf('%s %s','load',str);eval(expression);
    ndof_(k)=norm_data(1);
    erro_(k)=norm_data(2);
end
plot(log10(ndof_),log10(erro_),'d-');
leg=[leg ; 'a = 0.15'];

geo = 'shes_quad_mms_1_L1_nc7_emb1_a0.2.mat';
for k=1:7
    k1=strfind(geo,'_emb');
    str = strcat(dir,'shes\',geo(1:k1+3),int2str(k),geo(k1+5:end))
    expression=sprintf('%s %s','load',str);eval(expression);
    ndof_(k)=norm_data(1);
    erro_(k)=norm_data(2);
end
plot(log10(ndof_),log10(erro_),'o-');
leg=[leg ; 'a = 0.20'];

geo = 'shes_quad_mms_1_L1_nc7_emb1_a0.25.mat';
for k=1:7
    k1=strfind(geo,'_emb');
    str = strcat(dir,'shes\',geo(1:k1+3),int2str(k),geo(k1+5:end))
    expression=sprintf('%s %s','load',str);eval(expression);
    ndof_(k)=norm_data(1);
    erro_(k)=norm_data(2);
end
plot(log10(ndof_),log10(erro_),'h-');
leg=[leg ; 'a = 0.25'];

geo = 'shes_quad_mms_1_L1_nc7_emb1_a0.3.mat';
for k=1:7
    k1=strfind(geo,'_emb');
    str = strcat(dir,'shes\',geo(1:k1+3),int2str(k),geo(k1+5:end))
    expression=sprintf('%s %s','load',str);eval(expression);
    ndof_(k)=norm_data(1);
    erro_(k)=norm_data(2);
end
plot(log10(ndof_),log10(erro_),'^-');
leg=[leg ; 'a = 0.30'];

geo = 'shes_quad_mms_1_L1_nc7_emb1_a0.35.mat';
for k=1:7
    k1=strfind(geo,'_emb');
    str = strcat(dir,'shes\',geo(1:k1+3),int2str(k),geo(k1+5:end))
    expression=sprintf('%s %s','load',str);eval(expression);
    ndof_(k)=norm_data(1);
    erro_(k)=norm_data(2);
end
plot(log10(ndof_),log10(erro_),'h-');
leg=[leg ; 'a = 0.35'];

geo = 'shes_quad_mms_1_L1_nc7_emb1_a0.45.mat';
for k=1:7
    k1=strfind(geo,'_emb');
    str = strcat(dir,'shes\',geo(1:k1+3),int2str(k),geo(k1+5:end))
    expression=sprintf('%s %s','load',str);eval(expression);
    ndof_(k)=norm_data(1);
    erro_(k)=norm_data(2);
end
plot(log10(ndof_),log10(erro_),'v-');
leg=[leg ; 'a = 0.45'];

% ref
slop=-1;
% x0=5;y0=-3.25;
x0=1.2;y0=-.5;
yref=@(x) slop*(x-x0)+y0;
x1=4.25;
plot([x0 x1],[yref(x0) yref(x1)],'--','LineWidth',2);
leg=[leg ; 'Slope =1'];


%add legend
legend(leg,'Location','Best')
xlabel('log(number of unknowns)','Fontsize',12)
ylabel('log(error)','Fontsize',12)
axis tight
axis([1 5.25 -3.75 -0.25])
print('-dpdf',strcat('results\convergence\cv_shes_quad','.pdf'));
print('-dpng',strcat('results\convergence\cv_shes_quad','.png'));


% slop=-1;x0=5;y0=-3.5;
% yref=@(x) slop*(x-x0)+y0;
% x1=1.5;
% plot([x0 x1],[yref(x0) yref(x1)],'--','LineWidth',2);
% ref
% slop=-1;x0=5;y0=-3.;
% yref=@(x) slop*(x-x0)+y0;
% x1=1.5;
% plot([x0 x1],[yref(x0) yref(x1)],'--','LineWidth',2);
% ref
% slop=-0.5;x0=5;y0=-2;
% yref=@(x) slop*(x-x0)+y0;
% x1=1.5;
% plot([x0 x1],[yref(x0) yref(x1)],'--','LineWidth',2);
% ref
% slop=-0.25;x0=5;y0=-1;
% yref=@(x) slop*(x-x0)+y0;
% x1=1.5;
% plot([x0 x1],[yref(x0) yref(x1)],'--','LineWidth',2);

