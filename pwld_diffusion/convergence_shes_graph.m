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
% % % smooth
% % geo = 'smoo_quad_mms_1_L1_nc7_emb1_a0.15.mat';
% % for k=1:7
% %     k1=strfind(geo,'_emb');
% %     str = strcat(dir,'smooth\',geo(1:k1+3),int2str(k),geo(k1+5:end));
% %     expression=sprintf('%s %s','load',str);eval(expression);
% %     ndof_(k)=norm_data(1);
% %     erro_(k)=norm_data(2);
% % end
% % plot(log10(ndof_),log10(erro_),'+-');
% shestakov
geo = 'shes_quad_mms_1_L1_nc8_emb1_a0.1.mat';
for k=1:7
    k1=strfind(geo,'_emb');
    str = strcat(dir,'shes\',geo(1:k1+3),int2str(k),geo(k1+5:end))
    expression=sprintf('%s %s','load',str);eval(expression);
    ndof_(k)=norm_data(1);
    erro_(k)=norm_data(2);
end
plot(log10(ndof_),log10(erro_),'+-');

geo = 'shes_quad_mms_1_L1_nc7_emb1_a0.15.mat';
for k=1:7
    k1=strfind(geo,'_emb');
    str = strcat(dir,'shes\',geo(1:k1+3),int2str(k),geo(k1+5:end))
    expression=sprintf('%s %s','load',str);eval(expression);
    ndof_(k)=norm_data(1);
    erro_(k)=norm_data(2);
end
plot(log10(ndof_),log10(erro_),'+-');

geo = 'shes_quad_mms_1_L1_nc7_emb1_a0.2.mat';
for k=1:7
    k1=strfind(geo,'_emb');
    str = strcat(dir,'shes\',geo(1:k1+3),int2str(k),geo(k1+5:end))
    expression=sprintf('%s %s','load',str);eval(expression);
    ndof_(k)=norm_data(1);
    erro_(k)=norm_data(2);
end
plot(log10(ndof_),log10(erro_),'+-');

geo = 'shes_quad_mms_1_L1_nc7_emb1_a0.25.mat';
for k=1:7
    k1=strfind(geo,'_emb');
    str = strcat(dir,'shes\',geo(1:k1+3),int2str(k),geo(k1+5:end))
    expression=sprintf('%s %s','load',str);eval(expression);
    ndof_(k)=norm_data(1);
    erro_(k)=norm_data(2);
end
plot(log10(ndof_),log10(erro_),'+-');

geo = 'shes_quad_mms_1_L1_nc7_emb1_a0.3.mat';
for k=1:7
    k1=strfind(geo,'_emb');
    str = strcat(dir,'shes\',geo(1:k1+3),int2str(k),geo(k1+5:end))
    expression=sprintf('%s %s','load',str);eval(expression);
    ndof_(k)=norm_data(1);
    erro_(k)=norm_data(2);
end
plot(log10(ndof_),log10(erro_),'+-');

geo = 'shes_quad_mms_1_L1_nc7_emb1_a0.35.mat';
for k=1:7
    k1=strfind(geo,'_emb');
    str = strcat(dir,'shes\',geo(1:k1+3),int2str(k),geo(k1+5:end))
    expression=sprintf('%s %s','load',str);eval(expression);
    ndof_(k)=norm_data(1);
    erro_(k)=norm_data(2);
end
plot(log10(ndof_),log10(erro_),'+-');

geo = 'shes_quad_mms_1_L1_nc7_emb1_a0.45.mat';
for k=1:7
    k1=strfind(geo,'_emb');
    str = strcat(dir,'shes\',geo(1:k1+3),int2str(k),geo(k1+5:end))
    expression=sprintf('%s %s','load',str);eval(expression);
    ndof_(k)=norm_data(1);
    erro_(k)=norm_data(2);
end
plot(log10(ndof_),log10(erro_),'+-');

leg=['unif';'0.10';'0.15';'.020';'0.25';'0.30';'0.35';'0.45'];
legend(leg)
% % % z-mesh
% % str = strcat(dir,'z-mesh\zzzz_quad_mms_1.mat')
% % expression=sprintf('%s %s','load',str);eval(expression);
% % plot(log10(norm_data(:,1)),log10(norm_data(:,2)),'+-');

% --- POLY ---
% ------------
dir =strcat(basedir,'poly\');
