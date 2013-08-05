clear all; close all; clc;

% uniform
dir = '.\results\convergence\quads\';
str = strcat(dir,'uniform\mms-2\reg_quad_mms_2.mat');
expression=sprintf('%s %s','load',str);eval(expression);
plot(log10(norm_data(:,1)),log10(norm_data(:,2)),'+-');
hold all
clear norm_data;
leg=['unif'];

% amr threshold 0.2
dir = '.\results\amr\err_ind_sqrt_jump2\';
str = strcat(dir,'amr_mms_2_threshold0.2.mat');
expression=sprintf('%s %s','load',str);eval(expression);
plot(log10(norm_data(:,1)),log10(norm_data(:,2)),'+-');
hold all
clear norm_data;
leg=[leg ; ' 0.2'];

% amr threshold 0.8
dir = '.\results\amr\err_ind_sqrt_jump2\';
str = strcat(dir,'amr_mms_2_threshold0.8.mat');
expression=sprintf('%s %s','load',str);eval(expression);
plot(log10(norm_data(:,1)),log10(norm_data(:,2)),'+-');
hold all
clear norm_data;
leg=[leg ; ' 0.8'];

% amr threshold 0.6
dir = '.\results\amr\err_ind_sqrt_jump2\';
str = strcat(dir,'amr_mms_2_threshold0.6.mat');
expression=sprintf('%s %s','load',str);eval(expression);
plot(log10(norm_data(:,1)),log10(norm_data(:,2)),'+-');
hold all
clear norm_data;
leg=[leg ; ' 0.6'];

% amr threshold 0.6
dir = '.\results\amr\';
str = strcat(dir,'amr_mms_2_threshold0.6.mat');
expression=sprintf('%s %s','load',str);eval(expression);
plot(log10(norm_data(:,1)),log10(norm_data(:,2)),'+-');
hold all
clear norm_data;
leg=[leg ; 'n0.6'];

% amr threshold 0.8
dir = '.\results\amr\';
str = strcat(dir,'amr_mms_2_threshold0.8.mat');
expression=sprintf('%s %s','load',str);eval(expression);
plot(log10(norm_data(:,1)),log10(norm_data(:,2)),'+-');
hold all
clear norm_data;
leg=[leg ; 'n0.8'];

% amr threshold 0.9
dir = '.\results\amr\';
str = strcat(dir,'amr_mms_2_threshold0.9.mat');
expression=sprintf('%s %s','load',str);eval(expression);
plot(log10(norm_data(:,1)),log10(norm_data(:,2)),'+-');
hold all
clear norm_data;
leg=[leg ; 'n0.9'];


%add legend
legend(leg)



