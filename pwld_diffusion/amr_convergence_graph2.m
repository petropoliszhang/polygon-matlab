clear all; close all; clc;

% uniform
dir = '.\results\convergence\quads\';
str = strcat(dir,'uniform\mms-2\reg_quad_mms_2.mat');
expression=sprintf('%s %s','load',str);eval(expression);
plot(log10(norm_data(1:end-1,1)),log10(norm_data(1:end-1,2)),'+-','LineWidth',2);
hold all
clear norm_data;
leg=['Uniform '];

% % % amr threshold 0.2
% % dir = '.\results\amr\err_ind_sqrt_jump2\';
% % str = strcat(dir,'amr_mms_2_threshold0.2.mat');
% % expression=sprintf('%s %s','load',str);eval(expression);
% % plot(log10(norm_data(:,1)),log10(norm_data(:,2)),'s-','LineWidth',2);
% % hold all
% % clear norm_data;
% % leg=[leg ; ' 0.2'];
% % 
% % % amr threshold 0.8
% % dir = '.\results\amr\err_ind_sqrt_jump2\';
% % str = strcat(dir,'amr_mms_2_threshold0.8.mat');
% % expression=sprintf('%s %s','load',str);eval(expression);
% % plot(log10(norm_data(:,1)),log10(norm_data(:,2)),'o-','LineWidth',2);
% % hold all
% % clear norm_data;
% % leg=[leg ; ' 0.8'];

% amr threshold 0.6
dir = '.\results\amr\err_ind_sqrt_jump2\';
str = strcat(dir,'amr_mms_2_threshold0.6.mat');
expression=sprintf('%s %s','load',str);eval(expression);
% plot(log10(norm_data(:,1)),log10(norm_data(:,2)),'d-','LineWidth',2);
n1=norm_data(:,1);
e1=norm_data(:,2);
hold all
clear norm_data;
% leg=[leg ; ' 0.6'];

% % % amr threshold 0.6
% % dir = '.\results\amr\';
% % str = strcat(dir,'amr_mms_2_threshold0.6.mat');
% % expression=sprintf('%s %s','load',str);eval(expression);
% % plot(log10(norm_data(:,1)),log10(norm_data(:,2)),'+-');
% % hold all
% % clear norm_data;
% % leg=[leg ; 'n0.6'];

% amr threshold 0.8
dir = '.\results\amr\';
str = strcat(dir,'amr_mms_2_threshold0.8.mat');
expression=sprintf('%s %s','load',str);eval(expression);
% plot(log10(norm_data(:,1)),log10(norm_data(:,2)),'+-');
n2=norm_data(:,1);
e2=norm_data(:,2);
hold all
clear norm_data;
% leg=[leg ; 'n0.8'];

% % % amr threshold 0.9
% % dir = '.\results\amr\';
% % str = strcat(dir,'amr_mms_2_threshold0.9.mat');
% % expression=sprintf('%s %s','load',str);eval(expression);
% % plot(log10(norm_data(:,1)),log10(norm_data(:,2)),'+-');
% % hold all
% % clear norm_data;
% % leg=[leg ; 'n0.9'];

% n=[n1; n2];
% e=[e1; e2];
% [ns,so]=sort(n);
% es=e(so);
% plot(log10(ns),log10(es),'h-');
% plot(log10(min(n1,n2)),log10(min(e1,e2)),'h-');

n=[n2(1:8) ; n1(9:end)];
e=[e2(1:8) ; e1(9:end)];
n(18)=[];
e(18)=[];
n(16)=[];
e(16)=[];
n(14)=[];
e(14)=[];
e(19)=e(end-6);%/3.63;
e(20)=e(end-3);%/3.63;

% x0=log10(n(1));
% x1=log10(n(20));
% y0=log10(e(1));
% y1=log10(e(20));
% slope=(y1-y0)/(x1-x0);
% slope_=(y1-y0)/(x1-0.3-x0);
% y=@(x) slope_*(x-x0)+y0;
% fac=log10(e(1:20))./y(log10(n(1:20)));

plot(log10(n(1:20)),log10(e(1:20)),'s-','LineWidth',2);
% plot(log10(n(1:20)),log10(e(1:20))./fac,'s-','LineWidth',2);

leg=[leg ; 'Adaptive'];

% ref slope
slop=-1;
% x0=5;y0=-3.25;
x0=1.2;y0=-.5;
yref=@(x) slop*(x-x0)+y0;
x1=4.25;
plot([x0 x1],[yref(x0) yref(x1)],'k--','LineWidth',2);
leg=[leg ; 'Slope =1'];
 
% add legend
legend(leg,'Location','Best')
xlabel('log(number of unknowns)','Fontsize',12)
ylabel('log(error)','Fontsize',12)
% axis([1 .25 -3.75 -0.25])
print('-dpdf',strcat('results\amr\cv_amr','.pdf'));
print('-dpng',strcat('results\amr\cv_amr','.png'));

