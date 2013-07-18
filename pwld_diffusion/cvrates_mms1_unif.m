% mms-1 with uniform mesh refinement
clear all; close all; clc;

% mms 1
% ndof =    16, L2_error = +4.744057E-01 
% ndof =    64, L2_error = +2.676629E-01 
% ndof =   256, L2_error = +8.295995E-02 
% ndof =  1024, L2_error = +2.199696E-02 
% ndof =  4096, L2_error = +5.584659E-03 
% ndof = 16384, L2_error = +1.401626E-03 
% ndof = 65536, L2_error = +3.507498E-04  

ndof = [   16
           64
          256
         1024
         4096
        16384
        65536
		];
		
L2 = [	
+4.744057E-01
+2.676629E-01
+8.295995E-02
+2.199696E-02
+5.584659E-03
+1.401626E-03
+3.507498E-04
		];

norm_data = [ ndof L2 ]		

norm_data (1,:)=[];
norm_data (1,:)=[];
norm_data (1,:)=[];

figure(999); hold all;
plot(log10(norm_data(:,1)),log10(norm_data(:,2)),'+-')
nel_ = norm_data(:,1)/4;
h_ = 1./ sqrt(nel_);
plot(log10(1./h_),log10(norm_data(:,2)),'+-')



