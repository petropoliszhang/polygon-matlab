% mms2 with uniform mesh refinement
clear all; close all; clc;

% mms 2
% ndof =    16, L2_error = +5.499127E-01 
% ndof =    64, L2_error = +4.270107E-01 
% ndof =   256, L2_error = +1.982779E-01
% ndof =  1024, L2_error = +5.257954E-02 
% ndof =  4096, L2_error = +1.393577E-02 
% ndof = 16384, L2_error = +3.540492E-03 
% ndof = 65536, L2_error = +8.887954E-04 

ndof = [   16
           64
          256
         1024
         4096
        16384
		65536
		];
		
L2 = [	+5.499127E-01 	
        +4.270107E-01 
        +1.982779E-01
        +5.257954E-02 
        +1.393577E-02 
        +3.540492E-03
        +8.887954E-04		
		];

norm_data = [ ndof L2 ]		

norm_data (1,:)=[];
norm_data (1,:)=[];

figure(999); hold all;
plot(log10(norm_data(:,1)),log10(norm_data(:,2)),'+-')
nel_ = norm_data(:,1)/4;
h_ = 1./ sqrt(nel_);
plot(log10(1./h_),log10(norm_data(:,2)),'+-')



