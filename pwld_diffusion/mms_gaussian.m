function mms_gaussian
clear all; clc; close all
% % syms x y Lx Ly x0 y0 c_diff sigma_a varia
% % 
% % U = 100*x*(Lx-x)*y*(Ly-y)*exp(-((x-x0)^2+(y-y0)^2)/varia)/(Lx*Ly)^2
% % 
% % forcing_fn = -1.0*c_diff*( diff(diff(U,x),x) +  diff(diff(U,y),y) ) + sigma_a *U
% % 
% % codeU = matlabFunction(U)
% % 
% % coderhs= matlabFunction(forcing_fn)


Lx=10;
Ly=Lx;
x0=Lx/4;y0=x0;
varia=Lx^2/100;

exact=@(x,y) 100*x.*(Lx-x).*y.*(Ly-y).*exp(-((x-x0).^2+(y-y0).^2)/varia)/(Lx*Ly)^2;
xx=linspace(0,Lx);
yy=linspace(0,Ly);
[uu,vv]=meshgrid(xx,yy);
zz=exact(uu,vv);
figure(99)
surf(uu,vv,zz)

return
end