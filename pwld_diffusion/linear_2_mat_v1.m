% % analytical solutions
% % %%%%%%%%%%%%%%%%%%%%
% % 
% % discontinuous D:
% % 
% % u(0) -2D1 du/dx(0) =F
% % u(L/2)=v(L/2)
% % D1 du/dx(L/2) = D2 dv/dx(L/2)
% % v(L) +2D2 dv/dx(L) =0
% % 
% % u(x)=ax+b
% % v(x)=c(L-x)+d
% % 
% % b -2D1 a = F
% % aL/2 + b = cL/2 + d
% % D1 a = -D2 c
% % d -2D2 c = 0

figure(44)

F=9;
L=1;
D1=100;
D2=1;

a=-4*F/(L/2*(1+D1/D2)+4*D1)
b=4*F+2*D1*a
c=-D1/D2*a
d=-2*D1*a

y1=@(x) a*x+b;
y2=@(x) c*(L-x)+d;

x1=linspace(0,L/2);
x2=linspace(L/2,L);

y1_=y1(x1);
y2_=y2(x2);

plot([x1 x2],[y1_ y2_])


