function pwl_function

close all; clc

%            D
%   E                  C
%
%
%   A                  B
%
h=1.45;
A=[0 0]; B=[1 0]; C=[1 1]; D=[0.5 h]; E=[0 1];

% compute mid-point
M=(A+B+C+D+E)/5

% pick grid
n=151
xx=linspace(0,1,n);
yy=linspace(0,max([1 h]),n);

for i=1:length(xx)
    x=xx(i);
    for j=1:length(yy)
        y=yy(j);
        zA(i,j)=f_A(M,h,x,y);
        zB(i,j)=f_B(M,h,x,y);
        zC(i,j)=f_C(M,h,x,y);
        zE(i,j)=f_E(M,h,x,y);
        zD(i,j)=f_D(M,h,x,y);
        zM(i,j)=f_M(M,h,x,y);
    end
end
i=0;con=true; nc=69;

i=i+1;figure(i)
if(con),contour(xx,yy,zA',nc);else,surf(xx,yy,zA');end
xlabel('x');ylabel('y');
i=i+1;figure(i)
if(con),contour(xx,yy,zB',nc);else,surf(xx,yy,zB');end
xlabel('x');ylabel('y');
i=i+1;figure(i)
if(con),contour(xx,yy,zC',nc);else,surf(xx,yy,zC');end
xlabel('x');ylabel('y');
i=i+1;figure(i)
if(con),contour(xx,yy,zE',nc);else,surf(xx,yy,zE');end
xlabel('x');ylabel('y');
i=i+1;figure(i)
if(con),contour(xx,yy,zD',nc);else,surf(xx,yy,zD');end
xlabel('x');ylabel('y');

i=i+1;figure(i)
if(con),contour(xx,yy,zM',nc);else,surf(xx,yy,zM');end
xlabel('x');ylabel('y');

i=i+1;figure(i)
if(con),contour(xx,yy,(zA+zM/5)',nc);else,surf(xx,yy,(zA+zM/5)');end
xlabel('x');ylabel('y');
i=i+1;figure(i)
if(con),contour(xx,yy,(zD+zM/5)',nc);else,surf(xx,yy,(zD+zM/5)');end
xlabel('x');ylabel('y');
i=i+1;figure(i)
if(con),contour(xx,yy,(zE+zM/5)',nc);else,surf(xx,yy,(zE+zM/5)');end
xlabel('x');ylabel('y');

% i=i+1;figure(i)
% contour3(xx,yy,(zA+zM/5)',100)
% xlabel('x');ylabel('y');

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out=ff_A(x,y)
if(x<0|y<0|1-x-y<0)
    out=0;
else
    out=1-x-y;
end
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=ff_B(x,y)
if(x>1|y<0|-y+x<0)
    out=0;
else
    out=-y+x;
end
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=ff_C(h,x,y)
if(x>1|y-2*h*(1-x)<0|-y+2*(1-h)*x+2*h-1<0)
    out=0;
else
    out=1*(y-2*h*(1-x))+(-y+2*(1-h)*x+2*h-1)*0;
end
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=ff_E(h,x,y)
if(x<0|y-2*h*x<0|-y+2*(h-1)*x+1<0)
    out=0;
else
    out=1*(y-2*h*x)+(-y+2*(1-h)*x+2*h-1)*0;
end
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=ff_D(h,x,y)
if(h>=1)
    if(y<1|-y+2*(1-h)*x+2*h-1<0|-y+2*(h-1)*x+1<0)
        out=0;
    else
        out=1;%(y-1)/(h-1);
    end
else
    if(y>1|-y+2*(1-h)*x+2*h-1>0|-y+2*(h-1)*x+1>0)
%     if(-y+2*(h-1)*x+1<0)
        out=0;
    else
        out=(y-1)/(h-1);
    end
end
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=f_M(M,h,x,y)
xm=M(1);ym=M(2);
s=ym/xm;
if(y>0&y-s*x<0&y-s*xm/(1-xm)*(1-x)<0)
    out=y/ym;
elseif(x>0&y-s*x>0&y-s/ym*(ym-1)*x-1<0)
    out=x/xm;
elseif(x>0&y-s*xm/(1-xm)*(1-x)>0&y-(ym-1)/(xm-1)*(x-1)-1<0)
    out=(x-1)/(xm-1);
elseif(x>0.5 &-y+2*(1-h)*x+2*h-1>0&y-(ym-1)/(xm-1)*(x-1)-1>0)
    out=(-y+2*(1-h)*x+2*h-1)/(h-ym);
elseif(x<=0.5&-y+2*(h-1)*x+1>0&y-s/ym*(ym-1)*x-1>0)
    out=(-y+2*(h-1)*x+1)/(h-ym);
else
    out=0;
end
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=f_A(M,h,x,y)
xm=M(1);ym=M(2);
s=ym/xm;
if(y>0&y-s*x<0&y-s*xm/(1-xm)*(1-x)<0)
    out=-(y-s*xm/(1-xm)*(1-x))/s;
elseif(x>0&y-s*x>0&y-s/ym*(ym-1)*x-1<0)
    out=-(y-s/ym*(ym-1)*x-1);
else
    out=0;
end
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=f_B(M,h,x,y)
xm=M(1);ym=M(2);
s=ym/xm;
if(y>0&y-s*x<0&y-s*xm/(1-xm)*(1-x)<0)
    out=-(y-s*x)/s;
elseif(x>0&y-s*xm/(1-xm)*(1-x)>0&y-(ym-1)/(xm-1)*(x-1)-1<0)
    out=-(y-(ym-1)/(xm-1)*(x-1)-1);
else
    out=0;
end
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=f_C(M,h,x,y)
xm=M(1);ym=M(2);
s=ym/xm;
if(x>0&y-s*xm/(1-xm)*(1-x)>0&y-(ym-1)/(xm-1)*(x-1)-1<0)
    out=y-s*xm/(1-xm)*(1-x);
elseif(x>0.5 &-y+2*(1-h)*x+2*h-1>0&y-(ym-1)/(xm-1)*(x-1)-1>0)
    out=(x-xm)/(1-xm);
else
    out=0;
end
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=f_E(M,h,x,y)
xm=M(1);ym=M(2);
s=ym/xm;
if(x<=0.5&-y+2*(h-1)*x+1>0&y-s/ym*(ym-1)*x-1>0)
    out=(xm-x)/xm;
elseif(x>0&y-s*x>0&y-s/ym*(ym-1)*x-1<0)
    out=(y-s*x);
else
    out=0;
end
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=f_D(M,h,x,y)
xm=M(1);ym=M(2);
s=ym/xm;
if(x<=0.5&-y+2*(h-1)*x+1>0&y-s/ym*(ym-1)*x-1>0)
    out=(y-s/ym*(ym-1)*x-1)/(h-ym);
elseif(x>0.5 &-y+2*(1-h)*x+2*h-1>0&y-(ym-1)/(xm-1)*(x-1)-1>0)
    out=(y-(ym-1)/(xm-1)*(x-1)-1)/(h-ym);
else
    out=0;
end
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

