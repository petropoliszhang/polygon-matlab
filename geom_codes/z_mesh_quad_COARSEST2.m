function z_mesh_quad()
% Z-mesh generated for quadrilaterals
clear all; close all; clc;

logi_save_plots = false;
logi_write_file = true;

L=1;

fraction=0.4;
if(fraction<eps || fraction > 1-eps)
    error('wrong fraction');
end

slope1 = (1-2*fraction)/0.15;
slope2 = (2*fraction-1)/0.30;

y1=@(x) slope1*(x-0.20*L)+fraction*L;
y2=@(x) slope2*(x-0.35*L)+(1-fraction)*L;
y3=@(x) slope1*(x-0.65*L)+fraction*L;

%%%%%%
% nx=5;
% nsub=[1 2 3 4 5] + 1;
% xi=[0 0.2 0.35 0.65 0.8 1.0];
% ny=2;
%%%%%%
nx=6;
nsub=[1 2 4 5 6] + 1;
xi=[0 0.2 0.35 0.50 0.65 0.8 1.0];
ny=4;
%%%%%%
% nx=9;
% nsub=[2 3 6 7 9] + 1;
% xi=[0 0.1 0.2 0.35 0.45 0.55 0.65 0.8 0.9 1.0];
% ny=10;


x=zeros(nx+1,ny+1);
y=zeros(nx+1,ny+1);

for i=1:nx+1
    x(i,:) = xi(i);
end
% five different y-spacing, depending on the zone
for i=1:nx+1
    if(i<=nsub(1))
        yy=L*fraction;
    elseif(i<=nsub(2))
        yy=y1(x(i));
    elseif(i<=nsub(3))
        yy=y2(x(i));
    elseif(i<=nsub(4))
        yy=y3(x(i));
    elseif(i<=nsub(5))
        yy=L*(1-fraction);
    else
        error('not in nsub!!!');
    end
    yy
    aux1 = linspace(0,yy,ny/2+1);
    aux2 = linspace(yy,L,ny/2+1); aux2(1)=[];
    y(i,:) = [aux1 aux2];
end

figure(1)
surf(x,y,ones(nx+1,ny+1)*0,'FaceColor','white')
view(0,90)
output_file=strcat('.\figs\z_mesh_quad_L',int2str(L),'_n',int2str(nx),'_a',num2str(fraction,3))
if(logi_save_plots)
    print('-dpdf',strcat(output_file1,'.pdf'));
    print('-dpng',strcat(output_file1,'.png'));
    saveas(gcf,strcat(output_file1,'.fig'),'fig');
end
%
if(~logi_write_file), return; end;
%
%---------------------------------------------
% save txt file
matID=1;
srcID=1;

% date and time;
[yr, mo, da, hr, mi, s] = datevec(now);

%%%%%%%%%%%%%%%%%%%%%
output_file1=strcat(output_file,'.txt')
fid=fopen(output_file1,'w');

fprintf(fid,'# Date: %d/%d/%d   Time: %d:%d\n', mo, da, yr, hr, mi);

fprintf(fid,'# dimensions \n');
fprintf(fid,'%g %g \n',L,L);

ncells = nx*ny;

fprintf(fid,'# connectivity \n');
fprintf(fid,'%d\n',ncells);
for iel=1:ncells
    skip = 4*(iel-1);
    i1 = skip + 1;
    i2 = skip + 2;
    i3 = skip + 3;
    i4 = skip + 4;
    fprintf(fid,'%d %d %d %d %d %d %d \n',4,i1,i2,i3,i4,matID,srcID);
end

fprintf(fid,'# DG vertices (counter-clockwise) \n');
fprintf(fid,'%d\n',4*ncells);
for i=1:nx
    for j=1:ny
        fprintf(fid,' %g %g \n %g %g \n %g %g \n %g %g \n',x(i,j)    ,y(i,j)    ,...
            x(i+1,j)  ,y(i+1,j)  ,...
            x(i+1,j+1),y(i+1,j+1),...
            x(i,j+1)  ,y(i,j+1)  );
    end
end

fprintf(fid,'# grid vertices (counter-clockwise) \n');
fprintf(fid,'%d\n',(nx+1)*(ny+1));
for i=1:nx+1
    for j=1:ny+1
        fprintf(fid,'%g %g \n',x(i,j),y(i,j) );
    end
end

fclose(fid);
return
end
%%%%%%%%%%%%%%%%%%%%%
