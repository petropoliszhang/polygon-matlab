clear all; close all; clc
% random mesh generated for quadrilaterals

L=1;
n=10;
h=L/n;
xi=linspace(0,L,n+1);
eta=xi;
fraction=0.53;

for i=1:n+1
    ax=1;   if(i==1|i==n+1),ax=0;end
    for j=1:n+1
        ay=1;   if(j==1|j==n+1),ay=0;end
        x(i,j) = xi(i) + (2*rand(1,1)-1)*h*fraction*ax;
        y(i,j) = eta(j)+ (2*rand(1,1)-1)*h*fraction*ay;
    end
end

figure(1)
surf(x,y,ones(n+1,n+1))
view(0,90)
output_file1=strcat('.\figs\random_quad_mesh_L',int2str(L),'_n',int2str(n),'_a',num2str(fraction,3));
print('-dpdf',strcat(output_file1,'.pdf'));
print('-dpng',strcat(output_file1,'.png'));
saveas(gcf,strcat(output_file1,'.fig'),'fig');

% keep x=L/2 a straight line
if(mod(n,2)==0)
    xx=x;
    xx(n/2+1,:)=L/2;
    figure(2)
    surf(xx,y,ones(n+1,n+1))
    view(0,90)
end
output_file2=strcat('.\figs\random_quad_mesh_mid_x_L',int2str(L),'_n',int2str(n),'_a',num2str(fraction,3));
print('-dpdf',strcat(output_file2,'.pdf'));
print('-dpng',strcat(output_file2,'.png'));
saveas(gcf,strcat(output_file2,'.fig'),'fig');

% keep x=L/2  and y=L/2 as straight lines
if(mod(n,2)==0)
    yy=y;
    yy(:,n/2+1)=L/2;
    figure(3)
    surf(xx,yy,ones(n+1,n+1))
    view(0,90)
end
output_file3=strcat('.\figs\random_quad_mesh_mid_xy_L',int2str(L),'_n',int2str(n),'_a',num2str(fraction,3));
print('-dpdf',strcat(output_file3,'.pdf'));
print('-dpng',strcat(output_file3,'.png'));
saveas(gcf,strcat(output_file3,'.fig'),'fig');

%
% close all;
matID=0;
srcID=0;
%%%%%%%%%%%%%%%%%%%%%
output_file1=strcat(output_file1,'.txt')
fid=fopen(output_file1,'w');
fprintf(fid,'%s\n','polygon');
n2=n*n;
fprintf(fid,'%d\n',n2);
for i=1:n
    for j=1:n
        fprintf(fid,'%d %g %g %g %g %g %g %g %g  %d %d \n',4,x(i,j)    ,y(i,j)    ,...
                                                             x(i+1,j)  ,y(i+1,j)  ,...
                                                             x(i+1,j+1),y(i+1,j+1),...
                                                             x(i,j+1)  ,y(i,j+1)  ,matID,srcID);
    end
end 
fclose(fid)
%%%%%%%%%%%%%%%%%%%%%
output_file2=strcat(output_file2,'.txt')
fid=fopen(output_file2,'w');
fprintf(fid,'%s\n','polygon');
n2=n*n;
fprintf(fid,'%d\n',n2);
for i=1:n
    for j=1:n
        fprintf(fid,'%d %g %g %g %g %g %g %g %g  %d %d \n',4,xx(i,j)    ,y(i,j)    ,...
                                                             xx(i+1,j)  ,y(i+1,j)  ,...
                                                             xx(i+1,j+1),y(i+1,j+1),...
                                                             xx(i,j+1)  ,y(i,j+1)  ,matID,srcID);
    end
end 
fclose(fid)
%%%%%%%%%%%%%%%%%%%%%
output_file3=strcat(output_file3,'.txt')
fid=fopen(output_file3,'w');
fprintf(fid,'%s\n','polygon');
n2=n*n;
fprintf(fid,'%d\n',n2);
for i=1:n
    for j=1:n
        fprintf(fid,'%d %g %g %g %g %g %g %g %g  %d %d \n',4,xx(i,j)    ,yy(i,j)    ,...
                                                             xx(i+1,j)  ,yy(i+1,j)  ,...
                                                             xx(i+1,j+1),yy(i+1,j+1),...
                                                             xx(i,j+1)  ,yy(i,j+1)  ,matID,srcID);
    end
end 
fclose(fid)

