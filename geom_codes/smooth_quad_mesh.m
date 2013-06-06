clear all; close all; clc

L=1;
n=10;
h=L/n;
xi=linspace(0,L,n+1);
eta=xi;
fraction=0.12;

for i=1:n+1
%     ax=1;   if(i==1|i==n+1),ax=0;end
    for j=1:n+1
%         ay=1;   if(j==1|j==n+1),ay=0;end
        x(i,j) = xi(i) + fraction*sin(2*pi*xi(i)/L)*sin(2*pi*eta(j)/L);
        y(i,j) = eta(j)+ fraction*sin(2*pi*xi(i)/L)*sin(2*pi*eta(j)/L);
    end
end

figure(1)
surf(x,y,ones(n+1,n+1))
view(0,90)


figure(1)
surf(x,y,ones(n+1,n+1))
view(0,90)

output_file1=strcat('.\figs\smooth_quad_mesh_L',int2str(L),'_n',int2str(n),'_a',num2str(fraction,3));
print('-dpdf',strcat(output_file1,'.pdf'));
print('-dpng',strcat(output_file1,'.png'));
saveas(gcf,strcat(output_file1,'.fig'),'fig');

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