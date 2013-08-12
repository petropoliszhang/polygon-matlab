clear all; close all; clc

L=1;
n=60;
h=L/n;
xi=linspace(0,L,n+1);
eta=xi;
fraction=1/pi/2+.002;

logi_save_plots = false;
logi_write_file = false;

for i=1:n+1
%     ax=1;   if(i==1|i==n+1),ax=0;end
    for j=1:n+1
%         ay=1;   if(j==1|j==n+1),ay=0;end
        x(i,j) = xi(i) + L*fraction*sin(2*pi*xi(i)/L)*sin(2*pi*eta(j)/L);
        y(i,j) = eta(j)+ L*fraction*sin(2*pi*xi(i)/L)*sin(2*pi*eta(j)/L);
    end
end

figure(1)
surf(x,y,ones(n+1,n+1))
view(0,90)


figure(1)
surf(x,y,ones(n+1,n+1))
view(0,90)

output_file1=strcat('.\figs\smooth_quad_mesh_L',int2str(L),'_n',int2str(n),'_a',num2str(fraction,3));
if(logi_save_plots)
    print('-dpdf',strcat(output_file1,'.pdf'));
    print('-dpng',strcat(output_file1,'.png'));
    saveas(gcf,strcat(output_file1,'.fig'),'fig');
end

%---------------------------------------------

if(~logi_write_file), return; end

%---------------------------------------------
% save txt file
matID=1;
srcID=1;
% date and time;
[yr, mo, da, hr, mi, s] = datevec(now);
%%%%%%%%%%%%%%%%%%%%%
output_file1=strcat(output_file1,'.txt')
fid=fopen(output_file1,'w');

fprintf(fid,'# Date: %d/%d/%d   Time: %d:%d\n', mo, da, yr, hr, mi);

fprintf(fid,'# dimensions \n');
fprintf(fid,'%g %g \n',L,L);

ncells = n*n;

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
for i=1:n
    for j=1:n
        fprintf(fid,' %g %g \n %g %g \n %g %g \n %g %g \n',...
            x(i,j)    ,y(i,j)    ,...
            x(i+1,j)  ,y(i+1,j)  ,...
            x(i+1,j+1),y(i+1,j+1),...
            x(i,j+1)  ,y(i,j+1)  );
    end
end

fprintf(fid,'# grid vertices (counter-clockwise) \n');
fprintf(fid,'%d\n',(n+1)^2);
for i=1:n+1
    for j=1:n+1
        fprintf(fid,'%g %g \n',x(i,j),y(i,j) );
    end
end

fclose(fid)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% close all;
% % % matID=0;
% % % srcID=0;
% % % %%%%%%%%%%%%%%%%%%%%%
% % % output_file1=strcat(output_file1,'.txt')
% % % fid=fopen(output_file1,'w');
% % % fprintf(fid,'%s\n','polygon');
% % % n2=n*n;
% % % fprintf(fid,'%d\n',n2);
% % % for i=1:n
% % %     for j=1:n
% % %         fprintf(fid,'%d %g %g %g %g %g %g %g %g  %d %d \n',4,x(i,j)    ,y(i,j)    ,...
% % %                                                              x(i+1,j)  ,y(i+1,j)  ,...
% % %                                                              x(i+1,j+1),y(i+1,j+1),...
% % %                                                              x(i,j+1)  ,y(i,j+1)  ,matID,srcID);
% % %     end
% % % end 
% % % fclose(fid)
% % % %%%%%%%%%%%%%%%%%%%%%