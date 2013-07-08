clear all; close all; clc
% random mesh generated for quadrilaterals

L=100;
n=50;
h=L/n;
xi=linspace(0,L,n+1);
eta=xi;
fraction=0.53*0;

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
output_file1=strcat('.\figs\aa_random_quad_mesh_L',int2str(L),'_n',int2str(n),'_a',num2str(fraction,3));
print('-dpdf',strcat(output_file1,'.pdf'));
print('-dpng',strcat(output_file1,'.png'));
saveas(gcf,strcat(output_file1,'.fig'),'fig');

%---------------------------------------------
% keep x=L/2 a straight line
if(mod(n,2)==0)
    xx=x;
    xx(n/2+1,:)=L/2;
    figure(2)
    surf(xx,y,ones(n+1,n+1))
    view(0,90)
end
output_file2=strcat('.\figs\aa_random_quad_mesh_mid_x_L',int2str(L),'_n',int2str(n),'_a',num2str(fraction,3));
print('-dpdf',strcat(output_file2,'.pdf'));
print('-dpng',strcat(output_file2,'.png'));
saveas(gcf,strcat(output_file2,'.fig'),'fig');

%---------------------------------------------
% keep x=L/2  and y=L/2 as straight lines
if(mod(n,2)==0)
    yy=y;
    yy(:,n/2+1)=L/2;
    figure(3)
    surf(xx,yy,ones(n+1,n+1))
    view(0,90)
end
output_file3=strcat('.\figs\aa_random_quad_mesh_mid_xy_L',int2str(L),'_n',int2str(n),'_a',num2str(fraction,3));
print('-dpdf',strcat(output_file3,'.pdf'));
print('-dpng',strcat(output_file3,'.png'));
saveas(gcf,strcat(output_file3,'.fig'),'fig');

%
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
        fprintf(fid,' %g %g \n %g %g \n %g %g \n %g %g \n',x(i,j)    ,y(i,j)    ,...
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
%%%%%%%%%%%%%%%%%%%%%
output_file2=strcat(output_file2,'.txt')
fid=fopen(output_file2,'w');

fprintf(fid,'# Date: %d/%d/%d   Time: %d:%d\n', mo, da, yr, hr, mi);

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
        fprintf(fid,' %g %g \n %g %g \n %g %g \n %g %g \n',xx(i,j)    ,y(i,j)    ,...
                                                 xx(i+1,j)  ,y(i+1,j)  ,...
                                                 xx(i+1,j+1),y(i+1,j+1),...
                                                 xx(i,j+1)  ,y(i,j+1)  );
    end
end 

fprintf(fid,'# grid vertices (counter-clockwise) \n');
fprintf(fid,'%d\n',(n+1)^2);
for i=1:n+1
    for j=1:n+1
        fprintf(fid,'%g %g \n',xx(i,j),y(i,j) );
    end
end

fclose(fid)
%%%%%%%%%%%%%%%%%%%%%
output_file3=strcat(output_file3,'.txt')
fid=fopen(output_file3,'w');

fprintf(fid,'# Date: %d/%d/%d   Time: %d:%d\n', mo, da, yr, hr, mi);

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
        fprintf(fid,' %g %g \n %g %g \n %g %g \n %g %g \n',xx(i,j)    ,yy(i,j)    ,...
                                                 xx(i+1,j)  ,yy(i+1,j)  ,...
                                                 xx(i+1,j+1),yy(i+1,j+1),...
                                                 xx(i,j+1)  ,yy(i,j+1)  );
    end
end 

fprintf(fid,'# grid vertices (counter-clockwise) \n');
fprintf(fid,'%d\n',(n+1)^2);
for i=1:n+1
    for j=1:n+1
        fprintf(fid,'%g %g \n',xx(i,j),yy(i,j) );
    end
end

fclose(fid)


% % % for BT's code, not used anymore
% % % 
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
% % % output_file2=strcat(output_file2,'.txt')
% % % fid=fopen(output_file2,'w');
% % % fprintf(fid,'%s\n','polygon');
% % % n2=n*n;
% % % fprintf(fid,'%d\n',n2);
% % % for i=1:n
% % %     for j=1:n
% % %         fprintf(fid,'%d %g %g %g %g %g %g %g %g  %d %d \n',4,xx(i,j)    ,y(i,j)    ,...
% % %                                                              xx(i+1,j)  ,y(i+1,j)  ,...
% % %                                                              xx(i+1,j+1),y(i+1,j+1),...
% % %                                                              xx(i,j+1)  ,y(i,j+1)  ,matID,srcID);
% % %     end
% % % end 
% % % fclose(fid)
% % % %%%%%%%%%%%%%%%%%%%%%
% % % output_file3=strcat(output_file3,'.txt')
% % % fid=fopen(output_file3,'w');
% % % fprintf(fid,'%s\n','polygon');
% % % n2=n*n;
% % % fprintf(fid,'%d\n',n2);
% % % for i=1:n
% % %     for j=1:n
% % %         fprintf(fid,'%d %g %g %g %g %g %g %g %g  %d %d \n',4,xx(i,j)    ,yy(i,j)    ,...
% % %                                                              xx(i+1,j)  ,yy(i+1,j)  ,...
% % %                                                              xx(i+1,j+1),yy(i+1,j+1),...
% % %                                                              xx(i,j+1)  ,yy(i,j+1)  ,matID,srcID);
% % %     end
% % % end 
% % % fclose(fid)
