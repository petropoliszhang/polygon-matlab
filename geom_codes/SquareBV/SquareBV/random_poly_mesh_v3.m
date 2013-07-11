% clear all; 
close all; clc;

logi_save = true;

L=10;
n=2;
h=L/n;
eps=h/10*0;
xi=linspace(eps,L-eps,n+1);
eta=xi;
fraction=0.5;


ind=0;
for i=1:n+1
    ax=1;   if(i==1|i==n+1),ax=0;end
    for j=1:n+1
        ay=1;   if(j==1|j==n+1),ay=0;end
        ind=ind+1;
        x(ind) = xi(i) + (2*rand(1,1)-1)*h*fraction*ax;
        y(ind) = eta(j)+ (2*rand(1,1)-1)*h*fraction*ay;
    end
end

%%%%%%%%%%%%%
SquareBV(x,y,0)
%%%%%%%%%%%%%


return




%%%%%%%%%%%%%
output_file1=strcat('.\figs\random_poly_mesh_L',int2str(L),'_n',int2str(n),'_a',num2str(fraction,3));
%%%%%%%%%%%%%

%%%%%%%%%%%%%
% save_matlab_plot = input('save matlab plots? enter 1 for yes, 0 otherwise \n');
save_matlab_plot=0;
if save_matlab_plot~=0 || save_matlab_plot ~=1
    save_matlab_plot=0;
end
if save_matlab_plot
    print('-dpdf',strcat(output_file1,'.pdf'));
    print('-dpng',strcat(output_file1,'.png'));
    saveas(gcf,strcat(output_file1,'.fig'),'fig');
end
%%%%%%%%%%

%%%%%%%%%%%%%
% save txt file
matID=1;
srcID=1;
% date and time;
[yr, mo, da, hr, mi, s] = datevec(now);
%
output_file1=strcat(output_file1,'.txt')
fid=fopen(output_file1,'w');
fprintf(fid,'# Date: %d/%d/%d   Time: %d:%d\n', mo, da, yr, hr, mi);
fprintf(fid,'# dimensions \n');
fprintf(fid,'%g %g \n',L,L);

ncells = length(c);

fprintf(fid,'# connectivity \n');
fprintf(fid,'%d\n',ncells);
skip = 0 ;
for i = 1:length(c)
    nvert=length(c{i});
    fprintf(fid,'%d  ',nvert);
    for k=1:nvert
        fprintf(fid,'%d  ',skip+k);
    end
    fprintf(fid,'  %d %d \n',matID,srcID);
    skip = skip + nvert;
end

nvert_total = skip;

fprintf(fid,'# DG vertices (counter-clockwise) \n');
fprintf(fid,'%d\n',nvert_total);
for i = 1:length(c)
    % extract all coord
    aa=v(c{i},:);
    x=aa(:,1);
    y=aa(:,2);
    % Find the centroid:
    cx = mean(x);
    cy = mean(y);
    % Step 2: Find the angles:
    ang = atan2(y - cy, x - cx);
    % Step 3: Find the correct sorted order:
    [dummy, order] = sort(ang);
    % Step 4: Reorder the coordinates:
    x = x(order);
    y = y(order);
    aa(:,1)=x;
    aa(:,2)=y;
    fprintf(fid,'%g  %g  \n',aa');
end

fprintf(fid,'# grid vertices \n');
fprintf(fid,'%d\n',length(v(:,1)));
for k=1:length(v(:,1))
    fprintf(fid,'%g %g \n',v(k,1),v(k,2) );
end

fclose(fid);

