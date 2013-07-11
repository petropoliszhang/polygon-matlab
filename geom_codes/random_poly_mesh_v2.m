% clear all; 
close all; clc;

logi_save = true;

L=10;
n=2;
h=L/n;
eps=h/10*0;
xi=linspace(eps,L-eps,n+1);
eta=xi;
fraction=0.85;


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
x=[            0
            0
            0
       6.2601
       1.5599
       6.0209
           10
           10
           10];
y=[            0
       9.2129
           10
            0
       8.9765
           10
            0 
       5.3712
           10];
[v,c]=VoronoiLimit(x,y,[0 L 0 L]);
%%%%%%%%%%%%%
% clean up duplicated vertices
for i=1:length(c)
    nv=length(c{i});
    new_ci=c{i}(1);
    for k=2:nv
        vk=c{i}(k);
        ind=find(new_ci==vk);
        if(isempty(ind))
            new_ci = [ new_ci vk ];
        end
    end
    c{i} = new_ci;
    % ccw
    g=c{i}(:);
    vv=v(g,:);
    [x2, y2] = poly2ccw(vv(:,1), vv(:,2));

    new_ci=[];
    for k=1:length(x2)
        vk=[x2(k) y2(k)];
        dif = sum(sqrt( (vv-kron(ones(length(x2),1),vk)).^2 ) , 2);
        ind=find(dif<1e-10);
        if(length(ind)~=1)
            ind
            error('the vk vertex should have showed up once only....');
        end
        new_ci =[ new_ci c{i}(ind)];
    end
    c{i} = new_ci;
    
end

tot_area=0;
for iel=1:length(c)
    g=c{iel}(:);
    xx=v(g,1); yy=v(g,2);
    % check orientation, verify area
    [or,ar] = polyorient(xx,yy);
    ar
     if(or~=1), 
         error('orientation problem'); 
     end
    tot_area=tot_area+ar;

end
fprintf('total area read in geom = %g \n',tot_area);
return

%     for k=1:nv
%         vk=vv(k,:);
%         dif = sum(sqrt( (vv-kron(ones(nv,1),vk)).^2 ) , 2);
%         ind=find(dif<1e-10);
%         if(isempty(ind))
%             error('the vk vertex should have showed up....')
%         end
%         for q=1:length(ind)
%             if(ind(q)~=k)
%                 v(ind(q),:)=[];
%             end
%         end
%     end

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

