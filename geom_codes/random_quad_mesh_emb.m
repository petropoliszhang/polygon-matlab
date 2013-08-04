clear all; close all; clc
% random mesh generated for quadrilaterals

logi_save_plots = false;
logi_write_file = false;

L=1;
n=2^3; %= 2^4
h=L/n;
xi=linspace(0,L,n+1);
eta=xi;
fraction=0.50;
BT=false;

x=zeros(n+1,n+1);
y=x;
for i=1:n+1
    ax=1;   if(i==1||i==n+1),ax=0;end
    for j=1:n+1
        ay=1;   if(j==1||j==n+1),ay=0;end
        x(i,j) = xi(i) + (2*rand(1,1)-1)*h*fraction*ax;
        y(i,j) = eta(j)+ (2*rand(1,1)-1)*h*fraction*ay;
    end
end

figure(1)
surf(x,y,ones(n+1,n+1)*0,'FaceColor','white')
view(0,90)

% n = 2^nc --> log_2(n) = nc
nc = log2(n);
integerTest = ( nc == floor(nc) );
if ~integerTest
    error('must be embedded');
end
figure(2);
q1=(ceil(nc/2));
q2=ceil(nc/q1);
for k=1:nc
    skip=2^(k-1);
    subplot(q1,q2,k);
    nm_ = 2^(nc-k+1)+1;
    surf(x(1:skip:end,1:skip:end),y(1:skip:end,1:skip:end),ones(nm_,nm_),'FaceColor','white');
    set(gca,'PlotBoxAspectRatio',[5 5 1])
    view(0,90);
end

if(logi_save_plots)
    output_file1=strcat('.\figs\random_quad_mesh_L',int2str(L),'_n',int2str(n),'_emb','_a',num2str(fraction,3));
    print('-dpdf',strcat(output_file1,'.pdf'));
    print('-dpng',strcat(output_file1,'.png'));
    saveas(gcf,strcat(output_file1,'.fig'),'fig');
end

%---------------------------------------------
%
if(~logi_write_file), return; end;
%
%---------------------------------------------
% save txt file
matID=1;
srcID=1;

% date and time;
[yr, mo, da, hr, mi, s] = datevec(now);

% save embedded meshes
x_sav=x;
y_sav=y;
for k=1:nc

    nm_ = 2^(nc-k+1)+1;
    skip=2^(k-1);
    x=x_sav(1:skip:end,1:skip:end);
    y=y_sav(1:skip:end,1:skip:end);
    %%%%%%%%%%%%%%%%%%%%%
    output_file_k=strcat('.\figs\random_quad_mesh_L',int2str(L),...
        '_nc',int2str(nc),'_emb',int2str(nc-k+1),'_a',num2str(fraction,3),'.txt');
    fid=fopen(output_file_k,'w');


    fprintf(fid,'# Date: %d/%d/%d   Time: %d:%d\n', mo, da, yr, hr, mi);

    fprintf(fid,'# dimensions \n');
    fprintf(fid,'%g %g \n',L,L);

    n = (nm_-1);
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

    fclose(fid);
    %%%%
end % embedded meshes save


