clear all; close all; clc;

L=1;
n=100;
h=L/n;
xi=linspace(-h,L+h,n+1);
eta=xi;
fraction=0.1;

ind=0;
for i=1:n+1
    for j=1:n+1
        ind=ind+1;
        x(ind) = xi(i) +fraction*sin(2*pi*xi(i)/L)*sin(2*pi*eta(j)/L);
        y(ind) = eta(j)+fraction*sin(2*pi*xi(i)/L)*sin(2*pi*eta(j)/L);
    end
end

figure(1)
voronoi(x,y);
figure(2)
[vx,vy]=voronoi(x,y);
line(vx,vy,'Color','b'); % plot Voronoi edges

figure(3)
[v,c]=voronoin([x' y']);
for i = 1:length(c)
    if all(c{i}~=1)   % If at least one of the indices is 1,
        % then it is an open region and we can't
        % patch that.
        patch(v(c{i},1),v(c{i},2),i); % use color i.
    end
end

% copy vertex list
v2=v;

% loop through voronoi cells
for ic = 1:length(c)
    % skip cells with point at infinity
    if( any(c{ic}==1) )
        fprintf('Point at infinity for cell %d \n',ic);
        continue
    end
    outside=[]; inside=[];
    nc = length(c{ic});
    % loop over vertices of current voronoi cell
    for k=1:nc
        A=v(c{ic}(k),:);
        if( A(1)<0 | A(1)>L | A(2)<0 | A(2)>L )
            outside=[outside c{ic}(k)];
        else
            inside=[inside c{ic}(k)];
        end
    end
    % if all points are outside, skip
    if( length(outside)==nc )
        fprintf('All points are outside for cell %d \n',ic);
        continue
    end
    % if all points are inside, skip
    if( length(inside)==nc )
        fprintf('All points are inside for cell %d \n',ic);
        continue
    end

    % re-order the inside vertices to be in the same order as in c{ic};
%     inside
%     c{ic}
    ia=intersect_jcr(c{ic},inside,'sort');
    %     [C,ia,ib]=intersect(c{ic},inside)
    inside=c{ic}(ia); % same as inside(ib)

    % downselect to the inside points plus the first 2 outside points
    vertices_to_keep = find_vertices_to_keep(c{ic},inside);
    if(vertices_to_keep(1)==vertices_to_keep(end))
        vertices_to_keep(end) = length(v)+1;
        v(end+1,:) =v(vertices_to_keep(1),:);
        v2(end+1,:)=v(vertices_to_keep(1),:);
        ind=find(c{ic}(:)==vertices_to_keep(1));
        nc=nc+1;
        for i=nc:-1:ind+1
            c{ic}(i)=c{ic}(i-1);
        end
        c{ic}(ind)=vertices_to_keep(end);
    end

    ia=intersect_jcr(c{ic},vertices_to_keep);
    %     [C,ia]=intersect_jcr(c{ic},vertices_to_keep);
    ia=sort(ia);
    for k=nc:-1:1
        if( ~isempty(ia) & k==ia(end))
            ia(end)=[];
            continue;
        else
            v2(c{ic}(k),:)=NaN;
            c{ic}(k)=[];
        end
    end

    in = find( c{ic}(:) == vertices_to_keep(2) );
    ex = find( c{ic}(:) == vertices_to_keep(1) );
    interior_pt = v(c{ic}(in),:);
    exterior_pt = v(c{ic}(ex),:);
    [x,y]=find_intersect_with_square(interior_pt,exterior_pt,L);
    v2(c{ic}(ex),:)=[x y];

    in = find( c{ic}(:) == vertices_to_keep(end-1) );
    ex = find( c{ic}(:) == vertices_to_keep(end  ) );
    interior_pt = v(c{ic}(in),:);
    exterior_pt = v(c{ic}(ex),:);
    [x,y]=find_intersect_with_square(interior_pt,exterior_pt,L);
    v2(c{ic}(ex),:)=[x y];

%     fprintf('here -------------------------------------------\n\n');


end


for ic = length(c):-1:1
    if( any(c{ic}==1) )
        c{ic}=[];
        continue
    end
    outside=[]; inside=[];
    nc = length(c{ic});
    % loop over vertices of current voronoi cell
    for k=1:nc
        A=v(c{ic}(k),:);
        if( A(1)<0 | A(1)>L | A(2)<0 | A(2)>L )
            outside=[outside c{ic}(k)];
        else
            inside=[inside c{ic}(k)];
        end
    end
    % if all points are outside, skip
    if( length(outside)==nc )
        c{ic}=[];
        continue
    end
end

% fix corners
mTL=1e4; mBL=1e4; mTR=1e4; mBR=1e4;
TL=[L 0]; BL=[0 0];
TR=[L L]; BR=[0 L];

for ic = 1:length(c)
    if(length(c{ic})>0)
        d = norm(mean(v2(c{ic}(:),:))-TL);
        if(d<mTL)
            mTL=d; ID_TL=ic;
        end
        d = norm(mean(v2(c{ic}(:),:))-TR);
        if(d<mTR)
            mTR=d; ID_TR=ic;
        end
        d = norm(mean(v2(c{ic}(:),:))-BL);
        if(d<mBL)
            mBL=d; ID_BL=ic;
        end
        d = norm(mean(v2(c{ic}(:),:))-BR);
        if(d<mBR)
            mBR=d; ID_BR=ic;
        end
    end
end
mean(v2(c{ID_TL}(:),:))
mean(v2(c{ID_BL}(:),:))
mean(v2(c{ID_TR}(:),:))
mean(v2(c{ID_BR}(:),:))

% BL
ind=find( abs(v2(c{ID_BL}(:)                          ,1)-0)<1e-10 ); %x
if(       abs(v2(c{ID_BL}(mod(ind,length(c{ID_BL}))+1),2)-0)<1e-10 )  %y
    skip=1
else
    skip=0;
end
nc=length(c{ID_BL});
for i=nc+1:-1:ind+1+skip
    c{ID_BL}(i)=c{ID_BL}(i-1);
end
c{ID_BL}(ind+skip)=length(v2)+1;
v2(end+1,:)=[0 0];

% TL
ind=find( abs(v2(c{ID_TL}(:)                          ,1)-L)<1e-10 );
if(       abs(v2(c{ID_TL}(mod(ind,length(c{ID_TL}))+1),2)-0)<1e-10 )
    skip=1
else
    skip=0;
end
nc=length(c{ID_TL});
% for i=nc+1:-1:ind+2
%     c{ID_TL}(i)=c{ID_TL}(i-1);
% end
for i=nc+1:-1:ind+1+skip
    c{ID_TL}(i)=c{ID_TL}(i-1);
end
c{ID_TL}(ind+skip)=length(v2)+1;
v2(end+1,:)=[L 0];

% TR
ind=find( abs(v2(c{ID_TR}(:)                          ,2)-L)<1e-10 );
if(       abs(v2(c{ID_TR}(mod(ind,length(c{ID_TR}))+1),1)-L)<1e-10 )
    skip=1
else
    skip=0;
end
nc=length(c{ID_TR});
for i=nc+1:-1:ind+1+skip
    c{ID_TR}(i)=c{ID_TR}(i-1);
end
c{ID_TR}(ind+skip)=length(v2)+1;
v2(end+1,:)=[L L];

% BR
ind=find( abs(v2(c{ID_BR}(:)                          ,2)-L)<1e-10 );
if(       abs(v2(c{ID_BR}(mod(ind,length(c{ID_BR}))+1),1)-0)<1e-10 )
    skip=1
else
    skip=0;
end
nc=length(c{ID_BR});
for i=nc+1:-1:ind+1+skip
    c{ID_BR}(i)=c{ID_BR}(i-1);
end
c{ID_BR}(ind+skip)=length(v2)+1;
v2(end+1,:)=[0 L];


disp('plotting')

figure(4)
npoly=0;
for i = 1:length(c)
    if all(c{i}~=1)   % If at least one of the indices is 1,
        % then it is an open region and we can't
        % patch that.
        %         if(v2(c{i},1)==NaN)
        %             continue
        %         end
        [ii,jj]=size(v2(c{i},:));
        if(ii>0)
            patch(v2(c{i},1),v2(c{i},2),i); % use color i.
            npoly=npoly+1;
        end
    else
        disp('pbpbpb')
        i
        c{i}
        error('still cells that are infite regions ... %d',i);
    end
end

%%%%%%%%%%%%%
output_file1=strcat('.\figs\smooth_poly_mesh_L',int2str(L),'_n',int2str(n),'_a',num2str(fraction,3));
print('-dpdf',strcat(output_file1,'.pdf'));
print('-dpng',strcat(output_file1,'.png'));
saveas(gcf,strcat(output_file1,'.fig'),'fig');
%%%%%%%%%%
matID=0;
srcID=0;
%%%%%%%%%%
output_file1=strcat(output_file1,'.txt')
fid=fopen(output_file1,'w');
fprintf(fid,'%s\n','polygon');
fprintf(fid,'%d\n',npoly);
for i = 1:length(c)
    [ii,jj]=size(v2(c{i},:));
    if(ii>0)
        nvert=ii;
        fprintf(fid,'%d  ',nvert);
        % extract all coord
        aa=v2(c{i},:);
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
        fprintf(fid,'%g  ',aa');
        fprintf(fid,'  %d %d \n',matID,srcID);
    end
end
fclose(fid);
%%%%%%%%%%%%%

% clear all
% load clown
% image(X); colormap(map); axis off
% siz = size(X);
% N = 100; % number of points
% xp = ceil(rand(N,1)*(siz(2)-1)); % generate random x coordinate
% yp = ceil(rand(N,1)*(siz(1)-1)); % generate random y coordinate
% [vx,vy] = voronoi(xp,yp);
% line(vx,vy,'Color','y'); % plot Voronoi edges