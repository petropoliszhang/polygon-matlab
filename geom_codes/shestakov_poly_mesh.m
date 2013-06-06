clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number of subdivisions of the original rectangle
nc = 6;
% random parameter
a  = 0.25;
% rectangle dimensions
L = 1;
rm = L;
zm = L;

% allocate memory
nm = 2^nc + 1
ind = nm-1;
r=zeros(nm,nm);
z=zeros(nm,nm);

% initialize 4 corners
for i = 0:1
    k = 1 + i*ind;
    for j = 0:1
        l = 1 + j*ind;
        r(k,l) = L*(1-i);
        z(k,l) = L*j;
    end
end

d=0.09; % play with this
r(1,1)     = L+d;
r(1,end)   = L+d;
r(end,1)   =  -d;
r(end,end) =  -d;

z(1,1)     =  -d;
z(1,end)   = L+d;
z(end,1)   =  -d;
z(end,end) = L+d;

% fill in the rest of the points
for nl = 0:nc-1
    nn = 2^nl;
    inc = ind/nn;
    inch = inc/2;
    for k = 1:nn
        k1 = 1+(k-1)*inc;
        k2 = k1+inc;
        k3 = k1+inch;
        for l = 1:nn
            l1 = 1+(l-1)*inc;
            l2 = l1+inc;
            l3 = l1+inch;
            if (l==1)
                ar = a+rand(1,1)*(1-2*a);
                r(k3,l1) = ar*r(k1,l1)+(1-ar)*r(k2,l1);
                z(k3,l1) = ar*z(k1,l1)+(1-ar)*z(k2,l1);
            end
            ar = a+rand(1,1)*(1-2*a);
            r(k2,l3) = ar*r(k2,l1)+(1-ar)*r(k2,l2);
            z(k2,l3) = ar*z(k2,l1)+(1-ar)*z(k2,l2);
            ar = a+rand(1,1)*(1-2*a);
            r(k3,l2) = ar*r(k1,l2)+(1-ar)*r(k2,l2);
            z(k3,l2) = ar*z(k1,l2)+(1-ar)*z(k2,l2);
            if (k==1)
                ar = a+rand(1,1)*(1-2*a);
                r(k1,l3) = ar*r(k1,l1)+(1-ar)*r(k1,l2);
                z(k1,l3) = ar*z(k1,l1)+(1-ar)*z(k1,l2);
            end
            ar = a+rand(1,1)*(1-2*a);
            br = a+rand(1,1)*(1-2*a);
            r1 = r(k1,l1);
            r2 = r(k2,l1);
            r3 = r(k2,l2);
            r4 = r(k1,l2);
            z1 = z(k1,l1);
            z2 = z(k2,l1);
            z3 = z(k2,l2);
            z4 = z(k1,l2);
            % check for boomerang zones
            det2 = (r2-r1)*(z3-z2)-(r3-r2)*(z2-z1);
            det3 = (r3-r2)*(z4-z3)-(r4-r3)*(z3-z2);
            det4 = (r4-r3)*(z1-z4)-(r1-r4)*(z4-z3);
            det1 = (r1-r4)*(z2-z1)-(r2-r1)*(z1-z4);
            if (det2>0)
                d = (r4-r3)*(z2-z1)-(r2-r1)*(z4-z3);
                r3p = ((r2-r1)*(r4*z3-r3*z4)-(r4-r3)*(r2*z1-r1*z2))/d;
                z3p = ((z2-z1)*(r4*z3-r3*z4)-(z4-z3)*(r2*z1-r1*z2))/d;
                d = (r4-r1)*(z2-z3)-(r2-r3)*(z4-z1);
                r1p = ((r2-r3)*(r4*z1-r1*z4)-(r4-r1)*(r2*z3-r3*z2))/d;
                z1p = ((z2-z3)*(r4*z1-r1*z4)-(z4-z1)*(r2*z3-r3*z2))/d;
                r3 = r3p;
                z3 = z3p;
                r1 = r1p;
                z1 = z1p;
            elseif (det3>0)
                d = (r1-r4)*(z3-z2)-(r3-r2)*(z1-z4);
                r4p = ((r3-r2)*(r1*z4-r4*z1)-(r1-r4)*(r3*z2-r2*z3))/d;
                z4p = ((z3-z2)*(r1*z4-r4*z1)-(z1-z4)*(r3*z2-r2*z3))/d;
                d = (r1-r2)*(z3-z4)-(r3-r4)*(z1-z2);
                r2p = ((r3-r4)*(r1*z2-r2*z1)-(r1-r2)*(r3*z4-r4*z3))/d;
                z2p = ((z3-z4)*(r1*z2-r2*z1)-(z1-z2)*(r3*z4-r4*z3))/d;
                r4 = r4p;
                z4 = z4p;
                r2 = r2p;
                z2 = z2p;
            elseif (det4>0)
                d = (r2-r1)*(z4-z3)-(r4-r3)*(z2-z1);
                r1p = ((r4-r3)*(r2*z1-r1*z2)-(r2-r1)*(r4*z3-r3*z4))/d;
                z1p = ((z4-z3)*(r2*z1-r1*z2)-(z2-z1)*(r4*z3-r3*z4))/d;
                d = (r2-r3)*(z4-z1)-(r4-r1)*(z2-z3);
                r3p = ((r4-r1)*(r2*z3-r3*z2)-(r2-r3)*(r4*z1-r1*z4))/d;
                z3p = ((z4-z1)*(r2*z3-r3*z2)-(z2-z3)*(r4*z1-r1*z4))/d;
                r1 = r1p;
                z1 = z1p;
                r3 = r3p;
                z3 = z3p;
            elseif (det1>0)
                d = (r3-r2)*(z1-z4)-(r1-r4)*(z3-z2);
                r2p = ((r1-r4)*(r3*z2-r2*z3)-(r3-r2)*(r1*z4-r4*z1))/d;
                z2p = ((z1-z4)*(r3*z2-r2*z3)-(z3-z2)*(r1*z4-r4*z1))/d;
                d = (r3-r4)*(z1-z2)-(r1-r2)*(z3-z4);
                r4p = ((r1-r2)*(r3*z4-r4*z3)-(r3-r4)*(r1*z2-r2*z1))/d;
                z4p = ((z1-z2)*(r3*z4-r4*z3)-(z3-z4)*(r1*z2-r2*z1))/d;
                r2 = r2p;
                z2 = z2p;
                r4 = r4p;
                z4 = z4p;
            end
            r(k3,l3) = ar*br*r1 + ar*(1-br)*r4 + (1-ar)*br*r2 + (1-ar)*(1-br)*r3;
            z(k3,l3) = ar*br*z1 + ar*(1-br)*z4 + (1-ar)*br*z2 + (1-ar)*(1-br)*z3;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=reshape(r, nm^2,1);x=x';
y=reshape(z, nm^2,1);y=y';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
output_file1=strcat('.\figs\shestakov_poly_mesh_nc',int2str(nc),'_a',num2str(a,3));
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