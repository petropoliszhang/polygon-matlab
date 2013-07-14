clear all; close all; clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number of subdivisions of the original rectangle
nc = 5;
% random parameter
a  = 0.05;
% rectangle dimensions
L = 1;
rm = L;
zm = L;

% allocate memory
nm = 2^nc + 1;
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

d=0; % play with this
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
%%%%%%%%%%%%%
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
    % re-order poly ccw
    g=new_ci;
    vv=v(g,:);
    new_connectivity = re_order_ccw(vv,g);
    c{i} = new_connectivity;
end

clc

figure(9);
hold all
for id=1:length(c)
    vv= v(c{id},:);
    xc=mean(vv(:,1));
    yc=mean(vv(:,2));
    vv(end+1,:)=vv(1,:);
    plot(vv(:,1),vv(:,2),'+-')
    %plot(xc,yc,'x');
    str=sprintf('%d',id);
    text(xc,yc,str);
    %     pause
end


tot_area=0;
for iel=1:length(c)
    g=c{iel}(:);
    xx=v(g,1); yy=v(g,2);
    % check orientation, verify area
    [or,ar] = polyorient(xx,yy);
    if(or~=1),
        iel
        g
        or
        ar
        error('orientation problem');
    end
    tot_area=tot_area+ar;
end
fprintf('total area read in geom = %g \n',tot_area);

% decide whether some polygons overlap
for i=1:length(c)
    gi = c{i};
    for j=i+1:length(c)
        gj = c{j};
        % combine the vertex entries
        combined =[gi gj];
        uniq = unique(combined);
        nc=length(combined);
        nu=length(uniq);
        if ( nu> nc )
            error('nu cannot be > than nc');
        end
        
        if ( nu < nc-2)
            % need to investigate, some vertices are common
            vi=v(gi,:);
            vj=v(gj,:);
            % compute intersection of polygons
            [xb,yb]=polybool('intersection',vi(:,1),vi(:,2),vj(:,1),vj(:,2));
            
            % continue to investigate if nonzero intersection area
            if(~isempty(xb))
                % reorder ccw
                xb(end)=[];yb(end)=[];
                [xbb,ybb]=poly2ccw(xb,yb);
                % compute areas
                [orientation_i,area_i]=polyorient(vi(:,1),vi(:,2));
                [orientation_j,area_j]=polyorient(vj(:,1),vj(:,2));
                % compute soustraction of polygons
                if(area_j>area_i)
                    ID_to_modify = j;
                    connectivity_big_poly = gj;
                    vert_big_poly = vj;
                    connectivity_small_poly = gi;
                    vert_small_poly = vi;
                    [xc,yc]=polybool('-',vj(:,1),vj(:,2),vi(:,1),vi(:,2));
                else
                    ID_to_modify = i;
                    connectivity_big_poly = gi;
                    vert_big_poly = vi;
                    connectivity_small_poly = gj;
                    vert_small_poly = vj;
                    [xc,yc]=polybool('-',vi(:,1),vi(:,2),vj(:,1),vj(:,2));
                end
                % reorder ccw
                xc(end)=[];yc(end)=[];
                [xcc,ycc]=poly2ccw(xc,yc);
                vcc=[xcc ycc];
                % find the common points between the subtraction polygon
                % and the largest polygon
                new_connectivity = zeros(1,length(xcc));
                for k=1:length(xcc)
                    vk=vcc(k,:);
                    dif = sum(sqrt( (vert_big_poly-kron(ones(length(connectivity_big_poly),1),vk)).^2 ) , 2);
                    ind=find(dif<1e-10);
                    if(length(ind)==1)
                        new_connectivity(k) = connectivity_big_poly(ind);
                    else
                        new_connectivity(k) = -k;
                    end
                end
                % find the common points between the subtraction polygon
                % and the smallest polygon
                ind_neg = find( new_connectivity <0 );
                for k=1:length(ind_neg)
                    vk=vcc(ind_neg(k),:);
                    dif = sum(sqrt( (vert_small_poly-kron(ones(length(connectivity_small_poly),1),vk)).^2 ) , 2);
                    ind=find(dif<1e-10);
                    if(length(ind)==1)
                        new_connectivity(ind_neg(k)) = connectivity_small_poly(ind);
                    else
                        warning('missing points ...');
                    end
                end
                
                % overwrite the definition of the largest poly
                c{ID_to_modify} =   new_connectivity;
                
            end
        end
        
    end
end

% either one is fine
% result = c(~cellfun('isempty',c))
% cc(cellfun(@isempty,cc)) = []

figure(11)
hold all
for id=1:length(c)
    patch(v(c{id},1),v(c{id},2),id); % use color i.
end

figure(10)
hold all
for id=1:length(c)
    vv= v(c{id},:);
    xc=mean(vv(:,1));
    yc=mean(vv(:,2));
    vv(end+1,:)=vv(1,:);
    plot(vv(:,1),vv(:,2),'+-')
    plot(xc,yc,'x');
end


tot_area=0;
for iel=1:length(c)
    g=c{iel}(:);
    xx=v(g,1); yy=v(g,2);
    % check orientation, verify area
    [or,ar] = polyorient(xx,yy);
    if(or~=1),
        error('orientation problem');
    end
    tot_area=tot_area+ar;
end
fprintf('total area read in geom = %g \n',tot_area);

%%%%%%%%%%%%%
output_file1=strcat('.\figs\shestakov_poly_mesh_nc',int2str(nc),'_a',num2str(a,3));
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


%%%%%%%%%%
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


return