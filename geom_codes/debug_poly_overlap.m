% testing this
clc; clear all; close all;

load debug_poly_bad.mat;
c=c_ori;
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

for id=1:length(c)
    figure(10)
    hold all
    vv= v(c{id},:);
    xc=mean(vv(:,1));
    yc=mean(vv(:,2));
    vv(end+1,:)=vv(1,:);
    plot(vv(:,1),vv(:,2),'+-')
    plot(xc,yc,'x');
    %     pause
    figure(11)
    patch(v(c{id},1),v(c{id},2),id); % use color i.
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
