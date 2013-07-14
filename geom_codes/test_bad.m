clear all; close all; clc;

load bad.mat;

c=c_ori;

% clean up duplicated vertices in single polygon definition
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

%%%%%%%%%%%%%
% find Inf and NaN
nv=length(v(:,1));
TF=isfinite(v);
[ii,jj]=find(TF==0);
inf_nan=unique(ii);
% delete in descending order
inf_nan=sort(inf_nan,'descend');
for i=1:length(inf_nan)
    if(i==1)
        fprintf('the following point will be omitted b/c the are NOT finite \n');
    end
    fprintf(' vertex # %d, x-value=%g, y-value=%g \n',inf_nan(i),v(inf_nan(i),1),v(inf_nan(i),2) );
    %     v(del(i),:)=[];
end
%
% find duplicated vertices in v
del=[];
keep=[];
for k=1:nv
    if( ~isempty(find(inf_nan==k)) )
        fprintf('skipping inf_nan point %d \n',k);
        continue
    end
    vk=v(k,1:2);
    % get difference
    aux = v - kron(ones(nv,1),vk);
    % get vector of norm
    aa=sqrt(aux(:,1).^2+aux(:,2).^2);
    ind=find(aa<1e-12);
    len=length(ind);
    if(len==0),
        vk
        error('vk not found ????');
    end
    if(len >1),
        warning('vk duplicated');
        vk;
        v(ind,:);
        k
        ind
        ind2=find(ind>k);
        ind(ind2)
        del=[del ind(ind2)];
        if(~isempty(ind2))
            keep=[keep k];
        end
    end
end
% del=unique(del)
del
keep
% delete in descending order
[del,isort]=sort(del,'descend');
keep=keep(isort);
del
keep

ID_to_look_for=[];

for i=1:length(del)
    if(i==1)
        fprintf('the following point will be omitted b/c they are DUPLICATES \n');
    end
    fprintf(' vertex # %d, x-value=%g, y-value=%g \n',del(i),v(del(i),1),v(del(i),2) );
    %     v(del(i),:)=[];
    % update connectivity
    for id=1:length(c)
        g=c{id};
        ind=find(g==del(i));
        if(~isempty(ind))
            g;
            g(ind)=keep(i);
            g;
            ID_to_look_for=[ID_to_look_for id];
        end
        c{id}=g;
    end
end
ID_to_look_for

%%%%%%%%%%%%%
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

return

