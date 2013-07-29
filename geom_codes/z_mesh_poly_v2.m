% clear all;
close all; clc;
% Z-mesh generated for polygons

logi_save = true;

L=1;
n=20;
if(mod(n,20)~=0)
    error('n must be a multiple of 20 for the z-mesh')
else
    nn=n/20;
    nsub=[4 7 13 16 20]*nn +1;
end

fraction=0.5;
if(fraction<eps || fraction > 1-eps)
    error('wrong fraction');
end

slope1 = (1-2*fraction)/0.15;
slope2 = (2*fraction-1)/0.30;

y1=@(x) slope1*(x-0.20*L)+fraction*L;
y2=@(x) slope2*(x-0.35*L)+(1-fraction)*L;
y3=@(x) slope1*(x-0.65*L)+fraction*L;

% uniform subdivision along x
xi=linspace(0,L,n+1);

ind=0;
x=zeros((n+1)^2,1);
y=x;

for i=1:n+1
    % five different y-spacing, depending on the zone
    if(i<=nsub(1))
        yy=L*fraction;
    elseif(i<=nsub(2))
        yy=y1(xi(i));
    elseif(i<=nsub(3))
        yy=y2(xi(i));
    elseif(i<=nsub(4))
        yy=y3(xi(i));
    elseif(i<=nsub(5))
        yy=L*(1-fraction);
    else
        error('not in nsub!!!');
    end
    aux1 = linspace(0,yy,n/2+1);
    aux2 = linspace(yy,L,n/2+1); aux2(1)=[];
    aux = [aux1 aux2];
    
    % finish up
    for j=1:n+1
        ind=ind+1;
        % x-axis
        x(ind) = xi(i);
        y(ind) = aux(j);
    end
end

%%%%%%%%%%%%%
[v,c]=VoronoiLimit(x,y,[0 L 0 L]);
clc
%%%%%%%%%%%%%
% c_ori=c;
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
        k;
        ind;
        ind2=find(ind>k);
        ind(ind2);
        del;
        del=[del (ind(ind2))'];
        if(~isempty(ind2))
            keep=[keep k*ones(1,length(ind2))];
        end
        if(length(del)~=length(keep))
            error('length(del)~=length(keep)')
        end
    end
end
% del=unique(del)
% del
% keep
% delete in descending order
[del,isort]=sort(del,'descend');
keep=keep(isort);
% del
% keep

% ID_to_look_for=[];

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
            g(ind)=keep(i);
            %             ID_to_look_for=[ID_to_look_for id];
        end
        c{id}=g;
    end
end
% ID_to_look_for

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
        [xx yy]
        error('orientation problem');
    end
    tot_area=tot_area+ar;
end
fprintf('total area read in geom = %g \n',tot_area);

% decide whether some polygons overlap
for i=1:length(c)
    gi = c{i};
    if(i==420)
        disp('420')
    end
    for j=i+1:length(c)
        gj = c{j};
        if(j==441)
            disp('441')
        end
        % combine the vertex entries
        combined =[gi gj];
        uniq = unique(combined);
        nc=length(combined);
        nu=length(uniq);
        if ( nu> nc )
            error('nu cannot be > than nc');
        end
        
        if ( nu < nc-1)
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
                        error('missing points ...');
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

figure(19);
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

%%%%%%%%%%%%%
output_file1=strcat('.\figs\z_mesh_poly_L',int2str(L),'_n',int2str(n),'_a',num2str(fraction,3));
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
if(~logi_save), return; end
%%%%%%%%%%%%%

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

nv=length(v(:,1))-length(inf_nan)-length(del);
fprintf(fid,'# grid vertices \n');
fprintf(fid,'%d\n',nv);
for k=1:length(v(:,1))
    if( ~isempty(find(inf_nan==k)) )
        fprintf('skipping inf_nan point %d \n',k);
        continue
    end
    if( ~isempty(find(del==k)) )
        fprintf('skipping omitted point %d \n',k);
        continue
    end
    fprintf(fid,'%g %g \n',v(k,1),v(k,2) );
end

fclose(fid);


return


