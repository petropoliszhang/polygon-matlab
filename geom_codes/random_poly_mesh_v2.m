% clear all;
close all; clc;

logi_save = true;

L=10;
n=20;
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
% x=[            0
%     0
%     0
%     6.2601
%     1.5599
%     6.0209
%     10
%     10
%     10];
% y=[            0
%     9.2129
%     10
%     0
%     8.9765
%     10
%     0
%     5.3712
%     10];
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
    % re-order poly ccw
    g=c{i}(:);
    vv=v(g,:);
    [x2, y2] = poly2ccw(vv(:,1), vv(:,2));
    % re-order connectivity so that the poly is actually ccw
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

clc
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
del = [];
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
        flag=0;
        if ( nu < nc-2)
            % need to investigate
            %reverse order for gj
            n_gj=length(gj);
            for k=1:n_gj
                gjr(k)=gj(n_gj+1-k);
            end
            % find the intersection, but the result unfortunetaly comes out sorted
            [inter,ia,ib]=intersect(gi,gjr);
            n_gi=length(gi);
            % determine the sweeping order (postive or negastive) of the common elements
            val=inter(1);
            posi=ia(1);
            posj=ib(1);
            nexti=posi+1; if(nexti>n_gi), nexti=1; end
            previ=posi-1; if(previ<1), previ=n_gi; end
            % sign should be the same if the common vertices belong to
            % non-overlapping polygons
            sign=0;
            if(ia(2)==nexti), sign=+1; end
            if(ia(2)==previ), sign=-1; end
            if(sign==0)
                error('sign i =0');
            end
%             fprintf('sign = %d \n',sign);
            % now, find the position of the first common element in gi
            [dum,posi]=min(ia);
%             fprintf('posi = %d \n\n',posi);
%             gi
%             gjr
            indj=ib(posi);
            for k=1:length(ia)
                indi=ia(posi);
                vali=gi(indi);
%                 fprintf('indi = %d vali = %g\n',indi,vali);
                valj=gjr(indj);
%                 fprintf('indj = %d valj = %g\n\n',indj,valj);
                if(abs(vali-valj)>eps)
                    flag=1;
%                     warning(' aaa ' );
                end
                posi=posi+1;if(posi>length(ia)), posi=1; end
                indj=indj+sign;
                if(indj>n_gj), indj=1; end
                if(indj<1), indj=n_gi; end
            end
        end

        if(flag==1)
            % shortest poly scheduled for deletion
            if(n_gi<n_gj)
                del = [ del i ];
            else
                del = [ del j ];
            end
        end
    end
end
del
unique(del)
length(c)
for k=1:length(del)
    c{del(k)}=[];
end
% either one is fine 
% result = c(~cellfun('isempty',c))
% cc(cellfun(@isempty,cc)) = []
CC = c(~cellfun('isempty',c));
length(CC)

for id=1:length(CC)
    figure(10)
    hold all
    vv= v(CC{id},:);
    xc=mean(vv(:,1));
    yc=mean(vv(:,2));
    vv(end+1,:)=vv(1,:);
    plot(vv(:,1),vv(:,2),'+-')
    plot(xc,yc,'x');
    %     pause
    figure(11)
    patch(v(CC{id},1),v(CC{id},2),id); % use color i.
end


tot_area=0;
for iel=1:length(CC)
    g=CC{iel}(:);
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

