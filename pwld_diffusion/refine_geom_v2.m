function [nel,ndof,connectivity,corner_pos,vert,n_edge,edg2poly,edg2vert,edg_perp,i_mat,i_src,curr_ref_lev] =...
    refine_geom_v2(err_i,frac_ref,curr_ref_lev,old_nel,old_connectivity,old_corner_pos,old_vert,...
    old_edg2poly,old_edg2vert,old_i_mat,old_i_src)
% mesh refinement
insert_fn = @(val, x, pos)cat(2,  x(1:pos-1), val, x(pos:end));

t1=cputime;

global verbose

% make a copy of the current refinement level
next_ref_lev = curr_ref_lev;

% cells to refine
err_max=max(err_i);
err_i=err_i/err_max;
if(frac_ref<0 || frac_ref>1), error('frac ref must be between 0 and 1'); end
cell_refine = find( err_i >= frac_ref );

% increment the refinement level array
next_ref_lev(cell_refine) = next_ref_lev(cell_refine) + 1;
% implement a 1-level refinement difference rule
[next_ref_lev,cell_refine]=ref_difference_rule(next_ref_lev,cell_refine,...
    old_edg2poly,old_connectivity);
% sanity check
dmax = max(next_ref_lev - curr_ref_lev);
if(dmax>1)
    [next_ref_lev curr_ref_lev]
    error('max diff in ref levels >1');
end

% number of cells to refine
n_cell_to_ref = length(cell_refine);

%----------------------------------------
% the new mesh is created by generating the cells with the lowest
% refinement levels first

all_cells = 1:old_nel;
unchanged_cells = setxor(all_cells,cell_refine);
n_unchanged = length(unchanged_cells);
%sanity check
aux = intersect(all_cells,cell_refine);
aux = aux -sort(cell_refine);
if(any(aux))
    error('intersect(all_cells,cell_refine)');
end
%sanity check
if( (n_unchanged+n_cell_to_ref)~=old_nel)
    error('(n_unchanged+n_cell_to_ref)~=old_nel');
end

% refinement levels for the cells to be refined
next_ref_lev_to_process = next_ref_lev(cell_refine);
[sorted_next_ref_lev_to_process,ordered_cell_refine]=sort(next_ref_lev_to_process);
sorted_curr_ref_lev_to_process = curr_ref_lev(ordered_cell_refine);


level_values = unique(next_ref_lev_to_process);
how_many_levels=length(level_values);
ncells_per_lev=zeros(how_many_levels,1);
for k=1:how_many_levels
    ind = find( next_ref_lev == level_values(k) );
    ncells_per_lev(k) = length(ind);
end

%sanity check
if( length(ordered_cell_refine) ~= length(cell_refine) )
    error('length(ordered_cell_refine) ~= length(cell_refine)');
end

%----------------------------------------

% new nel
nel = (old_nel-n_cell_to_ref) + n_cell_to_ref*4;

% new connectivity, imat, isrc
connectivity=cell(nel,1);
corner_pos=cell(nel,1);
i_mat=zeros(nel,1);
i_src=zeros(nel,1);
n_vertices=zeros(nel,1);
curr_ref_lev=zeros(nel,1);

vert=zeros(4*nel,2); % 4*nel is only a lower bound ...

grid_vert_ID=[];
new_iel = 0;
new_dof = 0;


%%% cells to be refined
for k=1:n_cell_to_ref
    
    % get old cell ID
    iel = ordered_cell_refine(k);
    corners = old_corner_pos{iel};
    g = old_connectivity{iel}(corners);
    nedg = length(g);
    if(nedg~=4), error('refinement implemented for quads only as of now'); end
    v = old_vert(g,:);
    
    % get refinement levels of the old cell
    next_lev = sorted_next_ref_lev_to_process(k);
    curr_lev = sorted_curr_ref_lev_to_process(k);
    % sanity checks
    if(next_lev ~= next_ref_lev(iel)), error('next_lev ~= next_ref_lev(iel)'); end
    if(curr_lev ~= curr_ref_lev(iel)), error('curr_lev ~= curr_ref_lev(iel)'); end
    if(next_lev == curr_lev), error('next_lev == curr_lev'); end
    
    %%%%%%%% refine old cell first
    
    % new points
    centroid = mean(v);
    i=1; new_pt1 = (v(i,:)+v(i+1,:))/2;
    i=2; new_pt2 = (v(i,:)+v(i+1,:))/2;
    i=3; new_pt3 = (v(i,:)+v(i+1,:))/2;
    i=4; new_pt4 = (v(i,:)+v(1,:)  )/2;
    
    %  d---3---c
    %  |       |
    %  4   c   2
    %  |       |
    %  a---1---b
    %
    % ---- quad 1
    new_iel=new_iel+1;
    connectivity{new_iel}=zeros(nedg,1);
    i_mat(new_iel) = old_i_mat(iel);
    i_src(new_iel) = old_i_src(iel);
    curr_ref_lev(new_iel) = next_ref_lev(iel);
    
    local_con = (new_dof+1):(new_dof+nedg);
    connectivity{new_iel}(:) = local_con;
    corner_pos{new_iel} = 1:4;
    n_vertices(new_iel) = nedg;
    
    vert(new_dof+1,:) = v(1,:);
    vert(new_dof+2,:) = new_pt1;
    vert(new_dof+3,:) = centroid;
    vert(new_dof+4,:) = new_pt4;
    new_dof = new_dof + nedg;
    
    grid_vert_ID = [ grid_vert_ID local_con];
    
    %
    % ---- quad 2
    new_iel=new_iel+1;
    connectivity{new_iel}=zeros(nedg,1);
    i_mat(new_iel) = old_i_mat(iel);
    i_src(new_iel) = old_i_src(iel);
    curr_ref_lev(new_iel) = next_ref_lev(iel);
    
    local_con = (new_dof+1):(new_dof+nedg);
    connectivity{new_iel}(:) = local_con;
    corner_pos{new_iel} = 1:4;
    n_vertices(new_iel) = nedg;
    
    vert(new_dof+1,:) = new_pt1;
    vert(new_dof+2,:) = v(2,:);
    vert(new_dof+3,:) = new_pt2;
    vert(new_dof+4,:) = centroid;
    new_dof = new_dof + nedg;
    
    grid_vert_ID = [ grid_vert_ID local_con];
    
    %
    % ---- quad 3
    new_iel=new_iel+1;
    connectivity{new_iel}=zeros(nedg,1);
    i_mat(new_iel) = old_i_mat(iel);
    i_src(new_iel) = old_i_src(iel);
    curr_ref_lev(new_iel) = next_ref_lev(iel);
    
    local_con = (new_dof+1):(new_dof+nedg);
    connectivity{new_iel}(:) = local_con;
    corner_pos{new_iel} = 1:4;
    n_vertices(new_iel) = nedg;
    
    vert(new_dof+1,:) = centroid;
    vert(new_dof+2,:) = new_pt2;
    vert(new_dof+3,:) = v(3,:);
    vert(new_dof+4,:) = new_pt3;
    new_dof = new_dof + nedg;
    
    grid_vert_ID = [ grid_vert_ID local_con];
    
    %
    % ---- quad 4
    new_iel=new_iel+1;
    connectivity{new_iel}=zeros(nedg,1);
    i_mat(new_iel) = old_i_mat(iel);
    i_src(new_iel) = old_i_src(iel);
    curr_ref_lev(new_iel) = next_ref_lev(iel);
    
    local_con = (new_dof+1):(new_dof+nedg);
    connectivity{new_iel}(:) = local_con;
    corner_pos{new_iel} = 1:4;
    n_vertices(new_iel) = nedg;
    
    vert(new_dof+1,:) = new_pt4;
    vert(new_dof+2,:) = centroid;
    vert(new_dof+3,:) = new_pt3;
    vert(new_dof+4,:) = v(4,:);
    new_dof = new_dof + nedg;
    
    grid_vert_ID = [ grid_vert_ID local_con];
    
    
    %%%%%%%% then analyze what to do about the neighbors of that old cell
    
    % find the old neighbors to old iel
    list_edge_p = find( old_edg2poly(:,2) == iel );
    Kp = old_edg2poly(list_edge_p,2);
    % remove iel from list
    ind = find(Kp==iel); Kp(ind)=[];
    % remove boundary connections
    ind = find(Kp<0);
    ind = sort(ind,'descend');
    for ii=1:length(ind)
        Kp(ind(ii))=[];
    end
    list_edge_m = find( old_edg2poly(:,1) == iel );
    Km = old_edg2poly(list_edge_p,1);
    % remove iel from list
    ind = find(Km==iel); Km(ind)=[];
    % remove boundary connections
    ind = find(Km<0);
    ind = sort(ind,'descend');
    for ii=1:length(ind)
        Km(ind(ii))=[];
    end
    % generate lists of neighboring polygons and edges
    K_others=[Kp Km];
    list_edges = [list_edge_p list_edge_m];
    
    %  analyze each neighbor
    for kk=1:length(K_others);
        i_neigh = K_others(kk);
        next_lev_neigh = next_ref_lev(i_neigh);
        curr_lev_neigh = curr_ref_lev(i_neigh);
        
        if ( next_lev == next_lev_neigh)
            %%% do nothing more %%%
        elseif( next_lev - next_lev_neigh == 1)
            if( curr_lev == curr_lev_neigh )
                %%% add midpoint %%%
                % find edge
                [i_,j_]=find( old_edg2poly(list_edges,:) == i_neigh );
                if(length(i_)~=1 || length(j_)~=1)
                    error('i_neigh multiplicity in list_edges');
                end
                common_edge = list_edges(i_);
                if(j_2==1)
                    ed = old_edg2vert(common_edge,1:2);
                else
                    ed = old_edg2vert(common_edge,3:4);
                end
                v_neigh_ed = old_vert(ed,:);
                mid_pt = mean(v_neigh_ed);
                
                % add midpt to data from old neigh
                g_neigh = old_connectivity{i_neigh};
                i2 = find( g_neigh == ed(1) );
                if(length(i2)~=1), error('g_neigh'); end
                n_old_vert = length(old_vert(:,1));
                g_neigh = insert_fn(n_old_vert+1, g_neigh, i2+1);
                % deal with corners
                neigh_corners = old_corner_pos{i_neigh};
                c_neigh = old_connectivity{i_neigh}(neigh_corners);
                i2 = find( c_neigh == ed(1) );
                if(length(i2)~=1), error('c_neigh'); end
                neigh_corners(i2+1:end) = neigh_corners(i2+1:end) +1;
                % update struct
                old_connectivity{i_neigh} = g_neigh;
                old_vert(n_old_vert+1,:) = mid_pt;
                old_corner_pos{i_neigh} = neigh_corners;
                %             elseif( curr_lev - curr_lev_neigh == 1)
                %                 % second part
            else
                error(' wrong case for:  next_lev - next_lev_neigh == 1');
            end
        elseif( next_lev - next_lev_neigh == -1)
            if( curr_lev - curr_lev_neigh == -1 )
                %%% first part %%%
                % find edge
                [i_,j_]=find( old_edg2poly(list_edges,:) == i_neigh );
                if(length(i_)~=1 || length(j_)~=1)
                    error('i_neigh multiplicity in list_edges');
                end
                common_edge = list_edges(i_);
                if(j_2==1)
                    ed = old_edg2vert(common_edge,3:4); % flipped
                else
                    ed = old_edg2vert(common_edge,1:2); % flipped
                end
                v_iel_ed = old_vert(ed,:);
                mid_pt = mean(v_iel_ed);
                
                % which of the 4 quads that replace old iel should be
                % selected ?
                for ii=1:4
                    g = connectivity{new_iel-ii+1};
                    % not ed(1), we want a unique vertex, thus the corner of the new quads
                    i2 = find( g == ed(2) );
                    % exit small loop when found
                    if(~isempty(i2))
                        break;
                    end
                end
                if(length(i2)~=1), error('g new quads'); end
                new_iel_to_modify = new_iel-ii+1;
                
                % add midpt to data from newly created quads
                new_dof = new_dof + 1;
                vert(new_dof,:) = mid_pt;
                n_vertices(new_iel_to_modify) = n_vertices(new_iel_to_modify) + 1;
                corner_pos{new_iel_to_modify}(i2:end) =corner_pos{new_iel_to_modify}(i2:end) + 1;
                
                g = insert_fn(new_dof, g, i2);
                connectivity{new_iel_to_modify} = g;
                grid_vert_ID = [ grid_vert_ID new_dof];
                
            else
                error(' wrong case for:  next_lev - next_lev_neigh == -1');
            end
        else
            error('|next_lev - next_lev_neigh| >1')
        end
        
    end % end loop on neighbors
    
end

%%% unchanged cells
for k=1:n_unchanged
    
    % get old cell ID
    iel = unchanged_cells(k);
    g = old_connectivity{iel};
    nedg = length(g);
    v = old_vert(g,:);
    
    % ---- create simple quad
    new_iel=new_iel+1;
    connectivity{new_iel}=zeros(nedg,1);
    i_mat(new_iel) = old_i_mat(iel);
    i_src(new_iel) = old_i_src(iel);
    curr_ref_lev(new_iel) = next_ref_lev(iel);
    
    local_con = (new_dof+1):(new_dof+nedg);
    connectivity{new_iel} = local_con;
    corner_pos{new_iel} = old_corner_pos{iel};
    n_vertices(new_iel) = nedg;
    
    vert(local_con,1:2) = v(:,:);
    new_dof = new_dof + nedg;
    
    grid_vert_ID = [ grid_vert_ID local_con];
    
end

% sanity checks
test = (new_iel == nel);
if(~test), error('new_iel ~= nel'); end

% ndof
ndof=new_dof;
test = (ndof == length(vert(:,1)));
if(~test), error('new_dof'); end
test = (ndof == sum(n_vertices));
if(~test), error('new_dof n_vertices'); end

min_vert=min(n_vertices);
max_vert=max(n_vertices);
for k=min_vert:max_vert
    nnn=find(n_vertices==k);
    fprintf('Number of polygons with %d vertices = %d \n',k,length(nnn));
end
fprintf('Total number of polygons           = %d \n',nel);

% deal with grid_vert_ID
grid_vert_ID = unique(grid_vert_ID);
n_grid_vert = length(grid_vert_ID);
vert_grid = zeros(n_grid_vert,1);
vert_grid = vert(grid_vert_ID,:);

% find duplicated vertices in vert_grid
del=[];
for k=1:n_grid_vert
    vk=vert_grid(k,1:2);
    % get difference
    aux = vert_grid - kron(ones(n_grid_vert,1),vk);
    % get vector of norm
    aa=sqrt(aux(:,1).^2+aux(:,2).^2);
    ind=find(aa<1e-12);
    len=length(ind);
    if(len==0),
        vk
        error('vk not found ????');
    end
    if(len >1),
        %         warning('vk duplicated');
        vk;
        vert_grid(ind,:);
        k;
        ind;
        ind2=find(ind>k);
        ind(ind2);
        del=[del; ind(ind2)];
    end
end
% delete in descending order
del=unique(del);
[del,isort]=sort(del,'descend');
for i=1:length(del)
    if(verbose)
        if(i==1)
            fprintf('the following point will be omitted b/c they are DUPLICATES \n');
        end
        fprintf(' vertex # %d, x-value=%g, y-value=%g \n',del(i),vert_grid(del(i),1),vert_grid(del(i),2) );
    end
    vert_grid(del(i),:)=[];
end
n_grid_vert=length(vert_grid(:,1));

% complete mesh data
[n_edge,edg2poly,edg2vert,edg_perp ] = complete_mesh_data( nel,ndof,n_grid_vert,vert,vert_grid,connectivity );

t2=cputime;
fprintf('Mesh time     = %g \n\n',t2-t1);

return
end
