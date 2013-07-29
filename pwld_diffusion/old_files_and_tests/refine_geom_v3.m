function [nel,ndof,connectivity,vert,n_edge,edg2poly,edg2vert,edg_perp,i_mat,i_src,curr_ref_lev] =...
    refine_geom_v3(err_i,frac_ref,curr_ref_lev,old_nel,old_connectivity,old_vert,old_i_mat,old_i_src,old_edg2poly)
% mesh refinement

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

%----------------------------------------
% the new mesh is created by generating the cells with the lowest 
% refinement levels first

% all_cells = 1:old_nel;
% unchanged_cells = setxor(all_cells,cell_refine);
% n_unchanged = length(unchanged_cells);
% %sanity check
% aux = intersect(all_cells,cell_refine);
% aux = aux -sort(cell_refine);
% if(any(aux))
%     error('intersect(all_cells,cell_refine)');
% end
% %sanity check
% if( (n_unchanged+n_cell_to_ref)~=old_nel)
%     error('(n_unchanged+n_cell_to_ref)~=old_nel');
% end

% refinement levels 
[sorted_next_ref_lev,ordered_cells]=sort(next_ref_lev);
sorted_curr_ref_lev = curr_ref_lev(ordered_cells);

level_values = unique(next_ref_lev);
how_many_levels=length(level_values);
ncells_per_lev=zeros(how_many_levels,1);
for k=1:how_many_levels
    ind = find( next_ref_lev == level_values(k) );
    ncells_per_lev(k) = length(ind);
end


%----------------------------------------

% new nel
nel = (old_nel-n_cell_to_ref) + n_cell_to_ref*4; 

% new connectivity, imat, isrc
connectivity=cell(nel,1);
corners=cell(nel,1);
i_mat=zeros(nel,1);
i_src=zeros(nel,1);
n_vertices=zeros(nel,1);
curr_ref_lev=zeros(nel,1);

vert=zeros(4*nel,2); % 4*nel is only a lower bound ...

grid_vert_ID=[];
new_iel = 0;
new_dof = 0;

%%% loop over cells
for k=1:old_nel
    
    % get old cell ID
    iel = ordered_cells(k);
    g = old_corners{iel};
    nedg = length(g);
    if(nedg~=4), error('refinement implemented for quads only as of now'); end
    v = old_vert(g,:);
    
    % get refinement levels of the old cell 
    next_lev = sorted_next_ref_lev(k);
    curr_lev = sorted_curr_ref_lev(k);
    % sanity check
    if(next_lev ~= next_ref_lev(iel)), error('next_lev ~= next_ref_lev(iel)'); end
    
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
    n_vertices(new_iel) = nedg;

    vert(new_dof+1,:) = v(1,:);
    vert(new_dof+2,:) = new_pt1;
    vert(new_dof+3,:) = centroid;
    vert(new_dof+4,:) = new_pt4;
    new_dof = new_dof + nedg;

    
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
    n_vertices(new_iel) = nedg;

    vert(new_dof+1,:) = new_pt4;
    vert(new_dof+2,:) = centroid;
    vert(new_dof+3,:) = new_pt3;
    vert(new_dof+4,:) = v(4,:);
    new_dof = new_dof + nedg;
    
    grid_vert_ID = [ grid_vert_ID local_con];


    %%%%%%%%%%%%%%%%%%%%%
    % find the 4 (quads only!!!) old neighbors to old iel
    list_edge_p = find( old_edg2poly(:,2) == iel );
    Kp = old_edg2poly(list_edge_p,2);
    ind = find(Kp==iel); Kp(ind)=[];
    ind = find(Kp<0); 
    ind = sort(ind,'descend');
    for ii=1:length(ind)
        Kp(ind(ii))=[];
    end   
    list_edge_m = find( old_edg2poly(:,1) == iel );
    Km = old_edg2poly(list_edge_p,1);
    ind = find(Km==iel); Km(ind)=[];
    ind = find(Km<0); 
    ind = sort(ind,'descend');
    for ii=1:length(ind)
        Km(ind(ii))=[];
    end
    K_others=[Kp Km];
    
    % ------ below edge 1-2
    ed=g(1:2);
    [ old_connectivity,old_vert ] = insert_pt_quad( ed,K_others,old_connectivity,old_vert,cell_refine,new_pt1 );
    % ------ right of edge 2-3
    ed=g(2:3);
    [ old_connectivity,old_vert ] = insert_pt_quad( ed,K_others,old_connectivity,old_vert,cell_refine,new_pt2 );
    % ------ top of edge 3-4
    ed=g(3:4);
    [ old_connectivity,old_vert ] = insert_pt_quad( ed,K_others,old_connectivity,old_vert,cell_refine,new_pt3 );
    % ------ left of edge 4-1
    ed=[g(4); g(1)];
    [ old_connectivity,old_vert ] = insert_pt_quad( ed,K_others,old_connectivity,old_vert,cell_refine,new_pt4 );


end
 
% sanity check
test = (new_iel == 4*n_cell_to_ref);
if(~test), error('new_iel ~= 4*n_cell_to_ref'); end

% add the un-refined cells to the new data
for iel=1:old_nel
    ind = find (cell_refine ==iel);
    if(~isempty(ind))
        % that old cell has just been refined
        if(length(ind)~=1), error('old cell refined more than once ????'); end
        continue
    end
    
    % cell has not been refined. add it to the new list
    g = old_connectivity{iel};
    nedg = length(g); % not necessarily =4 because of AMR
    v = old_vert(g,:);

    new_iel=new_iel+1;
    connectivity{new_iel}=zeros(nedg,1);
    i_mat(new_iel) = old_i_mat(iel);
    i_src(new_iel) = old_i_src(iel);
    
    local_con = (new_dof+1):(new_dof+nedg);
    connectivity{new_iel} = local_con;
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