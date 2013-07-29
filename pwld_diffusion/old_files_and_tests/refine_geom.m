function [nel,ndof,connectivity,vert,n_edge,edg2poly,edg2vert,edg_perp,i_mat,i_src] =...
    refine_geom(err_i,frac_ref,old_nel,old_connectivity,old_vert,old_i_mat,old_i_src,old_edg2poly)
% mesh refinement

t1=cputime;

global verbose

% cells to refine
err_max=max(err_i);
err_i=err_i/err_max;
if(frac_ref<0 || frac_ref>1), error('frac ref must be between 0 and 1'); end
cell_refine = find( err_i > frac_ref );

% number of cells to refine
n_cell_to_ref = length(cell_refine);

% new nel
nel = (old_nel-n_cell_to_ref) + n_cell_to_ref*4; 

% new connectivity, imat, isrc
connectivity=cell(nel,1);
i_mat=zeros(nel,1);
i_src=zeros(nel,1);
n_vertices=zeros(nel,1);

vert=zeros(4*nel,2); % 4*nel is only a lower bound ...

grid_vert_ID=[];
new_iel = 0;
new_dof = 0;
for k=1:n_cell_to_ref
    
    % get old cell ID
    iel = cell_refine(k);
    g = old_connectivity{iel};
    nedg = length(g);
    if(nedg~=4), error('refinement implemented for quads only as of now'); end
    v = old_vert(g,:);
   
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
