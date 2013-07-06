function create_vtk_output(filename,ndof,nel,connectivity,vert,z)

t1=cputime;

% vtk output for mesh
%                ====
str = sprintf('%s_mesh.vtk',filename);
fid = fopen(str,'w');
fprintf(fid,'# vtk DataFile Version 3.0 \n');

% date and time;
[y, m, d, h, mi, s] = datevec(now);
fprintf(fid,'Date: %d/%d/%d   Time: %d:%d\n', m, d, y, h, mi);

% the file format: ASCII or BINARY
fprintf(fid,'ASCII \n');
% dataset structure
fprintf(fid,'DATASET UNSTRUCTURED_GRID \n');
fprintf(fid,' \n');

% all vertices
fprintf(fid,'POINTS  %d  double \n',ndof); % nbr mesh vertices = ndof
n_triangles=0;
for iel=1:nel
    % get vertices
    g=connectivity(iel,:);
    v=vert(g,:);
    nv=length(g);
    n_triangles = n_triangles + nv;
    for k=1:nv
        fprintf(fid,'  %E %E %E \n',v(k,1),v(k,2),0.);
    end
end
fprintf(fid,' \n');

% all cells
fprintf(fid,'CELLS %d %d \n',nel,nel+ndof);
for iel=1:nel
    % get vertices
    g=connectivity(iel,:);
    nv=length(g);
    fprintf(fid,' %d ',nv); % number of vertices/polygon
    for k=1:nv
        fprintf(fid,'%d ',g(k)-1); % vertex ID, numbering starts at 0!
    end
    fprintf(fid,' \n');
end
fprintf(fid,' \n');

% cell types
fprintf(fid,'CELL_TYPES %d \n',nel);
for iel=1:nel
    fprintf(fid,' %d \n',7); % polygon type=7
end

fclose(fid);

%------------------------------------------------
% vtk output for solution
%                ========
str = sprintf('%s_solution.vtk',filename);
fid = fopen(str,'w');
fprintf(fid,'# vtk DataFile Version 3.0 \n');

% date and time;
fprintf(fid,'Date: %d/%d/%d   Time: %d:%d\n', m, d, y, h, mi);

% the file format: ASCII or BINARY
fprintf(fid,'ASCII \n');
% dataset structure
fprintf(fid,'DATASET UNSTRUCTURED_GRID \n');
fprintf(fid,' \n');

% all vertices (based on the triangular mesh)
fprintf(fid,'POINTS  %d  double \n',n_triangles*3);
for iel=1:nel
    % get vertices
    g=connectivity(iel,:);
    v=vert(g,:);
    % compute element's centroid
    vC=mean(v);
    
    nv=length(g);
    for k=1:nv
        kp1=k+1;
        if(kp1>nv), kp1=1; end
        fprintf(fid,'  %E %E %E \n',v(k,1)  ,v(k,2)  ,0.);
        fprintf(fid,'  %E %E %E \n',v(kp1,1),v(kp1,2),0.);
        fprintf(fid,'  %E %E %E \n',vC(1)   ,vC(2)   ,0.);
    end
end
fprintf(fid,' \n');

% all cells
fprintf(fid,'CELLS %d %d \n',n_triangles,n_triangles*(1+3)); 
for iel=1:n_triangles
    ind=(iel-1)*3;
    fprintf(fid,' %d ',3); % number of vertices/triangle
    fprintf(fid,'%d %d %d \n',ind,ind+1,ind+2); % vertex ID
end
fprintf(fid,' \n');

% cell types
fprintf(fid,'CELL_TYPES %d \n',n_triangles);
for iel=1:n_triangles
    fprintf(fid,' %d \n',5); % triangle type=5
end
fprintf(fid,' \n');

% add solution data
fprintf(fid,'POINT_DATA %d %d \n',n_triangles*3);
fprintf(fid,'SCALARS flux_pwld double 1 \n');
fprintf(fid,'LOOKUP_TABLE   default \n');
for iel=1:nel
    % get connectivity
    g=connectivity(iel,:);
    nv=length(g);
    % get the dofs
    local_dof =  z(g);
    % compute center value
    zC=mean(local_dof);
    % loop over sides
    for k=1:nv
        kp1=k+1;
        if(kp1>nv), kp1=1; end
        fprintf(fid,'  %E %E %E \n',local_dof(k),local_dof(kp1),zC);
    end
end

fclose(fid);

t2=cputime;
fprintf('Vtk time      = %g \n',t2-t1);



return
end