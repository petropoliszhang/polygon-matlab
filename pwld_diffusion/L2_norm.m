function L2_error=L2_norm(ndof,nel,connectivity,vert,n_quad,z,exact)

t1=cputime;

% init L2-error
L2_error = 0;

% loop over elements
for iel=1:nel

    % get connectivity and vertices
    g=connectivity{iel}(:);
    v=vert(g,:);

    % get the dofs
    local_dof =  z(g);

    % compute element's centroid
    vC=mean(v);
    zC=mean(local_dof);

    % size of array to store local integral
    nv=length(g);
    Q=zeros(nv,3);

    % loop over sides
    for iside=1:nv
        % pick 1st vertex
        irow1=iside;
        % pick 2ndt vertex
        irow2=irow1+1; if(irow2>nv), irow2=1; end
        % assign A and B
        vA=v(irow1,:); vB=v(irow2,:);
        % create triangle vertex list
        triangle_vert=[vA; vB; vC];
        % get quadrature on that triangle
        [X,Y,Wx,Wy]=triquad(n_quad,triangle_vert);
        % create the 3 basis functions
        tf1=@(x,y) ( vC(1)*(y-vB(2)) + x*(vB(2)-vC(2)) + vB(1)*(-y+vC(2)) ) / ...
            ( vC(1)*(vA(2)-vB(2)) + vA(1)*(vB(2)-vC(2)) + vB(1)*(-vA(2)+vC(2)) );
        tf2=@(x,y) ( vC(1)*(-y+vA(2)) + vA(1)*(y-vC(2)) + x*(-vA(2)+vC(2)) ) / ...
            ( vC(1)*(vA(2)-vB(2)) + vA(1)*(vB(2)-vC(2)) + vB(1)*(-vA(2)+vC(2)) );
        tf3=@(x,y) ( vB(1)*(y-vA(2)) + x*(vA(2)-vB(2)) + vA(1)*(-y+vB(2)) ) / ...
            ( vC(1)*(vA(2)-vB(2)) + vA(1)*(vB(2)-vC(2)) + vB(1)*(-vA(2)+vC(2)) );

        % compute the numerical solution on each side
        local_num_solu = local_dof(irow1)*feval(tf1,X,Y) + ...
            local_dof(irow2)*feval(tf2,X,Y) + ...
            zC              *feval(tf3,X,Y) ;
        % compute the EXACT solution on each side
        local_exact = feval(exact,X,Y);

        % compute the contribution of that side to the L2-error

        % evaluate each integral per triangle for the 3 basis functions
        L2_error = L2_error + Wx'*(local_num_solu-local_exact).^2*Wy;

    end % end loop over sides
end % end loop over elements

t2=cputime;
fprintf('L2 norm time  = %g \n',t2-t1);

fprintf('\t ndof = %d, L2_error = %+E \n\n',ndof,L2_error);

return
end