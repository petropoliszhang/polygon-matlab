function z = DG_assemble_solve( ndof,nel,n_edge,vert,connectivity,edg2poly,edg2vert,edg_normal,C_pen,C_pen_bd,...
                                i_mat,c_diff,sigma_a,i_src,S_ext,logi_mms,mms,n_quad,bc_type,bc_val )


t1=cputime;

%------------------------------------------------
% DG assemble volumetric terms
A = spalloc(ndof,ndof,18*ndof); b=zeros(ndof,1);
for iel=1:nel
    g=connectivity{iel}(:);
    v=vert(g,:);
    mat = i_mat(iel);
    [M,K,f,grad{iel}]=build_pwld_local_matrices(g,v);
    A(g(:),g(:)) = A(g(:),g(:)) + c_diff(mat)*K +sigma_a(mat)*M;
    % rhs contribution
    if(logi_mms)
        % centroid
        vC=mean(v);
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
            % evaluate each integral per triangle for the 3 basis functions
            Q(iside,1)=Wx'*(feval(mms,X,Y).*feval(tf1,X,Y))*Wy;
            Q(iside,2)=Wx'*(feval(mms,X,Y).*feval(tf2,X,Y))*Wy;
            Q(iside,3)=Wx'*(feval(mms,X,Y).*feval(tf3,X,Y))*Wy;
        end
        % compute contribution to global rhs
        local_rhs=zeros(nv,1);
        alpha=1/nv;
        common = alpha*sum(Q(:,3));
        for iside=1:nv
            iside_m1 = iside-1;
            if(iside_m1==0), iside_m1=nv; end
            local_rhs(iside) = Q(iside,1) + Q(iside_m1,2) + common;
        end
        b(g(:)) = b(g(:)) + local_rhs;
    else
        b(g(:)) = b(g(:)) + S_ext(mat)*f;
    end
end

%------------------------------------------------
% DG assemble edge terms
%
%           v2 ^  w1
%              |
%              |  n_ed
%  Minus side  | --->    Plus side
%              |
%              |
%           v1 .  w2
%
m1d=[2 1 ; 1 2]/6;
m1d_mod=[1 2; 2 1]/6;

for ied=1:n_edge
    % get K-,K+ and their connectivities
    Kp = edg2poly(ied,2);
    Km = edg2poly(ied,1);
    % we want to loop only on INTERIOR edges
    if(Km<=0 || Kp==0), error('Km<=0 or Kp==0'); end
    if(Kp<0), continue; end
    % get the polygons' connectivities
    gp = connectivity{Kp}(:);
    gm = connectivity{Km}(:);
    % nbr of vertices in each poly
    nvp = length(gp);
    nvm = length(gm);
    % get normal
    ne = edg_normal(ied,1:2);
    % get edge vertices of K+ side
    W = edg2vert(ied,3:4);
    % get edge vertices of K- side
    V = edg2vert(ied,1:2);
    % get the local ID of the edge's first vertex in K+
    indp = find( gp == W(1) );
    if(length(indp) ~= 1), error('length(IDp) ~= 1'); end
    if(indp ~= nvp)
        IDp = [indp (indp+1) ];
    else
        IDp = [indp 1 ];
    end
    % get the local ID of the edge's first vertex in K-
    indm = find( gm == V(1) );
    if(length(indm) ~= 1), error('length(IDm) ~= 1'); end
    if(indm ~= nvm)
        IDm = [indm (indm+1) ];
    else
        IDm = [indm 1 ];
    end
    %     % skipping indices
    %     skip_p = [(indp:nvp) (1:indp-1)];
    %     skip_m = [(indm:nvm) (1:indm-1)];
    %     skip_pr= [(indp:-1:1) (nvp:-1:indp+1)];
    %     skip_mr= [(indm:-1:1) (nvm:-1:indm+1)];

    % length current edge
    Le = norm( diff(vert(V,:)) );
    % material properties
    Dp = c_diff(i_mat(Kp));
    Dm = c_diff(i_mat(Km));
    % penalty term
    h_perp=Le; % temporary!
    pen = C_pen * (Dp/h_perp + Dm/h_perp) /2;

    % build the local edge gradient matrices
    % [[phi]],{{D.grad(b).ne}}
    %      = (phi+ - phi-)(D+ grad(b+).ne + D- grad(b-).ne)/2
    %      =  (phi+) D+ grad(b+).ne/2
    %       + (phi+) D- grad(b-).ne/2
    %       - (phi-) D+ grad(b+).ne/2
    %       - (phi-) D- grad(b-).ne/2
    % [[b]],{{D.grad(phi).ne}}
    %      = (b+ - b-)(D+ grad(phi+).ne + D- grad(phi-).ne)/2
    %      =  (b+) D+ grad(phi+).ne/2
    %       + (b+) D- grad(phi-).ne/2
    %       - (b-) D+ grad(phi+).ne/2
    %       - (b-) D- grad(phi-).ne/2
    % compute row vector: n' * G (ne is already stored as a 1x2 row vector)

    % -/-
    row_grad_m = ne * grad{Km}(:,:,indm);
    % edge matrix for this side
    col_b_m = zeros(nvm,1); col_b_m(IDm) = Le/2;
    aux = -Dm/2 * (col_b_m * row_grad_m + row_grad_m' * col_b_m');
    aux(IDm,IDm) = aux(IDm,IDm) + pen * Le * m1d;
    A(gm(:),gm(:)) = A(gm(:),gm(:)) + aux;

    % +/+
    row_grad_p = ne * grad{Kp}(:,:,indp);
    % edge matrix for this side
    col_b_p = zeros(nvp,1); col_b_p(IDp) = Le/2;
    aux = +Dp/2 * (col_b_p * row_grad_p + row_grad_p' * col_b_p');
    aux(IDp,IDp) = aux(IDp,IDp) + pen * Le * m1d;
    A(gp(:),gp(:)) = A(gp(:),gp(:)) + aux;

    % -(test)/+(solution)
    aux = ( -col_b_m * Dp/2*row_grad_p + Dm/2*row_grad_m' * col_b_p');
    aux(IDm,IDp) = aux(IDm,IDp) - pen * Le * m1d_mod;
    A(gm(:),gp(:)) = A(gm(:),gp(:)) + aux;

    % +(test)/-(solution)
    aux = ( col_b_p * Dm/2*row_grad_m - Dp/2*row_grad_p' * col_b_m');
    aux(IDp,IDm) = aux(IDp,IDm) - pen * Le * m1d_mod;
    A(gp(:),gm(:)) = A(gp(:),gm(:)) + aux;

end

%------------------------------------------------
% boundary conditions
for ied=1:n_edge
    % get K-,K+ for that edge
    Kp = edg2poly(ied,2);
    Km = edg2poly(ied,1);

    % we want to loop only on BOUNDARY edges
    if(Kp>0), continue; end

    % homogeneous Neumann on the left:
    if(Kp==-10 && bc_type(1)==2), continue; end
    % homogeneous Neumann on the right:
    if(Kp==-20 && bc_type(2)==2), continue; end
    % homogeneous Neumann on the bottom:
    if(Kp==-30 && bc_type(3)==2), continue; end
    % homogeneous Neumann on the left:
    if(Kp==-40 && bc_type(4)==2), continue; end

    %     [ied Kp Km]
    % get the polygons' connectivities
    gm = connectivity{Km}(:);
    % nbr of vertices in each poly
    nvm = length(gm);
    % get normal
    ne = edg_normal(ied,1:2);
    % get edge vertices of K- side
    V = edg2vert(ied,1:2);
    % get the local ID of the edge's first vertex in K-
    indm = find( gm == V(1) );
    if(length(indm) ~= 1), error('length(IDm) ~= 1'); end
    if(indm ~= nvm)
        IDm = [indm (indm+1) ];
    else
        IDm = [indm 1 ];
    end
    % length current edge
    Le = norm( diff(vert(V,:)) );
    % material properties
    Dm = c_diff(i_mat(Km));
    % penalty term
    h_perp=Le; % temporary!
    pen = C_pen_bd * Dm/h_perp;


    % in/homogeneous Dirichlet
    if(     (Kp==-10 && (bc_type(1)==0||bc_type(1)==1)) || ...
            (Kp==-20 && (bc_type(2)==0||bc_type(2)==1)) || ...
            (Kp==-30 && (bc_type(3)==0||bc_type(3)==1)) || ...
            (Kp==-40 && (bc_type(4)==0||bc_type(4)==1)) )

        row_grad_m = ne * grad{Km}(:,:,indm);
        % edge matrix for this side
        col_b_m = zeros(nvm,1); col_b_m(IDm) = Le/2;
        aux = - Dm * (col_b_m * row_grad_m + row_grad_m' * col_b_m'); % do not divide by 2 here!!!!
        aux(IDm,IDm) = aux(IDm,IDm) + pen * Le * m1d;
        A(gm(:),gm(:)) = A(gm(:),gm(:)) + aux;
    end

    % inhomogeneous Dirichlet
    if(     (Kp==-10 && bc_type(1)==1) || ...
            (Kp==-20 && bc_type(2)==1) || ...
            (Kp==-30 && bc_type(3)==1) || ...
            (Kp==-40 && bc_type(4)==1) )

        switch(Kp)
            case{-10}
                val=bc_val.left;
                gm_=[gm(1) gm(4)];
            case{-20}
                val=bc_val.right;
                gm_=[gm(2) gm(3)];
            case{-30}
                val=bc_val.bottom;
                gm_=[gm(1) gm(2)];
            case{-40}
                val=bc_val.top;
                gm_=[gm(3) gm(4)];
            otherwise
                error('inhomogeneous Neumann');
        end
        b(gm_(:)) = b(gm_(:)) + pen*val*Le/2;

        row_grad_m = ne * grad{Km}(:,:,indm);
        % edge matrix for this side
        aux = - Dm * (val*Le) * row_grad_m' ;
        b(gm(:)) = b(gm(:)) + aux;
    end

    % inhomogeneous Neumann
    if(     (Kp==-10 && bc_type(1)==3) || ...
            (Kp==-20 && bc_type(2)==3) || ...
            (Kp==-30 && bc_type(3)==3) || ...
            (Kp==-40 && bc_type(4)==3) )

        switch(Kp)
            case{-10}
                val=bc_val.left;
                gm_=[gm(1) gm(4)];
            case{-20}
                val=bc_val.right;
                gm_=[gm(2) gm(3)];
            case{-30}
                val=bc_val.bottom;
                gm_=[gm(1) gm(2)];
            case{-40}
                val=bc_val.top;
                gm_=[gm(3) gm(4)];
            otherwise
                error('inhomogeneous Neumann');
        end
        b(gm_(:)) = b(gm_(:)) + val*Le/2;

    end

    % Robin
    if(     (Kp==-10 && bc_type(1)==4) || ...
            (Kp==-20 && bc_type(2)==4) || ...
            (Kp==-30 && bc_type(3)==4) || ...
            (Kp==-40 && bc_type(4)==4) )

        switch(Kp)
            case{-10}
                val=bc_val.left;
                gm_=[gm(1) gm(4)];
            case{-20}
                val=bc_val.right;
                gm_=[gm(2) gm(3)];
            case{-30}
                val=bc_val.bottom;
                gm_=[gm(1) gm(2)];
            case{-40}
                val=bc_val.top;
                gm_=[gm(3) gm(4)];
            otherwise
                error('Robin');
        end

        % add 0.5 (phi,v) to the lhs
        aux=zeros(nvm,nvm);
        aux(IDm,IDm) = aux(IDm,IDm) + 0.5 * Le * m1d;
        A(gm(:),gm(:)) = A(gm(:),gm(:)) + aux;
        % add 2(Jinc,v) to the rhs
        b(gm_(:)) = b(gm_(:)) + 2*val*Le/2;

    end

end

% spy(A)
t2=cputime;
fprintf('Assembly time = %g \n',t2-t1);


%------------------------------------------------
%solve

t1=cputime;

z=A\b;

t2=cputime;
fprintf('Solver time   = %g \n\n',t2-t1);

fprintf('Solution: Min = %+E \n          Max = %+E \n\n',min(z),max(z));
%------------------------------------------------

return
end