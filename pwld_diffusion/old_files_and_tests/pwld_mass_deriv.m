clear all; clc;

% vertices are entered anti-clockwise
poly=[1 2 3]
vert=[0 0; 1 0; 0 1]

% centroid
C=mean(vert);
% alpha coef
nv=length(poly);
alpha=1/nv;

% mass and stiffness matrices for the entire polygon
M=zeros(nv,nv);
K=zeros(nv,nv);
f=zeros(nv,1);
% mass integrals over one side:
M_side=zeros(nv,nv);
f_side=zeros(nv,1);
% integrals of b^i_1/2 b^i_1/2 :
m=[2 1;1 2]/24+alpha*(1+alpha)*[1 1;1 1]/12;
% integrals of b^i_1/2 b^i_c :
aux1=alpha/24*(1+2*alpha); % a/24 + a^2/12
% integrals of b^i_c b^i_c :
aux2=alpha^2/12;
% mass matrix
M_side([1 2],[1 2]) = m;
M_side([1 2],3:end) = aux1;
M_side(3:end,[1 2]) = aux1;
M_side(3:end,3:end) = aux2;
% list of vertices, when looping over sides
list_vert=1:nv;
f_side=[1+alpha; 1+alpha; alpha]/6;

% loop over sides
for iside=1:nv
    % pick 1st vertex
    irow1=iside;
    % pick 2ndt vertex
    irow2=irow1+1; if(irow2>nv), irow2=1; end
    % compute J=2.Area for that side
    A=vert(irow1,:); B=vert(irow2,:);
    Jac_i = cross([B-A 0]',[C-A 0]'); % AB ^ AC
    if(Jac_i(3)<0), error('not clockwise ordering in polygon'); end
    det_J_i=norm(Jac_i,2);
    % add contribution from M_side
    M(list_vert,list_vert) = M(list_vert,list_vert) + det_J_i*M_side;
    % add contribution from f_side
    f(list_vert) = f(list_vert) + det_J_i*f_side;
    % stiffness matrix
    a = alpha;
    L1=norm(B-A); L2=norm(C-B); L3=norm(A-C);
    r1 = L1^2/(2*det_J_i);
    r2 = L2^2/(2*det_J_i);
    r3 = L3^2/(2*det_J_i);
    [r1 r2 r3]
    % sanity check:
    test = 2*r1*(r2+r3)+2*r2*r3-(r1^2+r2^2+r3^2)-1;
    if(abs(test)>1e-13), test, error('r_i error in stiffness matrix coefficients'); end
    % stiffness matrix for side i
    K_side=[ ...
        (-1 + a)*a*r1 - (-1 + a)*r2 + a*r3,  ((1 - 2*a + 2*a^2)*r1 - r2 - r3)/2,  (a*((-1 + 2*a)*r1 - r2 + r3))/2;
        ((1 - 2*a + 2*a^2)*r1 - r2 - r3)/2,  (-1 + a)*a*r1 + a*r2 - (-1 + a)*r3,  (a*((-1 + 2*a)*r1 + r2 - r3))/2;
        (a*((-1 + 2*a)*r1 - r2 + r3))/2   ,  (a*((-1 + 2*a)*r1 + r2 - r3))/2,     a^2*r1
        ];
    % K_side=[...
    %         (-1 + a)*((-2 + a)*r1 + r2) - (-2 + a)*r3, ((-1 - 2*a + 2*a^2)*r1 + r2 - 3*r3)/2, (a*((-3 + 2*a)*r1 + r2 - r3))/2 ;
    %         ((-1 - 2*a + 2*a^2)*r1 + r2 - 3*r3)/2    , a*(1 + a)*r1 - a*r2 +  (1 + a)*r3    , (a*((1 + 2*a)*r1 - r2 + r3))/2 ;
    %         (a*((-3 + 2*a)*r1 + r2 - r3))/2          , (a*((1 + 2*a)*r1 - r2 + r3))/2       ,  a^2*r1 ];
    %     K_side=0.5*[ ...
    %         2*((1 - a)*r1 +  a*(r2 + (-1 + a)*r3)), -r1 - r2 + (1 - 2*a + 2*a^2)*r3 , a*(-r1 + r2 + (-1 + 2*a)*r3);
    %         -r1 - r2 + (1 - 2*a + 2*a^2)*r3       , 2*(a*r1 + (-1 + a)*(-r2 + a*r3)), a*(r1 - r2 - r3 + 2*a*r3);
    %         a*(-r1 + r2 + (-1 + 2*a)*r3)          , a*(r1 - r2 - r3 + 2*a*r3)       , 2*a^2*r3 ];
    K_side
    % add contribution from K_side
    K(list_vert,list_vert) = K(list_vert,list_vert) + K_side;
    % shift vertex IDs
    list_vert=[list_vert list_vert(1)];
    list_vert(1)=[];
end

% M*24
f
% K
% alpha=0;
% m=[2 1;1 2]/24+alpha*(1+alpha)*[1 1;1 1]/12;
% aux=alpha/24*(1+2*alpha);
% m
% for iside=1:nv
%     irow1=iside;
%     irow2=irow1+1; if(irow2>nv), irow2=1; end
%     A=vert(irow1,:);
%     B=vert(irow2,:);
%     J_i = cross([B-A 0]',[C-A 0]'); % AB ^ AC
%     J_i=norm(J_i,2);
%     M([i1 i2],[i1 i2]) = M([i1 i2],[i1 i2]) + J_i*m;
% %     M(i1,[i1 i2]) = M(i1,[i1 i2]) + J_i*m(1,:);
% %     M(i2,[i1 i2]) = M(i2,[i1 i2]) + J_i*m(2,:);
%     for j=1:nv
%         if(j==i1|j==i2),continue;end
%         j
%         M([i1 i2],j) = M([i1 i2],j) + J_i*aux*[1;1];
%     end
% end
% M