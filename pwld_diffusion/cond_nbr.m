function cond_nbr
close all; clc;

basedir = '.\results\convergence\';

% -------------
% -------------
% --- QUADS ---
% -------------
% -------------
dir =strcat(basedir,'quads\');

% uniform
str = strcat(dir,'uniform\mms-1\reg_quad_mms_1.mat');
iter_anal(str);

% random
% geo = 'rand_quad_mms_1_L1_nc8_emb1_a0.66.mat';
geo = 'rand_quad_mms_1_L1_nc_emb1_a0.95.mat';
for k=7:7
    k1=strfind(geo,'_emb');
    str = strcat(dir,'rand\',geo(1:k1+3),int2str(k),geo(k1+5:end));
end
iter_anal(str);

% smooth
geo = 'smoo_quad_mms_1_L1_nc7_emb1_a0.15.mat';
for k=7:7
    k1=strfind(geo,'_emb');
    str = strcat(dir,'smooth\',geo(1:k1+3),int2str(k),geo(k1+5:end));
end
iter_anal(str);

% shestakov
geo = 'shes_quad_mms_1_L1_nc8_emb1_a0.1.mat';
for k=7:7
    k1=strfind(geo,'_emb');
    str = strcat(dir,'shes\',geo(1:k1+3),int2str(k),geo(k1+5:end));
end
iter_anal(str);


geo = 'shes_quad_mms_1_L1_nc7_emb1_a0.15.mat';
for k=7:7
    k1=strfind(geo,'_emb');
    str = strcat(dir,'shes\',geo(1:k1+3),int2str(k),geo(k1+5:end));
end
iter_anal(str);


geo = 'shes_quad_mms_1_L1_nc7_emb1_a0.2.mat';
for k=7:7
    k1=strfind(geo,'_emb');
    str = strcat(dir,'shes\',geo(1:k1+3),int2str(k),geo(k1+5:end));
end
iter_anal(str);


geo = 'shes_quad_mms_1_L1_nc7_emb1_a0.25.mat';
for k=7:7
    k1=strfind(geo,'_emb');
    str = strcat(dir,'shes\',geo(1:k1+3),int2str(k),geo(k1+5:end));
end
iter_anal(str);


geo = 'shes_quad_mms_1_L1_nc7_emb1_a0.3.mat';
for k=7:7
    k1=strfind(geo,'_emb');
    str = strcat(dir,'shes\',geo(1:k1+3),int2str(k),geo(k1+5:end));
end
iter_anal(str);


geo = 'shes_quad_mms_1_L1_nc7_emb1_a0.35.mat';
for k=7:7
    k1=strfind(geo,'_emb');
    str = strcat(dir,'shes\',geo(1:k1+3),int2str(k),geo(k1+5:end));
end
iter_anal(str);


geo = 'shes_quad_mms_1_L1_nc7_emb1_a0.45.mat';
for k=7:7
    k1=strfind(geo,'_emb');
    str = strcat(dir,'shes\',geo(1:k1+3),int2str(k),geo(k1+5:end));
end
iter_anal(str);


% z-mesh 0.05
w=0.05;
str = strcat(dir,'z-mesh\zzzz_quad_mms_1_a0.',num2str(w*100,'%2.2d'),'.mat');
iter_anal(str);

w=0.1;
str = strcat(dir,'z-mesh\zzzz_quad_mms_1_L1_n20_a0.',num2str(w*10),'.mat');
iter_anal(str);

w=0.15;
str = strcat(dir,'z-mesh\zzzz_quad_mms_1_L1_n20_a0.',num2str(w*100,'%2.2d'),'.mat');
iter_anal(str);

w=0.2;
str = strcat(dir,'z-mesh\zzzz_quad_mms_1_L1_n20_a0.',num2str(w*10),'.mat');
iter_anal(str);

w=0.25;
str = strcat(dir,'z-mesh\zzzz_quad_mms_1_L1_n20_a0.',num2str(w*100),'.mat');
iter_anal(str);

w=0.30;
str = strcat(dir,'z-mesh\zzzz_quad_mms_1_L1_n20_a0.',num2str(w*10),'.mat');
iter_anal(str);

w=0.35;
str = strcat(dir,'z-mesh\zzzz_quad_mms_1_L1_n20_a0.',num2str(w*100),'.mat');
iter_anal(str);

w=0.40;
str = strcat(dir,'z-mesh\zzzz_quad_mms_1_L1_n20_a0.',num2str(w*10),'.mat');
iter_anal(str);



return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function iter_anal(str)

expression=sprintf('%s %s','load',str);eval(expression);

fprintf('\n\n\n===============================\n')
fprintf('===============================\n');
fprintf('%s \n',str);
fprintf('===============================\n');
fprintf('===============================\n');

[nel ndof]

[z,A,b] = DG_assemble_solve( ndof,nel,n_edge,vert,connectivity,edg2poly,edg2vert,edg_normal,edg_perp,C_pen,C_pen_bd,...
    i_mat,c_diff,sigma_a,i_src,S_ext,logi_mms,mms,n_quad,bc_type,bc_val );
% A=speye(ndof);b=ones(ndof,1);

fprintf('condest A = %g \n',condest(A));

mg=false;

if(mg)
    
    [x,flag,relres,iter,resvec]=agmg(A,b,1,1e-11,1000,1);
    fprintf('flag %d \n',flag);
    fprintf('nbr iter %d \n',iter);
    fprintf('relres %g \n',relres);
    fprintf('norm x-z = %g\n\n\n',norm(x-z));
    
else
    
    niter=1000;
    
    % A=D+E+F
    D=diag(diag(A));
    E=tril(A,-1);%spy(E);
    F=triu(A,+1);%spy(F);
    
    % fid=1;
    % eigplot(full(A),fid);
    % title(sprintf('No precond, cond=%g',condest(A)));
    % eval(sprintf('print -f%i -dpng %i.png',fid,fid))
    [x,flag,relres,iter,resvec]=pcg(A,b,1e-10,niter);
    flag
    iter
    relres
    fprintf('nbr iter unprec %d \n',iter);
    fprintf('norm x-z = %g\n',norm(x-z));
    % figure(fid);
    % semilogy(0:iter,resvec/norm(b),'-o');
    % xlabel('iteration number');ylabel('relative residual');
    
    
    % fid=fid+1; M=D\A;
    % condest(M)
    % eigplot(full(M),fid);
    % title(sprintf('Jacobi Precond, cond=%g',condest(M)));
    % eval(sprintf('print -f%i -dpng %i.png',fid,fid))
    [x,flag,relres,iter,resvec]=pcg(A,b,1e-10,niter,D);
    flag
    iter
    relres
    fprintf('nbr iter jac %d \n',iter);
    fprintf('norm x-z = %g\n',norm(x-z));
    % figure(fid);
    % semilogy(0:iter,resvec/norm(b),'-o');
    % xlabel('iteration number');ylabel('relative residual');
    
    w=1;
    GS=D+w*E;
    % fid=fid+1; M=GS\A;
    % condest(M)
    % eigplot(full(M),fid);
    % title(sprintf('Gauss-Seidel Precond, cond=%g',condest(M)));
    % eval(sprintf('print -f%i -dpng %i.png',fid,fid))
    [x,flag,relres,iter,resvec]=pcg(A,b,1e-10,niter,GS);
    flag
    iter
    relres
    fprintf('nbr iter GS %d \n',iter);
    fprintf('norm x-z = %g\n',norm(x-z));
    % figure(fid);
    % semilogy(0:iter,resvec/norm(b),'-o');
    % xlabel('iteration number');ylabel('relative residual');
    
    w=1.3;
    SOR=D+w*E;
    % fid=fid+1; M=SOR\A;
    % condest(M)
    % eigplot(full(M),fid);
    % title(sprintf('SOR w=1.3 Precond, cond=%g',condest(M)));
    % eval(sprintf('print -f%i -dpng %i.png',fid,fid))
    [x,flag,relres,iter,resvec]=pcg(A,b,1e-10,niter,SOR);
    flag
    iter
    relres
    fprintf('nbr iter SOR1.3 %d \n',iter);
    fprintf('norm x-z = %g\n',norm(x-z));
    % figure(fid);
    % semilogy(0:iter,resvec/norm(b),'-o');
    % xlabel('iteration number');ylabel('relative residual');
    
    w=1.;
    SSOR=(D+w*E)*inv(D)*(D+w*F);
    % fid=fid+1; M=SSOR\A;
    % condest(M)
    % eigplot(full(M),fid);
    % title(sprintf('SSOR w=1. Precond, cond=%g',condest(M)));
    % eval(sprintf('print -f%i -dpng %i.png',fid,fid))
    [x,flag,relres,iter,resvec]=pcg(A,b,1e-10,niter,SSOR);
    flag
    iter
    relres
    fprintf('nbr iter SSOR1 %d \n',iter);
    fprintf('norm x-z = %g\n',norm(x-z));
    % figure(fid);
    % semilogy(0:iter,resvec/norm(b),'-o');
    % xlabel('iteration number');ylabel('relative residual');
    
    w=1.3;
    SSOR=(D+w*E)*inv(D)*(D+w*F);
    % fid=fid+1; M=SSOR\A;
    % condest(M)
    % eigplot(full(M),fid);
    % title(sprintf('SSOR w=1.3 Precond, cond=%g',condest(M)));
    % eval(sprintf('print -f%i -dpng %i.png',fid,fid))
    [x,flag,relres,iter,resvec]=pcg(A,b,1e-10,niter,SSOR);
    flag
    iter
    relres
    fprintf('nbr iter SSOR1.3 %d \n',iter);
    fprintf('norm x-z = %g\n',norm(x-z));
    % figure(fid);
    % semilogy(0:iter,resvec/norm(b),'-o');
    % xlabel('iteration number');ylabel('relative residual');
    
end

return
end

