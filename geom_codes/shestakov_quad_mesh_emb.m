clear all; close all; clc;

% This routine generates an nm by nm mesh, commonly known as a
% Shestakov mesh, with 0 < r < rm and 0 < z < zm. The randomness
% is controlled by parameter "a", where 0 <= a <= 0.5. A value of 
% a = 0.5 gives a rectangular % mesh. Boomerang zones are allowed 
% but bowties are not allowed.
%
% Matlab version by Jean C. Ragusa, Texas A&M University, April 1, 2013
% (jean.ragusa@tamu.edu)
%
% based on the original Fortran program by Alek Sheshakov (LLNL); that 
% program is found in: Nuclear Science & Engineering, Vol 105, pp.88-104 
% (1990), "Test Problems in Radiative Transfer Calculations", by A. Shestakov, 
% D. Kershaw, and G. Zimmerman
%
% Many thanks to Michael Hall (LANL) for providing a copy of the F77 
% program


logi_save_plots = false;
logi_write_file = true;


% number of subdivisions of the original rectangle
nc = 7;
% random parameter
a  = 0.20;
% rectangle dimensions
L = 1;
rm = L;
zm = L;

% allocate memory
nm = 2^nc + 1;
ind = nm-1;
r=zeros(nm,nm);
z=zeros(nm,nm);

% initialize 4 corners
for i = 0:1
    k = 1 + i*ind;
    for j = 0:1
        l = 1 + j*ind;
        r(k,l) = L*(1-i);
        z(k,l) = L*j;
    end
end

% fill in the rest of the points
for nl = 0:nc-1
    nn = 2^nl;
    inc = ind/nn;
    inch = inc/2;
    for k = 1:nn
        k1 = 1+(k-1)*inc;
        k2 = k1+inc;
        k3 = k1+inch;
        for l = 1:nn
            l1 = 1+(l-1)*inc;
            l2 = l1+inc;
            l3 = l1+inch;
            if (l==1)
                ar = a+rand(1,1)*(1-2*a);
                r(k3,l1) = ar*r(k1,l1)+(1-ar)*r(k2,l1);
                z(k3,l1) = ar*z(k1,l1)+(1-ar)*z(k2,l1);
            end
            ar = a+rand(1,1)*(1-2*a);
            r(k2,l3) = ar*r(k2,l1)+(1-ar)*r(k2,l2);
            z(k2,l3) = ar*z(k2,l1)+(1-ar)*z(k2,l2);
            ar = a+rand(1,1)*(1-2*a);
            r(k3,l2) = ar*r(k1,l2)+(1-ar)*r(k2,l2);
            z(k3,l2) = ar*z(k1,l2)+(1-ar)*z(k2,l2);
            if (k==1)
                ar = a+rand(1,1)*(1-2*a);
                r(k1,l3) = ar*r(k1,l1)+(1-ar)*r(k1,l2);
                z(k1,l3) = ar*z(k1,l1)+(1-ar)*z(k1,l2);
            end
            ar = a+rand(1,1)*(1-2*a);
            br = a+rand(1,1)*(1-2*a);
            r1 = r(k1,l1);
            r2 = r(k2,l1);
            r3 = r(k2,l2);
            r4 = r(k1,l2);
            z1 = z(k1,l1);
            z2 = z(k2,l1);
            z3 = z(k2,l2);
            z4 = z(k1,l2);
            % check for boomerang zones
            det2 = (r2-r1)*(z3-z2)-(r3-r2)*(z2-z1);
            det3 = (r3-r2)*(z4-z3)-(r4-r3)*(z3-z2);
            det4 = (r4-r3)*(z1-z4)-(r1-r4)*(z4-z3);
            det1 = (r1-r4)*(z2-z1)-(r2-r1)*(z1-z4);
            if (det2>0)
                d = (r4-r3)*(z2-z1)-(r2-r1)*(z4-z3);
                r3p = ((r2-r1)*(r4*z3-r3*z4)-(r4-r3)*(r2*z1-r1*z2))/d;
                z3p = ((z2-z1)*(r4*z3-r3*z4)-(z4-z3)*(r2*z1-r1*z2))/d;
                d = (r4-r1)*(z2-z3)-(r2-r3)*(z4-z1);
                r1p = ((r2-r3)*(r4*z1-r1*z4)-(r4-r1)*(r2*z3-r3*z2))/d;
                z1p = ((z2-z3)*(r4*z1-r1*z4)-(z4-z1)*(r2*z3-r3*z2))/d;
                r3 = r3p;
                z3 = z3p;
                r1 = r1p;
                z1 = z1p;
            elseif (det3>0)
                d = (r1-r4)*(z3-z2)-(r3-r2)*(z1-z4);
                r4p = ((r3-r2)*(r1*z4-r4*z1)-(r1-r4)*(r3*z2-r2*z3))/d;
                z4p = ((z3-z2)*(r1*z4-r4*z1)-(z1-z4)*(r3*z2-r2*z3))/d;
                d = (r1-r2)*(z3-z4)-(r3-r4)*(z1-z2);
                r2p = ((r3-r4)*(r1*z2-r2*z1)-(r1-r2)*(r3*z4-r4*z3))/d;
                z2p = ((z3-z4)*(r1*z2-r2*z1)-(z1-z2)*(r3*z4-r4*z3))/d;
                r4 = r4p;
                z4 = z4p;
                r2 = r2p;
                z2 = z2p;
            elseif (det4>0)
                d = (r2-r1)*(z4-z3)-(r4-r3)*(z2-z1);
                r1p = ((r4-r3)*(r2*z1-r1*z2)-(r2-r1)*(r4*z3-r3*z4))/d;
                z1p = ((z4-z3)*(r2*z1-r1*z2)-(z2-z1)*(r4*z3-r3*z4))/d;
                d = (r2-r3)*(z4-z1)-(r4-r1)*(z2-z3);
                r3p = ((r4-r1)*(r2*z3-r3*z2)-(r2-r3)*(r4*z1-r1*z4))/d;
                z3p = ((z4-z1)*(r2*z3-r3*z2)-(z2-z3)*(r4*z1-r1*z4))/d;
                r1 = r1p;
                z1 = z1p;
                r3 = r3p;
                z3 = z3p;
            elseif (det1>0)
                d = (r3-r2)*(z1-z4)-(r1-r4)*(z3-z2);
                r2p = ((r1-r4)*(r3*z2-r2*z3)-(r3-r2)*(r1*z4-r4*z1))/d;
                z2p = ((z1-z4)*(r3*z2-r2*z3)-(z3-z2)*(r1*z4-r4*z1))/d;
                d = (r3-r4)*(z1-z2)-(r1-r2)*(z3-z4);
                r4p = ((r1-r2)*(r3*z4-r4*z3)-(r3-r4)*(r1*z2-r2*z1))/d;
                z4p = ((z1-z2)*(r3*z4-r4*z3)-(z3-z4)*(r1*z2-r2*z1))/d;
                r2 = r2p;
                z2 = z2p;
                r4 = r4p;
                z4 = z4p;
            end
            r(k3,l3) = ar*br*r1 + ar*(1-br)*r4 + (1-ar)*br*r2 + (1-ar)*(1-br)*r3;
            z(k3,l3) = ar*br*z1 + ar*(1-br)*z4 + (1-ar)*br*z2 + (1-ar)*(1-br)*z3;
        end
    end
end

% plotting in matlab
figure(1)
surf(r,z,ones(nm,nm),'FaceColor','white')
view(0,90)
% nc = 4;
% nm = 2^nc + 1;
figure(2);
q1=(ceil(nc/2));
q2=ceil(nc/q1);
for k=1:nc
    skip=2^(k-1);
    subplot(q1,q2,k);
    nm_ = 2^(nc-k+1)+1;
    surf(r(1:skip:end,1:skip:end),z(1:skip:end,1:skip:end),ones(nm_,nm_),'FaceColor','white');
     set(gca,'PlotBoxAspectRatio',[5 5 1])
     view(0,90);
end

output_file1=strcat('.\figs\shestakov_quad_L',int2str(L),'_nc',int2str(nc),'_emb','_a',num2str(a));
if(logi_save_plots)
    print('-dpdf',strcat(output_file1,'.pdf'));
    print('-dpng',strcat(output_file1,'.png'));
    saveas(gcf,strcat(output_file1,'.fig'),'fig');
end

%---------------------------------------------
%
if(~logi_write_file), return; end;
%
%---------------------------------------------
% save txt file
matID=1;
srcID=1;

% date and time;
[yr, mo, da, hr, mi, s] = datevec(now);

% save embedded meshes
r_sav=r;
z_sav=z;
for k=1:nc
    
    nm_ = 2^(nc-k+1)+1;
    skip=2^(k-1);
    r=r_sav(1:skip:end,1:skip:end);
    z=z_sav(1:skip:end,1:skip:end);

    %%%%%%%%%%%%%%%%%%%%%
    output_file_k=strcat('.\figs\shestakov_quad_L',int2str(L),...
        '_nc',int2str(nc),'_emb',int2str(nc-k+1),'_a',num2str(a),'.txt');
    fid=fopen(output_file_k,'w');
    
    fprintf(fid,'# Date: %d/%d/%d   Time: %d:%d\n', mo, da, yr, hr, mi);
    
    fprintf(fid,'# dimensions \n');
    fprintf(fid,'%g %g \n',L,L);
    
    n = (nm_-1);
    ncells = n*n;
    
    fprintf(fid,'# connectivity \n');
    fprintf(fid,'%d\n',ncells);
    for iel=1:ncells
        skip = 4*(iel-1);
        i1 = skip + 1;
        i2 = skip + 2;
        i3 = skip + 3;
        i4 = skip + 4;
        fprintf(fid,'%d %d %d %d %d %d %d \n',4,i1,i2,i3,i4,matID,srcID);
    end
    
    fprintf(fid,'# DG vertices (counter-clockwise) \n');
    fprintf(fid,'%d\n',4*ncells);
    for i=1:n
        for j=1:n
            fprintf(fid,' %g %g \n %g %g \n %g %g \n %g %g \n',...
                r(i,j)    ,z(i,j)    ,...
                r(i,j+1)  ,z(i,j+1)  ,...
                r(i+1,j+1),z(i+1,j+1),...
                r(i+1,j)  ,z(i+1,j)  );
        end
    end
    
    fprintf(fid,'# grid vertices (counter-clockwise) \n');
    fprintf(fid,'%d\n',(n+1)^2);
    for i=1:n+1
        for j=1:n+1
            fprintf(fid,'%g %g \n',r(i,j),z(i,j) );
        end
    end
    
    fclose(fid);
    %%%%
end % embedded meshes save


