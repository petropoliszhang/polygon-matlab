% testing this
clc;clear all;

del=[];
i=123;
j=36;

gi=      [5 6 100 7 ];
gj=[20 10 5 6 7 1 22 33 44 55 ];

% combine the vertex entries
combined =[gi gj];
uniq = unique(combined);
nc=length(combined);
nu=length(uniq);
if ( nu> nc )
    error('nu cannot be > than nc');
end
if ( nu < nc-2)
    % need to investigate
    %reverse order for gj
    n_gj=length(gj);
    for k=1:n_gj
        gjr(k)=gj(n_gj+1-k);
    end
    % find the intersection, but the result unfortunetaly comes out sorted
    [inter,ia,ib]=intersect(gi,gjr)
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
    fprintf('sign = %d \n',sign);
    % now, find the position of the first common element in gi
    [dum,posi]=min(ia);
    fprintf('posi = %d \n\n',posi);
    flag=0;
    gi
    gjr
    indj=ib(posi);
    for k=1:length(ia)
        indi=ia(posi);
        vali=gi(indi);
        fprintf('indi = %d vali = %g\n',indi,vali);
        valj=gjr(indj);
        fprintf('indj = %d valj = %g\n\n',indj,valj);
        if(abs(vali-valj)>eps)
            flag=1;
            warning(' aaa ' );
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
del