function convexity_util( connectivity, nel, vert)

nb_concave=0;
for iel=1:nel
    g=connectivity{iel};
    nsides=length(g);
    gg  = [g(end); g ; g(1)];
    v = vert(gg,:);
%     oo=ones(nsides,2);
%     oo(:,1)=0;
%     test = vert(g,:) - oo;
%     ind=find(sum(test.^2,2)<eps, 1);
%     if(~isempty(ind))
%         iel
%     end
    found_concave=false;
    for k=2:nsides+1
        %                +B
        %               /
        %     ccw      /
        %             /
        %  C+--------+ A
        %
        A=[v(k  ,:) 0];
        B=[v(k+1,:) 0];
        C=[v(k-1,:) 0];
        aux=cross(B-A,C-A);
        angle = mod(atan2(aux(3),dot(B-A,C-A)),2*pi);
        %         angle = mod(atan2(x1*y2-x2*y1,x1*x2+y1*y2),2*pi);
        if(angle>pi)
            found_concave=true;
        end
    end
    if found_concave
        nb_concave=nb_concave+1;
    end
end
fprintf('Nbr of concave polygons %d out of %d cells \n',nb_concave,nel);
return
end
