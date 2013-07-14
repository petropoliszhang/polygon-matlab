function [ new_ci2 ] = re_order_ccw( vert, old_connectivity )

% first, remove vertices that are the same
new_ci=[];
to_omit=[];
for k=1:length(old_connectivity)
    vk=vert(k,:);
    dif = sum(sqrt( (vert-kron(ones(length(old_connectivity),1),vk)).^2 ) , 2);
    % index in vert where vk is found
    ind=find(dif<1e-10);
    indk=find(ind==k);
    if(length(ind)~=1)
        warning('the vk vertex should have showed up once only....');
        ind
        ind2=find(ind>k);
        to_omit=[to_omit old_connectivity(ind(ind2)) ];
    end
    new_ci =[ new_ci old_connectivity(ind(indk))];
end

to_omit=unique(to_omit);
for k=1:length(to_omit)
    ind=find(new_ci==to_omit(k));
    ind=sort(ind,'descend');
    for n=1:length(ind)
        new_ci(ind(n))=[];
        vert(ind(n),:)=[];
    end
end

% then re-order poly ccw
[x2, y2] = poly2ccw(vert(:,1), vert(:,2));
% re-order connectivity so that the poly is actually ccw
new_ci2=[];
for k=1:length(x2)
    vk=[x2(k) y2(k)];
    dif = sum(sqrt( (vert-kron(ones(length(x2),1),vk)).^2 ) , 2);
    % index in vert where vk is found
    ind=find(dif<1e-10);
    if(length(ind)~=1)
        error('the vk vertex should have showed up once only....');
        ind
    end
    new_ci2 =[ new_ci2 new_ci(ind)];
end

return
end

