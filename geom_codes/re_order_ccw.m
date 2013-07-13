function [ new_ci ] = re_order_ccw( vert, old_connectivity )
% re-order poly ccw

[x2, y2] = poly2ccw(vert(:,1), vert(:,2));
% re-order connectivity so that the poly is actually ccw
new_ci=[];
for k=1:length(x2)
    vk=[x2(k) y2(k)];
    dif = sum(sqrt( (vert-kron(ones(length(x2),1),vk)).^2 ) , 2);
    ind=find(dif<1e-10);
    if(length(ind)~=1)
        ind
        error('the vk vertex should have showed up once only....');
    end
    new_ci =[ new_ci old_connectivity(ind)];
end

return
end

