function o = intersect_jcr(cell,inside,arg3)

if(nargin==3)
    to_sort=true;
else
    to_sort=false;
end


nc=length(cell);
ni=length(inside);
o=[];

for i=1:ni
    val=inside(i);
    ind=find(cell==val);
    if(~isempty(ind))
        o=[o ind];
    end
end

if(to_sort)
    ind=find(abs(diff(o))~=1);
    if(length(ind)>1), error('intersect_jcr');end
    if(length(ind)==1)
        oo=[ o(ind+1:end) o(1:ind)];
        o=oo;
    end
end

return
end
