function [edg2poly,edg2vert,new_edge]=...
    is_edge_already_recorded(ed,edg2poly,edg2vert,ipoly,vert_link,new_edge)
% ed: contains the 2 vertex IDs of the edge
% use the common grid vertices to find recorded edges

v1=vert_link(ed(1));
v2=vert_link(ed(2));

ind1 = find( vert_link(edg2vert(:,1))     == v1 );
ind2 = find( vert_link(edg2vert(ind1,2))  == v2 );
ind1_= find( vert_link(edg2vert(:,1))     == v2 );
ind2_= find( vert_link(edg2vert(ind1_,2)) == v1 );

if(length(ind2)>1), error('duplicated edges!'); end
if(length(ind2_)>1), error('duplicated edges!'); end
if(length(ind2)==1 & length(ind2_)==1), error('duplicated edges!'); end

if(length(ind2)==0 & length(ind2_)==0)
    new_edge=new_edge+1;
    edg2vert(new_edge,1:2)= ed;
    edg2poly(new_edge,1)=ipoly; % 1
    edg2poly(new_edge,2)=-1; % bc
else
    if(length(ind2)==1)
        work=ind1(ind2);
    else
        work=ind1_(ind2_);
    end
    ied=work(1);
    ied
    edg2vert(ied,3:4)= ed;
    edg2poly(ied,2)=ipoly
end

return
end