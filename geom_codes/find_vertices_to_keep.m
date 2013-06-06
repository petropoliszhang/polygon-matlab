function vertices_to_keep = find_vertices_to_keep(cell,inside)

nc=length(cell);

% downselect to the inside points plus the first 2 outside points
first_inside =( find(cell==inside(1)  ) );
last_inside  =( find(cell==inside(end)) );

first_before = first_inside-1;
if (first_before==0), first_before=nc; end

check = find( inside == cell(first_before) );
if(~isempty(check))
    first_before = first_inside+1;
    if (first_before==nc+1), first_before=1; end
end

last_after = last_inside+1;
if (last_after==nc+1), last_after=1; end

check = find( inside == cell(last_after) );
if(~isempty(check))
    last_after = last_inside-1;
    if (last_after==0), last_after=nc; end
end

vertices_to_keep = [cell(first_before) inside cell(last_after)];

return

end