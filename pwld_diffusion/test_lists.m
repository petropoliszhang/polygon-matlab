clear all; close all; clc;

cells=1:20;

next_ref_lev=zeros(20,1);

cell_refine = [1 11 6 3 19 18 14 15 16];

next_ref_lev(1)  = 1;
next_ref_lev(6)  = 2;
next_ref_lev(3)  = 3;
next_ref_lev(19) = 10;
next_ref_lev(18) = 5;
next_ref_lev(11) = 2;
next_ref_lev(14) = 4;
next_ref_lev(15) = 4;
next_ref_lev(16) = 4;

ref_to_process = next_ref_lev(cell_refine);
level_values = unique(ref_to_process);
how_many_levels=length(level_values);

ordered_cell_refine=[];

for k=1:how_many_levels
    ind = find( ref_to_process == level_values(k) );
    ref_batch{k}  =ind;
    ordered_cell_refine = [ ordered_cell_refine cell_refine(ind)];
%     ref_to_process(ind)
end
ordered_cell_refine


[sorted_next_ref_lev,sorting_order]=sort(next_ref_lev)
for k=1:how_many_levels
    ind = find( next_ref_lev == level_values(k) );
    ncells_per_lev(k) = length(ind);
end
[level_values ncells_per_lev']



[sorted_next_ref_lev,sorting_order]=sort(ref_to_process)

