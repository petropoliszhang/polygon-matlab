clear all; close all; clc;

cells=1:20;

next_ref_lev=zeros(20,1);

cell_refine = [1 11 6 3 19 18 ];

next_ref_lev(1)  = 1;
next_ref_lev(6)  = 2;
next_ref_lev(3)  = 3;
next_ref_lev(19) = 0;
next_ref_lev(18) = 5;
next_ref_lev(11) = 2;

ref_to_process = next_ref_lev(cell_refine);
level_values = unique(ref_to_process);
how_many_levels=length(level_values);
for k=1:how_many_levels
    ind= find( ref_to_process == level_values(k) );
    ref_batch{k}  =ind;
    cell_refine(ind)
%     ref_to_process(ind)
end
