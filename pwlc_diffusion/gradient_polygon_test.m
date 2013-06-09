% function [M,K,f]=gradient_polygon_test(poly,vert)
% vertices are entered anti-clockwise
clear all; close all; clc;

poly=[1 2 3]
vert=[0 0; 1 0; 0 1]

% centroid
C=mean(vert);
% alpha coef
nv=length(poly);
alpha=1/nv;

g=zeros(2,nv);

% loop over sides
for iside=1:nv
    % pick 1st vertex
    irow1=iside;
    % pick 2ndt vertex
    irow2=irow1+1; if(irow2>nv), irow2=1; end
    % compute J=2.Area for that side
    A=vert(irow1,:); B=vert(irow2,:);
    Jac=[(B-A)' (C-A)']
    iJt=inv(Jac')
    % stiffness matrix
    g_side=[-1 1; (alpha-1) alpha];
    g([1 2],[irow1 irow2]) = g([1 2],[irow1 irow2]) + iJt*g_side;
end
g
