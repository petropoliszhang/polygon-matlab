clear all; close all; clc;

A=[1 pi];
B=[5 5.33];
C=[exp(1) 10+1/7];
vert=[A; B; C]
% vert=[0 0; 1 0; 0 1]

f=@(x,y) 1+0*(x.*y);
N=8;
[X,Y,Wx,Wy]=triquad(N,vert);
Q=Wx'*feval(f,X,Y)*Wy
test_area=12.1223786902332

%%%%%
fprintf('\n \n now the other test \n');
g=@(x,y) cos(x/10).*sin(y/10);
test=6.65088376030847;
for i=1:8,
    [X,Y,Wx,Wy]=triquad(i,vert);
    Q=Wx'*feval(g,X,Y)*Wy;
    [Q Q-test]
end

%%%%%
fprintf('\n \n now the other test \n');
for i=1:8,
    [X,Y,Wx,Wy]=triquad(i,vert);
    Q=Wx'*(feval(g,X,Y)-feval(f,X,Y))*Wy;
    [Q Q-(test-test_area)]
end

%%%%%
fprintf('\n \n now the other test \n');
test=2.61709518928211;
for i=1:8,
    [X,Y,Wx,Wy]=triquad(i,vert);
    Q=Wx'*(feval(g,X,Y)-feval(f,X,Y)).^2*Wy;
    [Q Q-test]
end
