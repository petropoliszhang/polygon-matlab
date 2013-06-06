function [x,y]=find_intersect_with_square(interior_pt,exterior_pt,width)

I=interior_pt;
A=exterior_pt;
L=width;

aux=I;
I(1)=aux(2);
I(2)=aux(1);
aux=A;
A(1)=aux(2);
A(2)=aux(1);

IA=A-I;

% x-left
ind=1;
r(ind)=0;
t(ind)=IA(1)/(r(ind)-I(1));
s(ind)=IA(2)/t(ind)+I(2);
t(ind)=1/t(ind);
if(t(ind)<0),t(ind)=Inf;end

% x-right
ind=2;
r(ind)=L;
t(ind)=IA(1)/(r(ind)-I(1));
s(ind)=IA(2)/t(ind)+I(2);
t(ind)=1/t(ind);
if(t(ind)<0),t(ind)=Inf;end


% y-left
ind=3;
s(ind)=0;
t(ind)=IA(2)/(s(ind)-I(2));
r(ind)=IA(1)/t(ind)+I(1);
t(ind)=1/t(ind);
if(t(ind)<0),t(ind)=Inf;end

% y-right
ind=4;
s(ind)=L;
t(ind)=IA(2)/(s(ind)-I(2));
r(ind)=IA(1)/t(ind)+I(1);
t(ind)=1/t(ind);
if(t(ind)<0),t(ind)=Inf;end

%
[tt,ind]=min(t);
if(tt<0 |tt>1), error('tt<0 |tt>1');end

y = r(ind);
x = s(ind);

return
end
