clear all; close all; clc

% circle
ang=0:0.01:2*pi; 
% ang1=(0.6:0.01:pi+.3) ; 

% pentagon
a1=((0:4) +0.5)*2*pi/5 -pi/2;
P=zeros(5,2);
for i=1:5
    P(i,:) = [cos(a1(i)) sin(a1(i))];
end

figure()
Pent=[P; P(1,:)];
hold on
for i=1:6
    plot(Pent(:,1),Pent(:,2),'k-')
end
% daspect([1,1,1])
leng=norm(diff(P(1:2,:)),2);
rin=1/10*sqrt(25+10*sqrt(5))*leng;
xp=rin*cos(ang); yp=rin*sin(ang);
plot(xp,yp);
rcircum=1/10*sqrt(50+10*sqrt(5))*leng;
xp=rcircum*cos(ang); yp=rcircum*sin(ang);
plot(xp,yp);

A=Pent(1,:)
B=(Pent(3,:)+Pent(4,:))/2
plot([A(1,1),B(1,1)],[A(1,2),B(1,2)],'--k')

text(0,0.1,'h_\perp')
axis off
daspect([1,1,1])

% return

% hexagon
a1=((0:5) +0.5)*2*pi/6 -pi/2;
H=zeros(6,2);
for i=1:6
    H(i,:) = [cos(a1(i)) sin(a1(i))];
end

figure()
Hex=[H; H(1,:)];
hold on
for i=1:7
    plot(Hex(:,1),Hex(:,2),'k-')
end
daspect([1,1,1])

rin=abs(Hex(1,2));
xp=rin*cos(ang); yp=rin*sin(ang);
plot(xp,yp);

A=(Hex(1,:)+Hex(2,:))/2
B=(Hex(4,:)+Hex(5,:))/2
plot([A(1,1),B(1,1)],[A(1,2),B(1,2)],'--k')
text(0,0.1,'h_\perp')
axis off