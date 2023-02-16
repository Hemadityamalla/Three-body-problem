function main
%All the units are changed
global G;
global x1;
global x2;
global x3;
global y1;
global y2;
global y3;
global vx1;
global vx2;
global vx3;
global vy1;
global vy2;
global vy3;
global m;
m=[1;1;1];
% m=[3.0028e-6;1;0];
G=4*pi*pi;%Gravitational constant
delt=0.01;
tf=10.0;
N=fix(tf/delt);
t=linspace(0.0,tf,N);
%Intial conditions
%Body-1
x1=zeros(N,1);
x1(1,1)=3.3030197;
% x1(1,1)=1.0;
y1=zeros(N,1);
y1(1,1)=-0.82771837;
% y1(1,1)=0.0;
vx1=zeros(N,1);
vx1(1,1)=1.587433767;
% vx1(1,1)=0.0;
vy1=zeros(N,1);
vy1(1,1)=1.47221479;
% vy1(1,1)=6.1682;
%Body-2
x2=zeros(N,1);
x2(1,1)=-3.3030197;
% x2(1,1)=0.0;
y2=zeros(N,1);
y2(1,1)=0.82771837;
% y2(1,1)=0.0;
vx2=zeros(N,1);
vx2(1,1)=1.587433767;
% vx2(1,1)=0.0;
vy2=zeros(N,1);
vy2(1,1)=1.47221479;
% vy2(1,1)=0.0;
%Body-3
x3=zeros(N,1);
x3(1,1)=0.0;
y3=zeros(N,1);
y3(1,1)=0.0;
vx3=zeros(N,1);
vx3(1,1)=-3.174867535;
vy3=zeros(N,1);
vy3(1,1)=-2.94442961;
E1=zeros(N,1);
E2=zeros(N,1);
E3=zeros(N,1);
TE=zeros(N,1);


%Implementation of RK3 
for i=1:N-1
    
    %Determination of k1 which is a 12 X 1 array
    k1=[fx1(i);fy1(i);vx1(i);vy1(i);fx2(i);fy2(i);vx2(i);vy2(i);fx3(i);fy3(i);vx3(i);vy3(i)];
    %Determination of k2
    trvx1=vx1(i,1)+((1/3)*delt*k1(1,1));
    trvy1=vy1(i,1)+((1/3)*delt*k1(2,1));
    trx1=x1(i,1)+((1/3)*delt*k1(3,1));
    try1=y1(i,1)+((1/3)*delt*k1(4,1));
    
    trvx2=vx2(i,1)+((1/3)*delt*k1(5,1));
    trvy2=vy2(i,1)+((1/3)*delt*k1(6,1));
    trx2=x2(i,1)+((1/3)*delt*k1(7,1));
    try2=y2(i,1)+((1/3)*delt*k1(8,1));
    
    trvx3=vx3(i,1)+((1/3)*delt*k1(9,1));
    trvy3=vy3(i,1)+((1/3)*delt*k1(10,1));
    trx3=x3(i,1)+((1/3)*delt*k1(11,1));
    try3=y3(i,1)+((1/3)*delt*k1(12,1));
    
    trfx1=(-G)*((m(2,1)*(trx1-trx2)/r(trx1,trx2,try1,try2))+(m(3,1)*(trx1-trx3)/r(trx1,trx3,try1,try3)));
    trfy1=(-G)*((m(2,1)*(try1-try2)/r(trx1,trx2,try1,try2))+(m(3,1)*(try1-try3)/r(trx1,trx3,try1,try3)));
    trfx2=(-G)*((m(1,1)*(trx2-trx1)/r(trx1,trx2,try1,try2))+(m(3,1)*(trx2-trx3)/r(trx2,trx3,try2,try3)));
    trfy2=(-G)*((m(1,1)*(try2-try1)/r(trx1,trx2,try1,try2))+(m(3,1)*(try2-try3)/r(trx2,trx3,try2,try3)));
    trfx3=(-G)*((m(1,1)*(trx3-trx1)/r(trx1,trx3,try1,try3))+(m(2,1)*(trx3-trx2)/r(trx2,trx3,try2,try3)));
    trfy3=(-G)*((m(1,1)*(try3-try1)/r(trx1,trx3,try1,try3))+(m(2,1)*(try3-try2)/r(trx2,trx3,try2,try3)));
    
    k2=[trfx1;trfy1;trvx1;trvy1;trfx2;trfy2;trvx2;trvy2;trfx3;trfy3;trvx3;trvy3];
    
    %Determination of k3
    
    ttrvx1=vx1(i,1)+((2/3)*delt*k2(1,1));
    ttrvy1=vy1(i,1)+((2/3)*delt*k2(2,1));
    ttrx1=x1(i,1)+((2/3)*delt*k2(3,1));
    ttry1=y1(i,1)+((2/3)*delt*k2(4,1));
    
    ttrvx2=vx2(i,1)+((2/3)*delt*k2(5,1));
    ttrvy2=vy2(i,1)+((2/3)*delt*k2(6,1));
    ttrx2=x2(i,1)+((2/3)*delt*k2(7,1));
    ttry2=y2(i,1)+((2/3)*delt*k2(8,1));
    
    ttrvx3=vx3(i,1)+((2/3)*delt*k2(9,1));
    ttrvy3=vy3(i,1)+((2/3)*delt*k2(10,1));
    ttrx3=x3(i,1)+((2/3)*delt*k2(11,1));
    ttry3=y3(i,1)+((2/3)*delt*k2(12,1));
    
    ttrfx1=(-G)*((m(2,1)*(ttrx1-ttrx2)/r(ttrx1,ttrx2,ttry1,ttry2))+(m(3,1)*(ttrx1-ttrx3)/r(ttrx1,ttrx3,ttry1,ttry3)));
    ttrfy1=(-G)*((m(2,1)*(ttry1-ttry2)/r(ttrx1,ttrx2,ttry1,ttry2))+(m(3,1)*(ttry1-ttry3)/r(ttrx1,ttrx3,ttry1,ttry3)));
    ttrfx2=(-G)*((m(1,1)*(ttrx2-ttrx1)/r(ttrx1,ttrx2,ttry1,ttry2))+(m(3,1)*(ttrx2-ttrx3)/r(ttrx2,ttrx3,ttry2,ttry3)));
    ttrfy2=(-G)*((m(1,1)*(ttry2-ttry1)/r(ttrx1,ttrx2,ttry1,ttry2))+(m(3,1)*(ttry2-ttry3)/r(ttrx2,ttrx3,ttry2,ttry3)));
    ttrfx3=(-G)*((m(1,1)*(ttrx3-ttrx1)/r(ttrx1,ttrx3,ttry1,ttry3))+(m(2,1)*(ttrx3-ttrx2)/r(ttrx2,ttrx3,ttry2,ttry3)));
    ttrfy3=(-G)*((m(1,1)*(ttry3-ttry1)/r(ttrx1,ttrx3,ttry1,ttry3))+(m(2,1)*(ttry3-ttry2)/r(ttrx2,ttrx3,ttry2,ttry3)));
    
    k3=[ttrfx1;ttrfy1;ttrvx1;ttrvy1;ttrfx2;ttrfy2;ttrvx2;ttrvy2;ttrfx3;ttrfy3;ttrvx3;ttrvy3];
    
    %Final step
    vx1(i+1,1)=vx1(i,1)+(0.25*delt*(k1(1,1)+(3*k2(1,1))));
    vy1(i+1,1)=vy1(i,1)+(0.25*delt*(k1(2,1)+(3*k2(2,1))));
    x1(i+1,1)=x1(i,1)+(0.25*delt*(k1(3,1)+(3*k2(3,1))));
    y1(i+1,1)=y1(i,1)+(0.25*delt*(k1(4,1)+(3*k2(4,1))));
    vx2(i+1,1)=vx2(i,1)+(0.25*delt*(k1(5,1)+(3*k2(5,1))));
    vy2(i+1,1)=vy2(i,1)+(0.25*delt*(k1(6,1)+(3*k2(6,1))));
    x2(i+1,1)=x2(i,1)+(0.25*delt*(k1(7,1)+(3*k2(7,1))));
    y2(i+1,1)=y2(i,1)+(0.25*delt*(k1(8,1)+(3*k2(8,1))));
    vx3(i+1,1)=vx3(i,1)+(0.25*delt*(k1(9,1)+(3*k2(9,1))));
    vy3(i+1,1)=vy3(i,1)+(0.25*delt*(k1(10,1)+(3*k2(10,1))));
    x3(i+1,1)=x3(i,1)+(0.25*delt*(k1(11,1)+(3*k2(11,1))));
    y3(i+1,1)=y3(i,1)+(0.25*delt*(k1(12,1)+(3*k2(12,1))));
    %Calculation of energy
    E1(i,1)=(0.5*m(1,1)*vx1(i,1)*vx1(i,1))+(0.5*m(1,1)*vy1(i,1)*vy1(i,1))+(((-G)*m(1,1)*m(2,1))/r(x1(i,1),x2(i,1),y1(i,1),y2(i,1)))+(((-G)*m(1,1)*m(3,1))/r(x1(i,1),x3(i,1),y1(i,1),y3(i,1)));
    E2(i,1)=(0.5*m(2,1)*vx2(i,1)*vx2(i,1))+(0.5*m(2,1)*vy2(i,1)*vy2(i,1))+(((-G)*m(1,1)*m(2,1))/r(x1(i,1),x2(i,1),y1(i,1),y2(i,1)))+(((-G)*m(2,1)*m(3,1))/r(x2(i,1),x3(i,1),y2(i,1),y3(i,1)));
    E3(i,1)=(0.5*m(3,1)*vx3(i,1)*vx3(i,1))+(0.5*m(3,1)*vy3(i,1)*vy3(i,1))+(((-G)*m(3,1)*m(2,1))/r(x3(i,1),x2(i,1),y3(i,1),y2(i,1)))+(((-G)*m(1,1)*m(3,1))/r(x1(i,1),x3(i,1),y1(i,1),y3(i,1)));
%     TE=E1(i,1)+E2(i,1)+E3(i,1);
end
TE(1:N-1,1)=E1(1:N-1,1)+E2(1:N-1,1)+E3(1:N-1,1);
figure(1)
plot(t(1:N-1),TE(1:N-1,1));
for i=1:N
    figure(2)
plot(x1(i),y1(i),'b*',x2(i),y2(i),'r*',x3(i),y3(i),'g*');
pause(0.000001);
hold on;

end








end


function y=r(a,b,c,d)

y=(((a-b)^2)+((c-d)^2))^1.5;

end

function y=fx1(a)
global G;
global x1;
global x2;
global x3;
global y1;
global y2;
global y3;
global m;

y=-G*((m(2,1)*(x1(a,1)-x2(a,1))/r(x1(a,1),x2(a,1),y1(a,1),y2(a,1)))+(m(3,1)*(x1(a,1)-x3(a,1))/r(x1(a,1),x3(a,1),y1(a,1),y3(a,1))));


end

function y=fy1(a)
global G;
global x1;
global x2;
global x3;
global y1;
global y2;
global y3;
global m;
y=-G*((m(2,1)*(y1(a,1)-y2(a,1))/r(x1(a,1),x2(a,1),y1(a,1),y2(a,1)))+(m(3,1)*(y1(a,1)-y3(a,1))/r(x1(a,1),x3(a,1),y1(a,1),y3(a,1))));


end

function y=fx2(a)
global G;
global x1;
global x2;
global x3;
global y1;
global y2;
global y3;
global m;


y=-G*((m(1,1)*(x2(a,1)-x1(a,1))/r(x1(a,1),x2(a,1),y1(a,1),y2(a,1)))+(m(3,1)*(x2(a,1)-x3(a,1))/r(x2(a,1),x3(a,1),y2(a,1),y3(a,1))));


end
function y=fy2(a)
global G;
global x1;
global x2;
global x3;
global y1;
global y2;
global y3;
global m;
y=-G*((m(1,1)*(y2(a,1)-y1(a,1))/r(x1(a,1),x2(a,1),y1(a,1),y2(a,1)))+(m(3,1)*(y2(a,1)-y3(a,1))/r(x2(a,1),x3(a,1),y2(a,1),y3(a,1))));


end
function y=fx3(a)
global G;
global x1;
global x2;
global x3;
global y1;
global y2;
global y3;
global m;

y=-G*((m(1,1)*(x3(a,1)-x1(a,1))/r(x3(a,1),x1(a,1),y3(a,1),y1(a,1)))+(m(2,1)*(x3(a,1)-x2(a,1))/r(x2(a,1),x3(a,1),y2(a,1),y3(a,1))));


end
function y=fy3(a)
global G;
global x1;
global x2;
global x3;
global y1;
global y2;
global y3;
global m;
y=-G*((m(1,1)*(y3(a,1)-y1(a,1))/r(x1(a,1),x3(a,1),y1(a,1),y3(a,1)))+(m(2,1)*(y3(a,1)-y2(a,1))/r(x2(a,1),x3(a,1),y2(a,1),y3(a,1))));


end

