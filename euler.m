clear;
clc;
digits(4)
%% 机械臂参数?
n=7;%轴数
m1=3.76;m2=4.64;m3=3.02;m4=3.45;m5=1.15;m6=1.07;m7=0.89;
m=[m1 m2 m3 m4 m5 m6 m7];
I=zeros(3,3,n);
I(:,:,1)=Ine(0.0664,0.0635,0.0096,0,0,-0.0005);
I(:,:,2)=Ine(0.0697,0.011,0.0666,-0.0001,0,0.0006);
I(:,:,3)=Ine(0.0209,0.0284,0.0149,0.0031,-0.0001,-0.0001);
I(:,:,4)=Ine(0.0225,0.0184,0.0327,0.0135,0.0001,0);
I(:,:,5)=Ine(0.0276,0.024,0.0052,0,0,-0.0023);
I(:,:,6)=Ine(0.0034,0.008,0.0086,-0.0008,0.0011,-0.0001);
I(:,:,7)=Ine(0.0019,0.0019,0.0015,0.0002,-0.0001,0.0002);

Ic=zeros(3,3,n);
rc=[0 0.0159 -0.0957;0.0001 -0.0797 -0.0184;-0.036 0.0128 -0.0528;-0.0414 -0.0478 -0.0142;-0.0001 0.0413 -0.109;-0.064 0.0094 -0.0232;0.007 0.007 -0.0317].';
for i=1:n
    Ic(:,:,i)=I(:,:,i)-m(i)*(rc(:,i).'*rc(:,i)*eye(3)-rc(:,i)*rc(:,i).');
end
% MDH建立机械臂?
g=9.806;
a=[0 0 0 8.3 -8.3 0 -8.7 0]*1e-2;%m
d=[33.3 0 31.5 0 38.3 0 10.8 0]*1e-2;
alpha=[0 -pi/2 pi/2 -pi/2 pi/2 -pi/2 pi/2 0];
q=sym('q',[1,n]);
qq=[q 0];
T = sym('T', [4,4,n+1]);
R=sym('R',[3,3,n+1]);
p=sym('p',[3,n+1]);
for i=1:n+1
    T(:,:,i)=expand(Tr(a(i),d(i),sym(alpha(i)),qq(i)));
    R(:,:,i)=T(1:3,1:3,i);
    p(:,i)=T(1:3,4,i);
end
%% 动力学信息
v=sym('v',[3,n]);
vd=sym('vd',[3,n]);
vc=sym('vcd',[3,n]);
vcd=sym('vcd',[3,n]);
w=sym('qd',[3,n]);
wd=sym('qdd',[3,n]);
qd=sym('qd',[1,n]);
qdd=sym('qdd',[1,n]);
fc=sym('fc',[3,n]);
tao=sym('tao',[3,n]);
f=sym('f',[3,n]);
torque=sym('t',[3,n]);

%初始值?
w0=[0;0;0];
w=[w0 w];
wd0=[0;0;0];
wd=[wd0 wd];
v0=[0;0;0];
v=[v0 v];
vd0=[0;0;g];
vd=[vd0 vd];
e=[0;0;1];
f_end=[0;0;0];
f=[f f_end];
torque_end=[0;0;0];
torque=[torque torque_end];
%% 外推
for i=1:n
    w(:,i+1)=R(:,:,i).'*w(:,i)+qd(i)*e;
    wd(:,i+1)=R(:,:,i).'*wd(:,i)+cross(R(:,:,i).'*w(:,i),qd(i)*e)+qdd(i)*[0;0;1];
    v(:,i+1)=R(:,:,i).'*(v(:,i)+cross(w(:,i),p(:,i)));
    vd(:,i+1)=R(:,:,i).'*(vd(:,i)+cross(wd(:,i),p(:,i))+cross(w(:,i),cross(w(:,i),p(:,i))));
    vcd(:,i)=vd(:,i+1)+cross(wd(:,i+1),rc(:,i))+cross(w(:,i+1),cross(w(:,i+1),rc(:,i)));
    vc(:,i)=v(:,i+1)+cross(w(:,i+1),rc(:,i));
    fc(:,i)=m(i)*vcd(:,i);
    tao(:,i)=Ic(:,:,i)*wd(:,i+1)+cross(w(:,i+1),Ic(:,:,i)*w(:,i+1));
end
%% 内推
for i=n:-1:1
    f(:,i)=R(:,:,i+1)*f(:,i+1)+fc(:,i);
    torque(:,i)=R(:,:,i+1)*torque(:,i+1)+tao(:,i)+cross(rc(:,i),fc(:,i))+cross(p(:,i+1),R(:,:,i+1)*f(:,i+1));
end
tor=e.'*torque(:,1:7);

%% 力矩计算
ctrl_freq=100;
T=5;
total_step=T*ctrl_freq;
dt=1/ctrl_freq;

Q=zeros(7,total_step+1);
Qd=zeros(7,total_step+1);
Qdd=zeros(7,total_step+1);
t=zeros(1,total_step+1);
n1=zeros(1,total_step+1);
n2=zeros(1,total_step+1);
n3=zeros(1,total_step+1);
n4=zeros(1,total_step+1);
n5=zeros(1,total_step+1);
n6=zeros(1,total_step+1);
n7=zeros(1,total_step+1);
tic
WB=waitbar(0,'please wait');
for i=1:total_step+1
    t(1,i)=(i-1)*dt;
    Q(:,i)=[Track(pi/5,5,t(1,i));Track(-pi/5,5,t(1,i));Track(pi/5,5,t(1,i));Track(-pi/5,5,t(1,i));Track(pi/5,5,t(1,i));Track(-pi/5,5,t(1,i));Track(pi/5,5,t(1,i))];
    Qd(:,i)=[Vol(pi/5,5,t(1,i));Vol(-pi/5,5,t(1,i));Vol(pi/5,5,t(1,i));Vol(-pi/5,5,t(1,i));Vol(pi/5,5,t(1,i));Vol(-pi/5,5,t(1,i));Vol(pi/5,5,t(1,i))];
    Qdd(:,i)=[Acl(pi/5,5,t(1,i));Acl(-pi/5,5,t(1,i));Acl(pi/5,5,t(1,i));Acl(-pi/5,5,t(1,i));Acl(pi/5,5,t(1,i));Acl(-pi/5,5,t(1,i));Acl(pi/5,5,t(1,i))];
    N=subs(tor,[q qd qdd],[Q(:,i).' Qd(:,i).' Qdd(:,i).']);
    n1(i)=N(1);n2(i)=N(2);n3(i)=N(3);n4(i)=N(4);n5(i)=N(5);n6(i)=N(6);n7(i)=N(7);
    str=['运行中...',num2str(i/total_step*100),'%'];
    waitbar(i/total_step,WB,str)
end
delete(WB);
toc
%% 连杆变换矩阵
function Tran = Tr(a,d,alpha,theta)
Tran=[cos(theta) -sin(theta) 0 a;sin(theta)*cos(alpha) cos(theta)*cos(alpha) -sin(alpha) -d*sin(alpha);sin(theta)*sin(alpha) cos(theta)*sin(alpha) cos(alpha) d*cos(alpha);0 0 0 1];
end
%%
function Inertia=Ine(Ix,Iy,Iz,Ixy,Ixz,Iyz)
Inertia=[Ix -Ixy -Ixz ;-Ixy Iy -Iyz ;-Ixz -Iyz Iz];
end
%%
function track=Track(theta,time,t)
a_0=0;
a_1=0;
a_2=3*theta/time^2;
a_3=-2*theta/time^3;
track=a_0+a_1*t+a_2*t^2+a_3*t^3;
end
%%
function vol=Vol(theta,time,t)
a_1=0;
a_2=3*theta/time^2;
a_3=-2*theta/time^3;
vol=a_1+2*a_2*t+3*a_3*t^2;
end
%%
function acl=Acl(theta,time,t)
a_2=3*theta/time^2;
a_3=-2*theta/time^3;
acl=2*a_2+6*a_3*t;
end