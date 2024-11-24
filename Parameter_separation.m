clear;clc
digits 4
% 参数定义
%% 机械臂参数
n=7;%轴数
g=9.8;
% I=sym('I',[3,3,n]);%惯性矩阵
% for i=1:n
%     I(:,:,i)=symmetry(I(:,:,i));
% end
% rc=sym('rc',[3,n]);%连杆质心
% m=sym('m',[1,n]);%连杆质量
I=zeros(3,3,n);
I(:,:,1)=Ine(0.0664,0.0635,0.0096,0,0,-0.0005);
I(:,:,2)=Ine(0.0697,0.011,0.0666,-0.0001,0,0.0006);
I(:,:,3)=Ine(0.0209,0.0284,0.0149,0.0031,-0.0001,-0.0001);
I(:,:,4)=Ine(0.0225,0.0184,0.0327,0.0135,0.0001,0);
I(:,:,5)=Ine(0.0276,0.024,0.0052,0,0,-0.0023);
I(:,:,6)=Ine(0.0034,0.008,0.0086,-0.0008,0.0011,-0.0001);
I(:,:,7)=Ine(0.0019,0.0019,0.0015,0.0002,-0.0001,0.0002);
rc=[0 0.0159 -0.0957;0.0001 -0.0797 -0.0184;-0.036 0.0128 -0.0528;-0.0414 -0.0478 -0.0142;...
    -0.0001 0.0413 -0.109;-0.064 0.0094 -0.0232;0.007 0.007 -0.0317].';
m1=3.76;m2=4.64;m3=3.02;m4=3.45;m5=1.15;m6=1.07;m7=0.89;
m=[m1 m2 m3 m4 m5 m6 m7];

phi=sym('phi',[10*n,1]);%惯性参数集
for i=1:n
    phi(1+10*(i-1):10*i)=[m(i);m(i)*rc(:,i);I(1,1,i);I(1,2,i);I(1,3,i);I(2,2,i);I(2,3,i);I(3,3,i)];
end

% MDH建立机械臂
a=[0 0 0 8.3 -8.3 0 -8.7 0]*1e-2;%a=sym('a',[1,n]);
d=[33.3 0 31.5 0 38.3 0 10.8 0]*1e-2;
alpha=[0 -pi/2 pi/2 -pi/2 pi/2 -pi/2 pi/2 0];
q=sym('q',[1,n]);
q=[q 0];
T = sym('T', [4,4,n+1]);%定义符号空集
R=sym('R',[3,3,n+1]);%旋转矩阵
p=sym('p',[3,n+1]);%位移矢量
for i=1:n+1
    T(:,:,i)=expand(Tr(a(i),d(i),sym(alpha(i)),q(i)));
    R(:,:,i)=T(1:3,1:3,i);
    p(:,i)=T(1:3,4,i);
end
%% 动力学信息
v=sym('v',[3,n]);
vd=sym('vd',[3,n]);
w=sym('qd',[3,n]);
wd=sym('qdd',[3,n]);
qd=sym('qd',[1,n]);
qdd=sym('qdd',[1,n]);

%初始值
w0=[0;0;0];
w=[w0 w];
wd0=[0;0;0];
wd=[wd0 wd];
v0=[0;0;0];
v=[v0 v];
vd0=[0;0;g];
vd=[vd0 vd];
e=[0;0;1];
%% 外推
for i=1:n
    w(:,i+1)=R(:,:,i).'*w(:,i)+qd(i)*e;
    wd(:,i+1)=R(:,:,i).'*wd(:,i)+cross(R(:,:,i).'*wd(:,i),qd(i)*e)+qdd(i)*[0;0;1];
    v(:,i+1)=R(:,:,i).'*(v(:,i)+cross(w(:,i),p(:,i)));
    vd(:,i+1)=R(:,:,i).'*(vd(:,i)+cross(wd(:,i),p(:,i))+cross(w(:,i),cross(w(:,i),p(:,i))));
end                                               
1
%% A矩阵
A=sym('A',[6,10,n]);
for i=1:n
    A(1:3,1,i)=vd(:,i+1);A(1:3,2:4,i)=wedge(wd(:,i+1))+wedge(w(:,i+1))*wedge(w(:,i+1));A(1:3,5:10,i)=zeros(3,6);
    A(4:6,1,i)=zeros(3,1);A(4:6,2:4,i)=-wedge(vd(:,i+1));A(4:6,5:10,i)=kw(wd(:,i+1))+wedge(w(:,i+1))*kw(w(:,i+1));
end
2
%% Q矩阵
Q=sym('Q',[6,6,n+1]);
for i=1:n+1
    Q(:,:,i)=[R(:,:,i) zeros(3,3);wedge(p(:,i))*R(:,:,i) R(:,:,i)];
end
3
%% U矩阵
U=sym('U',[6*n,10*n]);
for i=1:n
    for j=i:n
        if i==j
            U(1+6*(i-1):6*i,1+10*(j-1):10*j)=A(:,:,i);
        else
            Q_all=eye(6);
            for k=i:j-1
            Q_all=Q_all*Q(:,:,k+1);
            end
            U(1+6*(i-1):6*i,1+10*(j-1):10*j)=Q_all*A(:,:,j);
        end
    end
end
for i=n:-1:2
    for j=i-1:-1:1
        U(1+6*(i-1):6*i,1+10*(j-1):10*j)=zeros(6,10);
    end
end
4
%% Y矩阵
sigma=[0,0,0,0,0,1];
Y=sym('Y',[1*n,10*n]);
for i=1:n
    for j=1:n
        Y(i,1+10*(j-1):10*j)=sigma*U(1+6*(i-1):6*i,1+10*(j-1):10*j);
    end
end
5
% Tor=Y*phi;
%% 最小惯性参数集
phi_min=sym('phi_min',[7*(n-1)+1,1]);
phi_re=phi;
for i=n:-1:2
    phi_re(1+10*(i-2):10*(i-1))=R1(phi_re(1+10*(i-2):10*(i-1)),phi_re(1+10*(i-1):10*i),a(i),sym(alpha(i)),d(i));
end
phi_min(1)=phi_re(10);
for i=2:n
    phi_min(2+7*(i-2):1+7*(i-1))=[phi_re(2+10*(i-1)),phi_re(3+10*(i-1)),phi_re(5+10*(i-1))-phi_re(8+10*(i-1)),phi_re(6+10*(i-1)),phi_re(7+10*(i-1)),phi_re(9+10*(i-1)),phi_re(10+10*(i-1))].';
end
phi_min=vpa(phi_min);
6
%% 观测矩阵Yr
Yr=sym('Yr',[n,1+7*(n-1)]);
Yr(:,1)=Y(:,10);
for i=2:n
    Yr(:,2+7*(i-2):1+7*(i-1))=[Y(:,2+10*(i-1)) Y(:,3+10*(i-1)) Y(:,5+10*(i-1)) Y(:,6+10*(i-1)) Y(:,7+10*(i-1)) Y(:,9+10*(i-1)) Y(:,10+10*(i-1))];
end
7
% matlabFunction(Yr,'File','Yr','Vars',[q(1:2) qd qdd g]);
    %% 对称矩阵
function inertia=symmetry(I)
I(2,1)=I(1,2);
I(3,1)=I(1,3);
I(3,2)=I(2,3);
inertia=I;
end
%% 连杆变换矩阵
function Tran = Tr(a,d,alpha,theta)
Tran=[cos(theta) -sin(theta) 0 a;sin(theta)*cos(alpha) cos(theta)*cos(alpha) -sin(alpha) -d*sin(alpha);sin(theta)*sin(alpha) cos(theta)*sin(alpha) cos(alpha) d*cos(alpha);0 0 0 1];
end
%% s（w）
function sw=wedge(w)
sw=[0 -w(3) w(2);w(3) 0 -w(1);-w(2) w(1) 0];
end
%% k(w)
function kw=kw(w)
kw=[w(1) w(2) w(3) 0 0 0;0 w(1) 0 w(2) w(3) 0;0 0 w(1) 0 w(2) w(3)];
end
%% 重组惯性参数R1
function p=R1(phi_1,phi_2,a,alpha,d)
M=phi_1(1)+phi_2(1);
MX=phi_1(2)+a*phi_2(1);
MY=phi_1(3)-sin(alpha)*(phi_2(4)+d*phi_2(1));
MZ=phi_1(4)+cos(alpha)*(phi_2(4)+d*phi_2(1));
XX=phi_1(5)+phi_2(8)+2*d*phi_2(4)+d*d*phi_2(1);
XY=phi_1(6)+a*sin(alpha)*(phi_2(4)+d*phi_2(1));
XZ=phi_1(7)-a*cos(alpha)*(phi_2(4)+d*phi_2(1));
YY=phi_1(8)+(cos(alpha))^2*phi_2(8)+2*d*(cos(alpha))^2*phi_2(4)+(d^2*(cos(alpha))^2+a^2)*phi_2(1);
YZ=phi_1(9)+cos(alpha)*sin(alpha)*(phi_2(8)+2*d*phi_2(4)+d^2*phi_2(1));
ZZ=phi_1(10)+(sin(alpha))^2*phi_2(8)+2*d*(sin(alpha))^2*phi_2(4)+(d^2*(sin(alpha))^2+a^2)*phi_2(1);
p=[M,MX,MY,MZ,XX,XY,XZ,YY,YZ,ZZ];
end
%%
function Inertia=Ine(Ix,Iy,Iz,Ixy,Ixz,Iyz)
Inertia=[Ix -Ixy -Ixz ;-Ixy Iy -Iyz ;-Ixz -Iyz Iz];
end