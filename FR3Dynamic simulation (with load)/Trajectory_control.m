clear;
clc;
addpath([cd '\Function']);
% digits(6);
%% 载荷参数
m_end=2;m_end_z=0.1;
%% MDH建立连杆矩阵
syms q1 q2 q3 q4 q5 q6 q7;
a=[0 0 0 8.3 -8.3 0 -8.9];%m
d=[33.3 0 31.5 0 38.3 0 15.9];
alpha=[0 -pi/2 pi/2 -pi/2 pi/2 -pi/2 pi/2];
Q=[q1 q2 q3 q4 q5 q6 q7];
T = sym('T', [4,28]);%定义符号空集
T=T-T;
for i=1:7
    T(:,4*i-3:4*i)=vpa(Tr(a(i),d(i),alpha(i),Q(i)));
end
T_01=T(:,1:4);
T_12=T(:,5:8);
T_23=T(:,9:12);
T_34=T(:,13:16);
T_45=T(:,17:20);
T_56=T(:,21:24);
T_67=T(:,25:28);
T_07=T_01*T_12*T_23*T_34*T_45*T_56*T_67;
%%  ode45仿真
% 仿真参数
ctrl_freq=1000;%控制频率
Time=3;%总模拟时长
total_step=Time*ctrl_freq; %控制步数（如果控制频率为1000Hz控制步数在1s的仿真中修改1000次）
dt=1/ctrl_freq;
%% 轨迹规划
q0=[0;0;0;0;0;0;0];
qd0=[0;0;0;0;0;0;0];

%% 计算st
P0=[-8.9 0  119].';
% P0=[-45 0 20].';
% P1=[-24.891 -54.142 64.66].';
P1=[-30 15 4].';
P2=[-30 -15 36].';
P3=[-60 -15 36].';
P4=[-60 15 4].';
P0_4=[P0 P1 P2 P3 P4];%目标位置

cycle_num=2;%循环次数 
total=4*(total_step+1);
P=zeros(4,cycle_num*total);
for i=1:cycle_num
    [P(1:4,1+(i-1)*total:i*total),x]=cal(P0_4(1:3,i),P0_4(1:3,i+1),total_step,dt);
end

%% 数据存储
qe=zeros(7,cycle_num*(total_step+1));%期望位置
qed=zeros(7,cycle_num*(total_step+1));
qedd=zeros(7,cycle_num*(total_step+1));
q=zeros(7,cycle_num*(total_step+1));%实际位置
qd=zeros(7,cycle_num*(total_step+1));
qdd=zeros(7,cycle_num*(total_step+1));
t=zeros(1,cycle_num*(total_step+1));
for i=1:cycle_num*(total_step+1)
    t(i)=(i-1)*dt;
end
%% 计算逆运动学
G=1e-4;
iter=100;
qe(:,1)=q0;
T_0=subs(T_07,Q,q0.');
fq_0=trans(T_0);
fq_now=fq_0;
%% 关节角度
for k=1:cycle_num*(total_step+1)
    T_end=P(1:4,4*k-3:4*k);
    fq_end=trans(T_end);
    i=1;
    while i<=iter
        e=fq_end-fq_now;
        g=One_norm(e);
        if g<=G
            break;
        end
        jcb=JACOB(qe(1,k),qe(2,k),qe(3,k),qe(4,k),qe(5,k),qe(6,k),qe(7,k));
        deta_q=jcb.'*pinv(jcb*jcb.')*e;
        qe(:,k)=qe(:,k)+deta_q;
        T_now=subs(T_07,Q,qe(:,k).');
        fq_now=trans(T_now);
        i=i+1;
    end
    if k==1
        for m=1:7
            j=1;
            while j<=100
                tt=double(qe(m,k));
                if tt<-3.14159
                    qe(m,k)=qe(m,k)+2*3.14159;
                end
                if tt>3.14159
                    qe(m,k)=qe(m,k)-2*3.14159;
                end
                if tt<=3.14159 && tt>=-3.14159
                    break;
                end
            end
        end
    end
    qe(:,k+1)=qe(:,k);
end
%% 关节速度
qed(:,:)=diff(qe,1,2)/dt;
qedd(:,1:cycle_num*(total_step+1)-1)=diff(qed,1,2)/dt;
qedd(:,cycle_num*(total_step+1))=qedd(:,cycle_num*(total_step+1)-1);
qe=qe(:,1:cycle_num*(total_step+1));
%% 初始条件
tau_action=zeros(7,cycle_num*(total_step+1));
step=zeros(1,cycle_num*(total_step+1));
q(:,1)=q0;%机构初始位置
qd(:,1)=qd0;%机构初始速度
options = odeset('AbsTol',5e-14,'Reltol',5e-14);
tspan = [0 1/ctrl_freq];
D=diag([7 9 6 8 5 4 3]);% 调整不同关节的增益比例
y0=[q(:,1);qd(:,1)];
k1=15;%位置增益
k2=1;%速度增益
y_last=y0;
tic
WB=waitbar(0,'please wait'); 
for i=2:cycle_num*(total_step+1)
%     仿真一个控制步
    [t_temp,y_temp]=ode45(@(t,y) myode(y_last,tau_action(:,i-1),m_end,m_end_z),tspan,y_last,options); %仿真过程为一小个控制步
    y_last=y_temp(end,:).';
    q(:,i)=y_last(1:7,:);
    qd(:,i)=y_last(8:14,:);
    [dydt,M,C,G] = myode_output(y_last,tau_action(:,i-1),m_end,m_end_z);
%% 控制律选择
% PD Control 
%     tau_action(:,i)=k1*D*(qe(:,i)-q(:,i))+k2*D*(qed(:,i)-qd(:,i));
    tau_action(:,i)=k1/k2*M*(qe(:,i)-q(:,i))+C*qd(:,i)+G+M*(qed(:,i)-qd(:,i))+M*qedd(:,i);
%% 
    str=['运行中...',num2str(i/total_step/cycle_num*100),'%'];
    waitbar(i/total_step/cycle_num,WB,str)
end
delete(WB);
toc
%% MDH连杆矩阵
function Tran = Tr(a,d,alpha,theta)
theta=theta/pi;
alpha=alpha/pi;
Tran=[cospi(theta) -sinpi(theta) 0 a;sinpi(theta)*cospi(alpha) cospi(theta)*cospi(alpha) -sinpi(alpha) -d*sinpi(alpha);sinpi(theta)*sinpi(alpha) cospi(theta)*sinpi(alpha) cospi(alpha) d*cospi(alpha);0 0 0 1];
end
%% 1norm
function norm=One_norm(M)
norm=0;
for i=1:6
    norm=norm+M(i)^2;
end
norm=sqrt(norm);
end
%% 连杆矩阵转换笛卡尔坐标
function X=trans(T)
X=zeros(6,1);
X(1)=T(1,4);
X(2)=T(2,4);
X(3)=T(3,4);
X(4)=atan2(T(3,2),T(3,3));
X(5)=atan2(-T(3,1),sqrt(T(1,1)^2+T(2,1)^2));
X(6)=atan2(T(2,1),T(1,1));
end
%% 梯形加减速
function [P,deta_P]=cal(P0,P1,total_step,dt)

st=zeros(1,total_step+1);
vt=zeros(1,total_step+1);
at=zeros(1,total_step+1);

deta_P=P1-P0;
L=sqrt(deta_P(1)^2+deta_P(2)^2+deta_P(3)^2);

t0=0;
t1=3;
Ta=0.75;
Vv=L/(t1-Ta);

t=zeros(1,total_step+1);
for i=1:total_step+1
    t(1,i)=(i-1)*dt;
end
for i=1:total_step+1
    if(t0<=t(i)&&t(i)<t0+Ta)
        st(i)=Vv/(2*Ta)*(t(i)-t0)^2;
        vt(i)=Vv/(Ta)*(t(i)-t0);
        at(i)=Vv/Ta;
    end
    if(t0+Ta<=t(i)&&t(i)<t1-Ta)
        st(i)=Vv*(t(i)-t0-Ta/2);
        vt(i)=Vv;
        at(i)=0;
    end
    if(t1-Ta<=t(i)&&t(i)<=t1)
        st(i)=L-Vv/(2*Ta)*(t1-t(i))^2;
        vt(i)=Vv/(Ta)*(t1-t(i));
        at(i)=-Vv/Ta;
    end
end
%% 计算插补点
p=zeros(3,total_step+1);
for i=1:total_step+1
    p(:,i)=P0+deta_P*st(i)/L;
end
P=zeros(4,4*(total_step+1));
for i=1:total_step+1
    P(:,4*i-3:4*i)=[1 0 0 p(1,i);0 1 0 p(2,i);0 0 1 p(3,i);0 0 0 1];
end
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