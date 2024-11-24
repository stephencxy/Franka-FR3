 clc;clear;
% digits(5);
% syms pi
syms q1 q2 q3 q4 q5 q6 q7;
%% MDH建立机械臂
a=[0 0 0 8.3 -8.3 0 -8.9];%m
d=[33.3 0 31.5 0 38.3 0 15.9];
alpha=[0 -pi/2 pi/2 -pi/2 pi/2 -pi/2 pi/2];
q=[q1 q2 q3 q4 q5 q6 q7];
T = sym('T', [4,28]);%定义符号空集
T=T-T;
for i=1:7
    T(:,4*i-3:4*i)=vpa(Tr(a(i),d(i),alpha(i),q(i)));
end
T_01=T(:,1:4);
T_12=T(:,5:8);
T_23=T(:,9:12);
T_34=T(:,13:16);
T_45=T(:,17:20);
T_56=T(:,21:24);
T_67=T(:,25:28);
T_07=T_01*T_12*T_23*T_34*T_45*T_56*T_67;

%%  
ctrl_freq=1000;%控制频率
T=3;%总模拟时长
total_step=T*ctrl_freq; %控制步数（如果控制频率为1000Hz控制步数在1s的仿真中修改1000次）
dt=1/ctrl_freq;
%% 计算st
num=11;
t=zeros(1,num*(total_step+1));
for i=1:num*(total_step+1)
    t(i)=(i-1)*dt;
end
tic
P_0=[-8.9 0  119].';
P0=[-45 0 20].';
P1=[-30 15 4].';
P2=[-30 -15 36].';
P3=[-60 -15 36].';
P4=[-60 15 4].';
P0_4=[P_0 P0 P1 P2 P3 P4 P0 P1 P2 P3 P4 P0];


total=4*(total_step+1);
P=zeros(4,num*total);
for i=1:num
    [P(1:4,1+(i-1)*total:i*total),x]=cal(P0_4(1:3,i),P0_4(1:3,i+1),total_step,dt);
end
%% 计算逆运动学
G=1e-4;
iter=100;
% num=xlsread('C:\bishe\Tk151_A.xlsx');
% q0=[num(1,end) num(2,end) num(3,end) num(4,end) num(5,end) num(6,end) num(7,end)];
q0=[0 0 0 0 0 0 0];
qc=zeros(7,num*(total_step+1)+1);
qc(:,1)=q0.';
T_0=subs(T_07,q,q0);
fq_0=trans(T_0);
fq_now=fq_0;
tic
WB=waitbar(0,'please wait');
for k=1:num*total_step+num
    T_end=P(1:4,4*k-3:4*k);
    fq_end=trans(T_end);
    i=1;
    while i<=iter
        e=fq_end-fq_now;
        g=One_norm(e);
        if g<=G
            break;
        end
        jcb=JACOB(qc(1,k),qc(2,k),qc(3,k),qc(4,k),qc(5,k),qc(6,k),qc(7,k));
        deta_q=jcb.'*pinv(jcb*jcb.')*e;
        qc(:,k)=qc(:,k)+deta_q;
        T_now=subs(T_07,q,qc(:,k).');
        fq_now=trans(T_now);
        i=i+1;
    end
    if k==1
        for m=1:7
            j=1;
            while j<=100
                tt=double(qc(m,k));
                if tt<-3.14159
                    qc(m,k)=qc(m,k)+2*3.14159;
                end
                if tt>3.14159
                    qc(m,k)=qc(m,k)-2*3.14159;
                end
                if tt<=3.14159 && tt>=-3.14159
                    break;
                end
            end
        end
    end
    qc(:,k+1)=qc(:,k);
    %% 
    str=['运行中...',num2str(k/total_step/num*100),'%'];
    waitbar(k/total_step/num,WB,str)
end
delete(WB);
toc
%% 1norm
function norm=One_norm(M)
norm=0;
for i=1:6
    norm=norm+M(i)^2;
end
norm=sqrt(norm);
end
%%
function X=trans(T)
X=zeros(6,1);
X(1)=T(1,4);
X(2)=T(2,4);
X(3)=T(3,4);
X(4)=atan2(T(3,2),T(3,3));
X(5)=atan2(-T(3,1),sqrt(T(1,1)^2+T(2,1)^2));
X(6)=atan2(T(2,1),T(1,1));
end
%%
function Tran = Tr(a,d,alpha,theta)
theta=theta/pi;
alpha=alpha/pi;
Tran=[cospi(theta) -sinpi(theta) 0 a;sinpi(theta)*cospi(alpha) cospi(theta)*cospi(alpha) -sinpi(alpha) -d*sinpi(alpha);sinpi(theta)*sinpi(alpha) cospi(theta)*sinpi(alpha) cospi(alpha) d*cospi(alpha);0 0 0 1];
end
%%
function [P,deta_P]=cal(P0,P1,total_step,dt)
tic
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
toc
end