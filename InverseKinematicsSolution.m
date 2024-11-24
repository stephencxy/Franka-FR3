clc;clear;
digits(4);
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
    T(:,4*i-3:4*i)=vpa(Tr(a(i),d(i),sym(alpha(i)),q(i)));
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
P1=[1 0 0 0;0 1 0 0;0 0 1 0;-40 0 20 1].';
P2=[1 0 0 0;0 1 0 0;0 0 1 0;-20 15 0 1].';
P3=[1 0 0 0;0 1 0 0;0 0 1 0;-20 -15 40 1].';
P4=[1 0 0 0;0 1 0 0;0 0 1 0;-60 -15 40 1].';
P5=[1 0 0 0;0 1 0 0;0 0 1 0;-60 15 0 1].';
P=[P1 P2 P3 P4 P5];
%% 计算逆运动学
G=1e-2;
iter=100;
% num=xlsread('C:\bishe\Tk_a.xlsx');
% q0=[num(1,end) num(2,end) num(3,end) num(4,end) num(5,end) num(6,end) num(7,end)];
q0=[0 0 0 0 0 0 0];
qc=zeros(7,5);
T_0=subs(T_07,q,q0);
fq_0=trans(T_0);
fq_now=fq_0;
for k=1:5
    T_end=P(1:4,4*k-3:4*k);
    fq_end=trans(T_end);
    i=1;
    tic
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
    for m=1:7
        j=1;
        while j<=100
            tt=double(qc(m,k));
            if tt<-3.14159-1
                qc(m,k)=qc(m,k)+2*3.14159;
            end
            if tt>3.14159+1
                qc(m,k)=qc(m,k)-2*3.14159;
            end
            if tt<=3.14159+1 && tt>=-3.14159-1
                break;
            end
        end
    end
qc(:,k+1)=qc(:,k);
end
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
Tran=[cos(theta) -sin(theta) 0 a;sin(theta)*cos(alpha) cos(theta)*cos(alpha) -sin(alpha) -d*sin(alpha);sin(theta)*sin(alpha) cos(theta)*sin(alpha) cos(alpha) d*cos(alpha);0 0 0 1];
end