clear;
clc;
addpath([cd '\Function']);
% syms pi;
%% 载荷参数
m_end=1:0.5:5;m_end_z=0.1;
%%  ode45仿真
% 仿真参数
ctrl_freq=1000;%控制频率
T=2;%总模拟时长
total_step=T*ctrl_freq; %控制步数（如果控制频率为1000Hz控制步数在1s的仿真中修改1000次）
dt=1/ctrl_freq;
%% 数据存储
qe=zeros(7,total_step+1);%期望位置，坐标选择为末端执行器坐标（q1,q2,q3,q4,q5）
qed=zeros(7,total_step+1);
% qedd=zeros(7,total_step+1);
q=zeros(7,total_step+1);%实际位置，坐标选择为末端执行器坐标（q1,q2,q3,q4,q5）
qd=zeros(7,total_step+1);
qdd=zeros(7,total_step+1);
t=zeros(1,total_step+1);
%% 轨迹规划(暂时不需要)
q0=[0.3141
   -1.5618
   -0.2750
   -1.2527
    0.7200
    2.7172
    0.3627];
% q0=[0;0;0;0;0;0;0];
qd0=[0;0;0;0;0;0;0];
% 路径插值规划
for i=1:total_step+1
    t(1,i)=(i-1)*dt;
    qe(:,i)=[-0.2311
   -2.0299
   -0.0628
   -1.0510
    0.7551
    3.0594
    0.9567];
%     qe(:,i)=[pi/4;pi/4;pi/4;pi/4;pi/4;pi/4;pi/4];
end
%% 初始条件
KK=zeros(5,9);
for h=1:5
    kp=15:0.5:20;
    k=zeros(1,9);
%     e=zeros(7*9,total_step+1);
    for j=1:9
        tau_action=zeros(7,total_step+1);
        step=zeros(1,total_step+1);
        p=zeros(9,2001);%末端准确度
        q(:,1)=q0;%机构初始位置
        qd(:,1)=qd0;%机构初始速度
        options = odeset('AbsTol',5e-14,'Reltol',5e-14);
        tspan = [0 1/ctrl_freq];
        D=diag([7 9 6 8 5 4 3]);% 调整不同关节的增益比例
        y0=[q(:,1);qd(:,1)];
        k1=kp(2*h);%位置增益
        k2=2;%速度增益
        y_last=y0;
        tic
        WB=waitbar(0,'please wait');

        for i=2:total_step+1
            %             仿真一个控制步
            [t_temp,y_temp]=ode45(@(t,y) myode(y_last,tau_action(:,i-1),m_end(j),m_end_z),tspan,y_last,options); %仿真过程为一小个控制步
            y_last=y_temp(end,:).';
            q(:,i)=y_last(1:7,:);
            qd(:,i)=y_last(8:14,:);
            [dydt,M,C,G] = myode_output(y_last,tau_action(:,i-1),m_end(j),m_end_z);
            f_l=stribeck(y_last(8:14,:));
            % 控制律选择
            % PD Control
            tau_action(:,i)=k1*D*(qe(:,i)-q(:,i))+k2*D*(-qd(:,i))+G+f_l;
            % tau_action(:,i)=k1*D*(qe(:,i)-q(:,i))+k2*D*(-qd(:,i));
            str=['运行中...',num2str(i/total_step*100),'%'];
            waitbar(i/total_step,WB,str)
        end
        delete(WB);
        toc
        for ii=1:2001
            Q=q(:,ii);
            P=Trans_07(Q(1),Q(2),Q(3),Q(4),Q(5),Q(6),Q(7));
            p(j,ii)=sqrt((P(1,4)+30)^2+(P(2,4)-15)^2+(P(3,4)-4)^2);
        end
        x=fun(p(j,:));
        k(j)=x;
        %         e(7*j-6:7*j,:)=qe-q;
        %         j
        % e=qe-q;
    end
    KK(h,:)=k;
end
%%
function eta=fun(P)
n=0;m=0;
p1=diff(P);
for i=2:1999
    if p1(i)<=0 && p1(i+1)>=0
        n=i+1;
        break;
    end
end
for j=n:1999
    if p1(j)<=0 && p1(j+1)>=0
        m=j+1;
        break;
    end
end

eta=log(P(n)/P(m));
end