clear;
clc;
addpath([cd '\Function']);
% syms pi
%% 载荷参数
m_end=2;m_end_z=0.1;
%%  ode45仿真
% 仿真参数
ctrl_freq=100;%控制频率
T=5;%总模拟时长
total_step=T*ctrl_freq; %控制步数（如果控制频率为1000Hz控制步数在1s的仿真中修改1000次）
dt=1/ctrl_freq;
%% 数据存储
qe=zeros(7,total_step+1);%期望位置，坐标选择为末端执行器坐标（q1,q2,q3,q4,q5）
qed=zeros(7,total_step+1);
qedd=zeros(7,total_step+1);
q=zeros(7,total_step+1);%实际位置，坐标选择为末端执行器坐标（q1,q2,q3,q4,q5）
qd=zeros(7,total_step+1);
qdd=zeros(7,total_step+1);
t=zeros(1,total_step+1);
%% 轨迹规划(暂时不需要)
q0=[0;0;0;0;0;0;0];
qd0=[0;0;0;0;0;0;0];
% 路径插值规划
for i=1:total_step+1
    t(1,i)=(i-1)*dt;
    qe(:,i)=[Track(pi/5,5,t(1,i));Track(-pi/5,5,t(1,i));Track(pi/5,5,t(1,i));Track(-pi/5,5,t(1,i));Track(pi/5,5,t(1,i));Track(-pi/5,5,t(1,i));Track(pi/5,5,t(1,i))];
    qed(:,i)=[Vol(pi/5,5,t(1,i));Vol(-pi/5,5,t(1,i));Vol(pi/5,5,t(1,i));Vol(-pi/5,5,t(1,i));Vol(pi/5,5,t(1,i));Vol(-pi/5,5,t(1,i));Vol(pi/5,5,t(1,i))];
    qedd(:,i)=[Acl(pi/5,5,t(1,i));Acl(-pi/5,5,t(1,i));Acl(pi/5,5,t(1,i));Acl(-pi/5,5,t(1,i));Acl(pi/5,5,t(1,i));Acl(-pi/5,5,t(1,i));Acl(pi/5,5,t(1,i))];
end
%% 初始条件
tau_action=zeros(7,total_step+1);
step=zeros(1,total_step+1);
q(:,1)=q0;%机构初始位置
qd(:,1)=qd0;%机构初始速度
options = odeset('AbsTol',5e-14,'Reltol',5e-14);
tspan = [0 1/ctrl_freq];
D=diag([7 9 6 8 5 4 3]);% 调整不同关节的增益比例
y0=[q(:,1);qd(:,1)];
k1=14;%位置增益
k2=1.6;%速度增益
y_last=y0;
% tao=zeros(7,total_step+1);
tic
WB=waitbar(0,'please wait'); 
for i=2:total_step+1
    %     仿真一个控制步
    [t_temp,y_temp]=ode45(@(t,y) myode(y_last,tau_action(:,i-1),m_end,m_end_z),tspan,y_last,options); %仿真过程为一小个控制步
    y_last=y_temp(end,:).';
    q(:,i)=y_last(1:7,:);
    qd(:,i)=y_last(8:14,:);
    [dydt,M,C,G] = myode_output(y_last,tau_action(:,i-1),m_end,m_end_z);
    %     tao(:,i)=M*dydt(8:14)+C*dydt(1:7)+G;
    %% 控制律选择
    % PD Control
    %     tau_action(:,i)=k1*D*(qe(:,i)-q(:,i))+k2*D*(qed(:,i)-qd(:,i));
    tau_action(:,i)=k1/k2*M*(qe(:,i)-q(:,i))+C*qd(:,i)+G+M*(qed(:,i)-qd(:,i))+M*qedd(:,i);
    %%
    str=['运行中...',num2str(i/total_step*100),'%'];
    waitbar(i/total_step,WB,str)
end
delete(WB);
toc
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