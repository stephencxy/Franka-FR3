clear;
clc;
addpath([cd '\Function']);
% syms pi
%% �غɲ���
m_end=2;m_end_z=0.1;
%%  ode45����
% �������
ctrl_freq=100;%����Ƶ��
T=5;%��ģ��ʱ��
total_step=T*ctrl_freq; %���Ʋ������������Ƶ��Ϊ1000Hz���Ʋ�����1s�ķ������޸�1000�Σ�
dt=1/ctrl_freq;
%% ���ݴ洢
qe=zeros(7,total_step+1);%����λ�ã�����ѡ��Ϊĩ��ִ�������꣨q1,q2,q3,q4,q5��
qed=zeros(7,total_step+1);
qedd=zeros(7,total_step+1);
q=zeros(7,total_step+1);%ʵ��λ�ã�����ѡ��Ϊĩ��ִ�������꣨q1,q2,q3,q4,q5��
qd=zeros(7,total_step+1);
qdd=zeros(7,total_step+1);
t=zeros(1,total_step+1);
%% �켣�滮(��ʱ����Ҫ)
q0=[0;0;0;0;0;0;0];
qd0=[0;0;0;0;0;0;0];
% ·����ֵ�滮
for i=1:total_step+1
    t(1,i)=(i-1)*dt;
    qe(:,i)=[Track(pi/5,5,t(1,i));Track(-pi/5,5,t(1,i));Track(pi/5,5,t(1,i));Track(-pi/5,5,t(1,i));Track(pi/5,5,t(1,i));Track(-pi/5,5,t(1,i));Track(pi/5,5,t(1,i))];
    qed(:,i)=[Vol(pi/5,5,t(1,i));Vol(-pi/5,5,t(1,i));Vol(pi/5,5,t(1,i));Vol(-pi/5,5,t(1,i));Vol(pi/5,5,t(1,i));Vol(-pi/5,5,t(1,i));Vol(pi/5,5,t(1,i))];
    qedd(:,i)=[Acl(pi/5,5,t(1,i));Acl(-pi/5,5,t(1,i));Acl(pi/5,5,t(1,i));Acl(-pi/5,5,t(1,i));Acl(pi/5,5,t(1,i));Acl(-pi/5,5,t(1,i));Acl(pi/5,5,t(1,i))];
end
%% ��ʼ����
tau_action=zeros(7,total_step+1);
step=zeros(1,total_step+1);
q(:,1)=q0;%������ʼλ��
qd(:,1)=qd0;%������ʼ�ٶ�
options = odeset('AbsTol',5e-14,'Reltol',5e-14);
tspan = [0 1/ctrl_freq];
D=diag([7 9 6 8 5 4 3]);% ������ͬ�ؽڵ��������
y0=[q(:,1);qd(:,1)];
k1=14;%λ������
k2=1.6;%�ٶ�����
y_last=y0;
% tao=zeros(7,total_step+1);
tic
WB=waitbar(0,'please wait'); 
for i=2:total_step+1
    %     ����һ�����Ʋ�
    [t_temp,y_temp]=ode45(@(t,y) myode(y_last,tau_action(:,i-1),m_end,m_end_z),tspan,y_last,options); %�������ΪһС�����Ʋ�
    y_last=y_temp(end,:).';
    q(:,i)=y_last(1:7,:);
    qd(:,i)=y_last(8:14,:);
    [dydt,M,C,G] = myode_output(y_last,tau_action(:,i-1),m_end,m_end_z);
    %     tao(:,i)=M*dydt(8:14)+C*dydt(1:7)+G;
    %% ������ѡ��
    % PD Control
    %     tau_action(:,i)=k1*D*(qe(:,i)-q(:,i))+k2*D*(qed(:,i)-qd(:,i));
    tau_action(:,i)=k1/k2*M*(qe(:,i)-q(:,i))+C*qd(:,i)+G+M*(qed(:,i)-qd(:,i))+M*qedd(:,i);
    %%
    str=['������...',num2str(i/total_step*100),'%'];
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