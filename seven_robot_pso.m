clear;
clc;
addpath([cd '\Function']);
% syms pi
%%  ode45����
% �������
ctrl_freq=1000;%����Ƶ��
T=20;%��ģ��ʱ��
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
x=load('C:\Users\17881\Desktop\frank_pso\data\x.mat','globalBestPosition');
x=x.globalBestPosition;
%% �켣�滮(��ʱ����Ҫ)
q0=[0;0;0;0;0;0;0];
qd0=[0;0;0;0;0;0;0];
% ·����ֵ�滮
for i=1:total_step+1
    t(1,i)=(i-1)*dt;
   for j=1:7
       [qe(j,i) ,qed(j,i) ,qedd(j,i)]=Tra(x(1+8*(j-1):8*j),t(i));
   end
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
k1=60;%λ������
k2=6;%�ٶ�����
y_last=y0;
% tao=zeros(7,total_step+1);
tic
WB=waitbar(0,'please wait'); 
for i=2:total_step+1
    %     ����һ�����Ʋ�
    [t_temp,y_temp]=ode45(@(t,y) myode(y_last,tau_action(:,i-1)),tspan,y_last,options); %�������ΪһС�����Ʋ�
    y_last=y_temp(end,:).';
    q(:,i)=y_last(1:7,:);
    qd(:,i)=y_last(8:14,:);
    [dydt,M,C,G] = myode_output(y_last,tau_action(:,i-1));
    qdd(:,i)=dydt(8:14);
    %% ������ѡ��
    % PD Control
    %     tau_action(:,i)=k1*D*(qe(:,i)-q(:,i))+k2*D*(qed(:,i)-qd(:,i));
    tau_action(:,i)=M*(qedd(:,i)+k1*D*(qe(:,i)-q(:,i))+k2*D*(qed(:,i)-qd(:,i)))+C*qed(:,i)+G;
    %%
    str=['������...',num2str(i/total_step*100),'%'];
    waitbar(i/total_step,WB,str)
end
delete(WB);
toc
qs=q(:,10001:20001);
qds=qd(:,10001:20001);
qdds=qdd(:,10001:20001);
taus=tau_action(:,10001:20001);
save('C:\Users\17881\Desktop\frank_pso\data\q.mat','qs')
save('C:\Users\17881\Desktop\frank_pso\data\qd.mat','qds')
save('C:\Users\17881\Desktop\frank_pso\data\qdd.mat','qdds')
save('C:\Users\17881\Desktop\frank_pso\data\tau.mat','taus')
%% �����켣
function [q ,qd ,qdd]=Tra(x,t)
w=2*pi/10;
%�����켣
a=zeros(1,5);
b=zeros(1,5);
a(1:4)=x(1:4);
b(1:3)=x(5:7);
q0=x(8);
a(5)=-sum(a(1:4));%��ʽԼ��
b(4)=-(2/27)*(144*b(1) + 63*b(2) + 32*b(3) - 150*q0*w);
b(5)=5/27*(45*b(1) + 18*b(2) + 7*b(3) - 48*q0*w);
q=q0;
qd=0;
qdd=0;
for n=1:5
    q=q+(a(n)/w/n*sin(w*n*t)-b(n)/w/n*cos(w*n*t));
    qd=qd+(a(n)*cos(w*n*t)+b(n)*sin(w*n*t));
    qdd=qdd-(a(n)*w*n*sin(w*n*t)-b(n)*w*n*cos(w*n*t));
end
end