clear;
clc;
addpath([cd '\Function']);
% syms pi;
%% �غɲ���
m_end=1:0.5:5;m_end_z=0.1;
%%  ode45����
% �������
ctrl_freq=1000;%����Ƶ��
T=2;%��ģ��ʱ��
total_step=T*ctrl_freq; %���Ʋ������������Ƶ��Ϊ1000Hz���Ʋ�����1s�ķ������޸�1000�Σ�
dt=1/ctrl_freq;
%% ���ݴ洢
qe=zeros(7,total_step+1);%����λ�ã�����ѡ��Ϊĩ��ִ�������꣨q1,q2,q3,q4,q5��
qed=zeros(7,total_step+1);
% qedd=zeros(7,total_step+1);
q=zeros(7,total_step+1);%ʵ��λ�ã�����ѡ��Ϊĩ��ִ�������꣨q1,q2,q3,q4,q5��
qd=zeros(7,total_step+1);
qdd=zeros(7,total_step+1);
t=zeros(1,total_step+1);
%% �켣�滮(��ʱ����Ҫ)
q0=[0.3141
   -1.5618
   -0.2750
   -1.2527
    0.7200
    2.7172
    0.3627];
% q0=[0;0;0;0;0;0;0];
qd0=[0;0;0;0;0;0;0];
% ·����ֵ�滮
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
%% ��ʼ����
KK=zeros(5,9);
for h=1:5
    kp=15:0.5:20;
    k=zeros(1,9);
%     e=zeros(7*9,total_step+1);
    for j=1:9
        tau_action=zeros(7,total_step+1);
        step=zeros(1,total_step+1);
        p=zeros(9,2001);%ĩ��׼ȷ��
        q(:,1)=q0;%������ʼλ��
        qd(:,1)=qd0;%������ʼ�ٶ�
        options = odeset('AbsTol',5e-14,'Reltol',5e-14);
        tspan = [0 1/ctrl_freq];
        D=diag([7 9 6 8 5 4 3]);% ������ͬ�ؽڵ��������
        y0=[q(:,1);qd(:,1)];
        k1=kp(2*h);%λ������
        k2=2;%�ٶ�����
        y_last=y0;
        tic
        WB=waitbar(0,'please wait');

        for i=2:total_step+1
            %             ����һ�����Ʋ�
            [t_temp,y_temp]=ode45(@(t,y) myode(y_last,tau_action(:,i-1),m_end(j),m_end_z),tspan,y_last,options); %�������ΪһС�����Ʋ�
            y_last=y_temp(end,:).';
            q(:,i)=y_last(1:7,:);
            qd(:,i)=y_last(8:14,:);
            [dydt,M,C,G] = myode_output(y_last,tau_action(:,i-1),m_end(j),m_end_z);
            f_l=stribeck(y_last(8:14,:));
            % ������ѡ��
            % PD Control
            tau_action(:,i)=k1*D*(qe(:,i)-q(:,i))+k2*D*(-qd(:,i))+G+f_l;
            % tau_action(:,i)=k1*D*(qe(:,i)-q(:,i))+k2*D*(-qd(:,i));
            str=['������...',num2str(i/total_step*100),'%'];
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