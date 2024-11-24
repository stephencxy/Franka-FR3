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
T=3;%总模拟时长
total_step=T*ctrl_freq; %控制步数（如果控制频率为1000Hz控制步数在1s的仿真中修改1000次）
dt=1/ctrl_freq;
num=2;

%% 目标位置
P_0=[-8.9 0  119].';
P0=[-45 0 20].';
P1=[-30 15 4].';
P2=[-30 -15 36].';
P3=[-60 -15 36].';
P4=[-60 15 4].';
P0_4=[P_0 P0 P1 P2 P3 P4 P0 P1 P2 P3 P4 P0];
l_i=zeros(1,num*5);
l_all=0;
s_l=0;

for i=1:num*5
    q_actual=q(:,(i+1)*(total_step+1));%P_0
    T_actual=subs(T_07,Q,q_actual.');
    X=T_actual(1,4);
    Y=T_actual(2,4);
    Z=T_actual(3,4);
    l_i(i)=sqrt((P0_4(1,i+2)-X)^2+(P0_4(2,i+2)-Y)^2+(P0_4(3,i+2)-Z)^2);
    l_all=l_i(i)+l_all;
end

l_ave=l_all/num/5;
for i=1:num*5
    s_l=s_l+(l_i(i)-l_ave)^2;
end
s_l=sqrt(s_l/(num*5-1));
rp_l=l_ave+3*s_l;
%% MDH连杆矩阵
function Tran = Tr(a,d,alpha,theta)
theta=theta/pi;
alpha=alpha/pi;
Tran=[cospi(theta) -sinpi(theta) 0 a;sinpi(theta)*cospi(alpha) cospi(theta)*cospi(alpha) -sinpi(alpha) -d*sinpi(alpha);sinpi(theta)*sinpi(alpha) cospi(theta)*sinpi(alpha) cospi(alpha) d*cospi(alpha);0 0 0 1];
end