clear
clc
% 参数设置
numParticles = 50;          % 粒子数量
numIterations = 50;        % 迭代次
w = 0.8;                    % 惯性权重
w_min = 0.4;                % 最小惯性权重
w_max = 0.9;                % 最大惯性权重
w_damp = 0.98;              % 惯性权重的衰减因子
cognitiveWeight = 1.5;      % 认知权重
socialWeight = 2.0;         % 社会权重
maxVelocity = 0.1;          % 最大速度
sigma = 1;                % 惩罚因子

% 约束设置
q_max=[pi/3 pi/3 pi/3 pi/3 pi/3 pi/3 pi/6];
qd_max=[1.5*pi 1.5*pi 1.5*pi 1.5*pi 1.5*pi 1.5*pi 1.5*pi];
qdd_max=[1.5*pi 1.5*pi 1.5*pi 1.5*pi 1.5*pi 1.5*pi 1.5*pi];

% 初始化粒子位置和速度
positions = rand(numParticles, 56); % 随机初始化位置
velocities = rand(numParticles, 56) * maxVelocity; % 随机初始化速度
fitness = zeros(numParticles, 1);
constraint_all = zeros(numParticles, 1);
constraint=zeros(1,7);


% 计算初始适应度
for i=1:numParticles
    for j = 1:7
        constraint(j)=max(non_q(positions(i,1+8*(j-1):8*j),q_max(j)),0)+3*max(non_qd(positions(i,1+8*(j-1):8*j),qd_max(j)),0)+5*max(non_qdd(positions(i,1+8*(j-1):8*j),qdd_max(j)),0);
    end
    constraint_all(i)=sum(constraint);
    fitness(i) = fun(positions(i,:))+ sigma*constraint_all(i);
end
personalBestPositions = positions;
personalBestFitness = fitness;

% 找到全局最优解和适应度
[globalBestFitness, globalBestIdx] = min(personalBestFitness);
globalBestPosition = personalBestPositions(globalBestIdx, :);
constraint_best = constraint_all(globalBestIdx);

% 迭代
for iter = 1:numIterations
    tic
    % 更新粒子位置和速度
    r1 = rand(numParticles, 56);
    r2 = rand(numParticles, 56);
    velocities = w' * velocities ...
        + cognitiveWeight * r1 .* (personalBestPositions - positions) ...
        + socialWeight * r2 .* (repmat(globalBestPosition, numParticles, 1) - positions);
    
    % 限制速度
    velocities = min(max(velocities, -maxVelocity), maxVelocity);
    
    % 更新位置
    positions = positions + velocities;
    
    % 修正超出约束的位置
    for i=1:56
        positions(:,i) = min(positions(:,i), 1); % 限制定义域内
        positions(:,i) = max(positions(:,i), -1);
    end
    
    % 计算适应度
    for i=1:numParticles
        for j = 1:7
            constraint(j)=max(non_q(positions(i,1+8*(j-1):8*j),q_max(j)),0)+3*max(non_qd(positions(i,1+8*(j-1):8*j),qd_max(j)),0)+5*max(non_qdd(positions(i,1+8*(j-1):8*j),qdd_max(j)),0);
        end
        constraint_all(i)=sum(constraint);
        fitness(i) = fun(positions(i,:))+ sigma*constraint_all(i);
    end
    % 更新个体最优解和适应度
    updateIdx = fitness < personalBestFitness;
    personalBestPositions(updateIdx, :) = positions(updateIdx, :);
    personalBestFitness(updateIdx) = fitness(updateIdx);
    
    % 更新全局最优解和适应度
    [iterBestFitness, iterBestIdx] = min(personalBestFitness);
    if iterBestFitness < globalBestFitness
        globalBestFitness = iterBestFitness;
        globalBestPosition = personalBestPositions(iterBestIdx, :);
        constraint_best = constraint_all(iterBestIdx);
    end
    ong(iter)=globalBestFitness;

    % 自适应调整惯性权重
    w = w_min + (w_max - w_min) * rand() + 0.25*randn();
    w = min(w, w_max);

    % 惯性权重衰减
    w = w * w_damp;

    % 输出每次迭代的结果
    fprintf('Iteration %d: Best Fitness = %.4f\n Constraintt_all: %.4f\n', iter, globalBestFitness,constraint_best);
    toc
    if iter>1
        if abs(ong(iter)-ong(iter-1))<1e-4
%              break;
        end
    end
end

% 输出最终结果
fprintf('\nFinal Result:\n');
globalBestPosition;
fprintf('Best Fitness: %.4f\n', globalBestFitness);
save('C:\Users\17881\Desktop\frank_pso\data\x.mat','globalBestPosition');
%% fitness
function [f]=fun(X)
ctrl_freq=1000;%频率步数
T=10;%总时间?
total_step=T*ctrl_freq; %总步数
dt=1/ctrl_freq;
t=zeros(1,total_step+1);
observe_matrix=zeros(7*(total_step+1),43);
w=2*pi/10;
%激励轨迹
a=zeros(7,5);
b=zeros(7,5);
q0=zeros(7,1);
x=zeros(7,8);
for i=1:7
    x(i,:)=X(1+8*(i-1):8*i);
    a(i,1:4)=x(i,1:4);
    b(i,1:3)=x(i,5:7);
    q0(i)=x(i,8);
    a(i,5)=-sum(a(i,1:4));%等式约束
    b(i,4)=-(2/27)*(144*b(i,1) + 63*b(i,2) + 32*b(i,3) - 150*q0(i)*w);
    b(i,5)=5/27*(45*b(i,1) + 18*b(i,2) + 7*b(i,3) - 48*q0(i)*w);
end
q=q0;
qd=zeros(7,1);
qdd=zeros(7,1);
for i=1:total_step+1
    t(i)=(i-1)*dt;
    for j=1:7
        for n=1:5
            q(j)=q(j)+a(j,n)/w/n*sinpi(w/pi*n*t(i))-b(j,n)/w/n*cospi(w/pi*n*t(i));
            qd(j)=qd(j)+a(j,n)*cospi(w/pi*n*t(i))+b(j,n)*sinpi(w/pi*n*t(i));
            qdd(j)=qdd(j)-a(j,n)*w*n*sinpi(w/pi*n*t(i))+b(j,n)*w*n*cospi(w/pi*n*t(i));
        end
    end
     observe_matrix(1+7*(i-1):7*i,:)=Yr_7(q(1),q(2),q(3),q(4),q(5),q(6),q(7),qd(1),qd(2),qd(3),qd(4),qd(5),qd(6),qd(7),qdd(1),qdd(2),qdd(3),qdd(4),qdd(5),qdd(6),qdd(7));
end
f=cond(observe_matrix,2)/1e15;
end
%% nonliner
%% q
function qmax=non_q(x,q)
a=zeros(1,5);
b=zeros(1,5);
a(1:4)=x(1:4);
b(1:3)=x(5:7);
q0=x(8);
w=2*pi/10;
a(5)=-sum(a(1:4));%等式约束
b(4)=-(2/27)*(144*b(1) + 63*b(2) + 32*b(3) - 150*q0*w);
b(5)=5/27*(45*b(1) + 18*b(2) + 7*b(3) - 48*q0*w);

qmax=0;
for i=1:5
    qmax=qmax+sqrt(a(i)^2+b(i)^2)/i/w;
end
qmax=qmax+abs(q0)-q;
end
%% qd
function qdmax=non_qd(x,q)
a=zeros(1,5);
b=zeros(1,5);
a(1:4)=x(1:4);
b(1:3)=x(5:7);
q0=x(8);
w=2*pi/10;
a(5)=-sum(a(1:4));%等式约束
b(4)=-(2/27)*(144*b(1) + 63*b(2) + 32*b(3) - 150*q0*w);
b(5)=5/27*(45*b(1) + 18*b(2) + 7*b(3) - 48*q0*w);


qdmax=0;
for i=1:5
    qdmax=qdmax+sqrt(a(i)^2+b(i)^2);
end
qdmax=qdmax-q;
end
%% qdd
function qddmax=non_qdd(x,q)
a=zeros(1,5);
b=zeros(1,5);
a(1:4)=x(1:4);
b(1:3)=x(5:7);
q0=x(8);
w=2*pi/10;
a(5)=-sum(a(1:4));%等式约束
b(4)=-(2/27)*(144*b(1) + 63*b(2) + 32*b(3) - 150*q0*w);
b(5)=5/27*(45*b(1) + 18*b(2) + 7*b(3) - 48*q0*w);

qddmax=0;
for i=1:5
    qddmax=qddmax+sqrt(a(i)^2+b(i)^2)*i*w;
end
qddmax=qddmax-q;
end