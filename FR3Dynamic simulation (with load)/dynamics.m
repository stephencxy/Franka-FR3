clear
clc
digits(4)
% syms pi;
syms m_end m_end_x m_end_y m_end_z
syms q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d
assume([q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z],'real');
m_end=2;
m_end_z=0.1;
q=[q1;q2;q3;q4;q5;q6;q7];
q_d=[q1_d;q2_d;q3_d;q4_d;q5_d;q6_d;q7_d];
%% ��е������
m1=3.75;m2=4.64;m3=3.02;m4=3.45;m5=1.15;m6=1.07;m7=0.89;g=9.8;
%% MDH������е��
a=[0 0 8.3 0 -8.3 0 -8.9]*1e-2;%m
d=[19.2 0 31.5 0 38.3 0 15.9]*1e-2;
alpha=[0 -pi/2 pi/2 -pi/2 pi/2 -pi/2 pi/2];
T = sym('T', [4,28]);%������ſռ�
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
T_78=[1 0 0 0;0 1 0 0;0 0 1 m_end_z/2;0 0 0 1];
%% �ؽ�1��Link1���ܱ�ʾ
I1=I(0.0033,0.0062,0.0602,0,0,-0.0003,m1,0,0.0154,-0.0957);
T_1=0;
for j=1:7
    for k=1:7
        T_1=T_1+vpa(trace(diff(T_01,q(j))*I1*diff(T_01.',q(k))*q_d(j)*q_d(k)));
    end
end
%% �ؽ�2��Link2���ܱ�ʾ
I2=I(0.004,0.0626,0.007,-0.0001,0,0.0004,m2,0.0001,-0.0797,-0.0179);
T_02=T_01*T_12;
T_2=0;
for j=1:7
    for k=1:7
        T_2=T_2+vpa(trace(diff(T_02,q(j))*I2*diff(T_02.',q(k))*q_d(j)*q_d(k)));
    end
end
%% �ؽ�3��Link3���ܱ�ʾ
I3=I(0.0138,0.0037,0.0172,0,0.013,0,m3,-0.0465,0.0123,-0.0528);
T_03=T_01*T_12*T_23;
T_3=0;
for j=1:7
    for k=1:7
        T_3=T_3+vpa(trace(diff(T_03,q(j))*I3*diff(T_03.',q(k))*q_d(j)*q_d(k)));
    end
end
%% �ؽ�4��Link4���ܱ�ʾ
I4=I(0.0143,0.0185,0.004,0.0135,0.0001,0,m4,-0.0414,-0.0478,-0.0142);
T_04=T_01*T_12*T_23*T_34;
T_4=0;
for j=1:7
    for k=1:7
        T_4=T_4+vpa(trace(diff(T_04,q(j))*I4*diff(T_04.',q(k))*q_d(j)*q_d(k)));
    end
end
%% �ؽ�5��Link5���ܱ�ʾ
I5=I(0.0008,0.0045,0.0232,0.0002,-0.0003,-0.0024,m5,0.0029,0.0424,-0.1088);
T_05=T_01*T_12*T_23*T_34*T_45;
T_5=0;
for j=1:7
    for k=1:7
        T_5=T_5+vpa(trace(diff(T_05,q(j))*I5*diff(T_05.',q(k))*q_d(j)*q_d(k)));
    end
end
%% �ؽ�6��Link6���ܱ�ʾ
I6=I(0.007,0.002,0.0014,-0.0008,0.0011,-0.0002,m6,-0.064,0.0094,-0.0232);
T_06=T_01*T_12*T_23*T_34*T_45*T_56;
T_6=0;
for j=1:7
    for k=1:7
        T_6=T_6+vpa(trace(diff(T_06,q(j))*I6*diff(T_06.',q(k))*q_d(j)*q_d(k)));
    end
end
%% �ؽ�7��Link7���ܱ�ʾ
I7=I(0.0005,0.001,0.0011,0,0,0.0002,m7,0.0001,-0.0099,-0.0317);
T_07=T_01*T_12*T_23*T_34*T_45*T_56*T_67;
T_7=0;
for j=1:7
    for k=1:7
        T_7=T_7+vpa(trace(diff(T_07,q(j))*I7*diff(T_07.',q(k))*q_d(j)*q_d(k)));
    end
end
%% ĩ��ִ��������
Ix=m_end*m_end_z^2/6;
Iy=m_end*m_end_z^2/6;
Iz=m_end*m_end_z^2/6;
I8=I(Ix,Iy,Iz,0,0,0,m_end,0,0,0);
T_08=T_01*T_12*T_23*T_34*T_45*T_56*T_67*T_78;
T_end=0;
for j=1:7
    for k=1:7
        T_end=T_end+vpa(trace(diff(T_08,q(j))*I8*diff(T_08.',q(k))*q_d(j)*q_d(k)));
    end
end
%% ϵͳ�ܶ���
T_total=(T_1+T_2+T_3+T_4+T_5+T_6+T_7+T_end)/2;
%% ���������
M = sym('M', [7,7]);%������ſռ�
M=M-M;
for i=1:7
    for j=1:7
        temp=diff(T_total,q_d(i));
        M(i,j)=diff(temp,q_d(j));
    end
end
% tic
% matlabFunction(M(1,:),'File','Calc_M1','Vars',[q1 q2 q3 q4 q5 q6 q7 m_end m_end_z]);
% toc
% tic
% matlabFunction(M(2,:),'File','Calc_M2','Vars',[q1 q2 q3 q4 q5 q6 q7 m_end m_end_z]);
% toc
% tic
% matlabFunction(M(3,:),'File','Calc_M3','Vars',[q1 q2 q3 q4 q5 q6 q7 m_end m_end_z]);
% toc
% tic
% matlabFunction(M(4,:),'File','Calc_M4','Vars',[q1 q2 q3 q4 q5 q6 q7 m_end m_end_z]);
% toc
% tic
% matlabFunction(M(5,:),'File','Calc_M5','Vars',[q1 q2 q3 q4 q5 q6 q7 m_end m_end_z]);
% toc
% tic
% matlabFunction(M(6,:),'File','Calc_M6','Vars',[q1 q2 q3 q4 q5 q6 q7 m_end m_end_z]);
% toc
% tic
% matlabFunction(M(7,:),'File','Calc_M7','Vars',[q1 q2 q3 q4 q5 q6 q7 m_end m_end_z]);
% toc
%% ���������
% C = sym('C', [7,7]);%������ſռ�
% C=C-C;
% for i=1:7
%    for j=1:7
%       for k=1:7
%         C(i,j)=C(i,j)+(diff(M(i,j),q(k))-1/2*diff(M(j,k),q(i)))*q_d(k);%(diff(M(i,j),q(k,1))-1/2*diff(M(j,k),q(i,1)))*qd(k,1);
%       end
%    end
% end
% tic
% matlabFunction(C(1,1),'File','Calc_C11','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(1,2),'File','Calc_C12','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(1,3),'File','Calc_C13','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(1,4),'File','Calc_C14','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(1,5),'File','Calc_C15','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(1,6),'File','Calc_C16','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(1,7),'File','Calc_C17','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(2,1),'File','Calc_C21','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(2,2),'File','Calc_C22','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(2,3),'File','Calc_C23','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(2,4),'File','Calc_C24','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(2,5),'File','Calc_C25','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(2,6),'File','Calc_C26','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(2,7),'File','Calc_C27','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% matlabFunction(C(3,1),'File','Calc_C31','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(3,2),'File','Calc_C32','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(3,3),'File','Calc_C33','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(3,4),'File','Calc_C34','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(3,5),'File','Calc_C35','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(3,6),'File','Calc_C36','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(3,7),'File','Calc_C37','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(4,1),'File','Calc_C41','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(4,2),'File','Calc_C42','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(4,3),'File','Calc_C43','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(4,4),'File','Calc_C44','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(4,5),'File','Calc_C45','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(4,6),'File','Calc_C46','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(4,7),'File','Calc_C47','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(5,1),'File','Calc_C51','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(5,2),'File','Calc_C52','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(5,3),'File','Calc_C53','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(5,4),'File','Calc_C54','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(5,5),'File','Calc_C55','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(5,6),'File','Calc_C56','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(5,7),'File','Calc_C57','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(6,1),'File','Calc_C61','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(6,2),'File','Calc_C62','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(6,3),'File','Calc_C63','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(6,4),'File','Calc_C64','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(6,5),'File','Calc_C65','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(6,6),'File','Calc_C66','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(6,7),'File','Calc_C67','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(7,1),'File','Calc_C71','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(7,2),'File','Calc_C72','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(7,3),'File','Calc_C73','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(7,4),'File','Calc_C74','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(7,5),'File','Calc_C75','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(7,6),'File','Calc_C76','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
% tic
% matlabFunction(C(7,7),'File','Calc_C77','Vars',[q1 q2 q3 q4 q5 q6 q7 q1_d q2_d q3_d q4_d q5_d q6_d q7_d m_end m_end_z]);
% toc
%% ϵͳ����U����o0-x0-y0��Ϊ��������
% Link1�͹ؽ�1����
C_L1=T_01*[0 0.0154 -0.0957 1].';
U_L1=m1*g*C_L1(3);
% Link2�͹ؽ�2����
C_L2=T_01*T_12*[0.0001 -0.07974 -0.01791 1].';
U_L2=m2*g*C_L2(3);
% Link3�͹ؽ�3����
C_L3=T_01*T_12*T_23*[-0.04655 0.01232 -0.05283 1].';
U_L3=m3*g*C_L3(3);
% Link4�͹ؽ�4����
C_L4=T_01*T_12*T_23*T_34*[-0.0414 -0.04783 -0.01418 1].';
U_L4=m4*g*C_L4(3);
% Link5�͹ؽ�5����
C_L5=T_01*T_12*T_23*T_34*T_45*[0.00292 0.04235 -0.10879 1].';
U_L5=m5*g*C_L5(3);
% Link6�͹ؽ�6����
C_L6=T_01*T_12*T_23*T_34*T_45*T_56*[-0.06399 0.00944 -0.02321 1].';
U_L6=m6*g*C_L6(3);
% Link7�͹ؽ�7����
C_L7=T_01*T_12*T_23*T_34*T_45*T_56*T_67*[0.00005 -0.00989 -0.03166 1].';
U_L7=m7*g*C_L7(3);
% ĩ��ִ��������
C_end=T_08*[0 0 0 1].';
U_end=m_end*g*C_end(3);
% ϵͳ������
U=(U_L1+U_L2+U_L3+U_L4+U_L5+U_L6+U_L7+U_end);
%% ����������
G = sym('G', [7,1]);%������ſռ�
G=G-G;
for i=1:7
    G(i)=diff(U,q(i));
end
% matlabFunction(G,'File','Calc_G','Vars',[q1 q2 q3 q4 q5 q6 q7 m_end m_end_z]);
%%
function Tran = Tr(a,d,alpha,theta)
Tran=[cos(theta) -sin(theta) 0 a;sin(theta)*cos(alpha) cos(theta)*cos(alpha) -sin(alpha) -d*sin(alpha);sin(theta)*sin(alpha) cos(theta)*sin(alpha) cos(alpha) d*cos(alpha);0 0 0 1];
end
%%
function Inertia=I(Ixx,Iyy,Izz,Ixy,Ixz,Iyz,m,x,y,z)
Inertia=[Ixx Ixy Ixz m*x;Ixy Iyy Iyz m*y;Ixz Iyz Izz m*z;m*x m*y m*z m];
end