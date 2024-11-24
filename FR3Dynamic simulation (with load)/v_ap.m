% clear 
% clc
digits(5)
% q0=[0.5413;-0.9199;-0.0508;-1.89;0.1233;2.8069;-0.3938];
% P1=[-30 15 4].';
qe=[-0.2311
   -2.0299
   -0.0628
   -1.0510
    0.7551
    3.0594
    0.9567];
p=zeros(9,2001);
t=zeros(1,2001);
for i=1:2001
    t(1,i)=(i-1)*0.001;
end
q=zeros(63,2001);
for j=1:9
    tic
    for i=1:2001
        q(7*(j-1)+1:7*j,i)= qe-e(7*(j-1)+1:7*j,i);
        Q= qe-e(7*(j-1)+1:7*j,i);
        P=Trans_07(Q(1),Q(2),Q(3),Q(4),Q(5),Q(6),Q(7));
        p(j,i)=sqrt((P(1,4)+30)^2+(P(2,4)-15)^2+(P(3,4)-4)^2);
    end
    toc
end
k=zeros(1,9);
for i=1:9
    x=fun(p(i,:));
    k(i)=x;
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