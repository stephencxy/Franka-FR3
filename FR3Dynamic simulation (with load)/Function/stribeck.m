function [f_l]=stribeck(q_d)

m1=3.75;m2=4.64;m3=3.02;m4=3.45;m5=1.15;m6=1.07;m7=0.89;
deta_m=m2-m7;
m=[m1 m2 m3 m4 m5 m6 m7];

sigma_0=zeros(1,7);
sigma_1=zeros(1,7);
Fc_0=zeros(1,7);
Fc_1=zeros(1,7);
Fs_0=zeros(1,7);
Fs_1=zeros(1,7);
vs_0=zeros(1,7);
vs_1=zeros(1,7);

sigma0_min=0.06988;
sigma1_min=0.07609;
Fc0_min=0.0978;
Fc1_min=0.0890;
Fs0_min=0.06483;
Fs1_min=0.05764;
vs0_max=+0.30838;
vs1_max=0.43791;

deta_sigma1=-0.06988+0.1210;
deta_sigma0=-0.07609+0.2345;
deta_Fc1=-0.0978+0.2117;
deta_Fc0=-0.08902+0.1393;
deta_Fs1=-0.06483+0.4635;
deta_Fs0=-0.05764+0.4531;
deta_vs1=-0.2530+0.43791;
deta_vs0=-0.1896+0.30838;

for i=1:7
    sigma_0(i)=(m(i)-min(m))*deta_sigma0/deta_m+sigma0_min;
    sigma_1(i)=(m(i)-min(m))*deta_sigma1/deta_m+sigma1_min;
    Fc_0(i)=(m(i)-min(m))*deta_Fc0/deta_m+Fc0_min;
    Fc_1(i)=(m(i)-min(m))*deta_Fc1/deta_m+Fc1_min;
    Fs_0(i)=(m(i)-min(m))*deta_Fs0/deta_m+Fs0_min;
    Fs_1(i)=(m(i)-min(m))*deta_Fs1/deta_m+Fs1_min;
    vs_0(i)=-(m(i)-min(m))*deta_vs0/deta_m+vs0_max;
    vs_1(i)=-(m(i)-min(m))*deta_vs1/deta_m+vs1_max;
end

f_l=zeros(7,1);
for i=1:7
    if q_d(i)<0
        f_l(i) = lugref_ss(q_d(i), Fc_0(i) ,Fs_0(i) ,vs_0(i), sigma_0(i));
    end
    f_l(i) = lugref_ss(q_d(i), Fc_1(i) ,Fs_1(i) ,vs_1(i), sigma_1(i));
end
%%
function Fss = lugref_ss(v, Fc, Fs, vs, sigma_2)
    r = -(v/vs).^2;
    Fss = Fc*sign(v)+ (Fs - Fc) * exp(r) .* sign(v) + sigma_2 * v;      
end
end