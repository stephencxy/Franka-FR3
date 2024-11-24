y=1:0.5:5;
% x=15:0.5:20;%kp
x=1.5:0.05:2;%kd
[x,y]=meshgrid(x,y);
%z=(0.001148*x.^2-0.04802*x+0.5329).*y.^2+(-0.01076*x.^2+0.4567*x-5.236).*y+(0.03197*x.^2-1.43*x+18.77);%kp
 z=(0.1809*x.^2-0.5284*x+0.40033).*y.^2+(-1.651*x.^2+4.665*x-3.529).*y+(4.462*x.^2-10.82*x+8.293);%kd
surf(x,y,z),xlabel('{\it{{K}_p}}','FontName','Times New Roman','FontSize',12),ylabel('{\it{m}}(kg)','FontName','Times New Roman','FontSize',12)
shading flat
shading interp
cmap=colorbar;