%%model-1 �⾲��ѧ����
clear 
close all
clc
global alfa A Zn F K_oj K_ij f D omega_i m I dm Q_1j...
    Q_2j alfa_i alfa_o Ic omega_Rj omega_cj
%��������
f=[0.519 0.519];%��������
dm=56.6;%��λ��mm
D=13.5;%��λ��mm
alfa=25;%��ʼ�Ӵ���/���
Zn=11;%���Ӹ���
den=7800;%�����ܶ� kg/m^3
m=pi*den*(D*1e-3)^3/6;%�������� kg
I=m*(D*1e-3)^2/10;%�����ת������ kg*m^2
Ic=I+0.25*m*(dm*1e-3)^2;%������������ߵ�ת������
K_ij =1.0856e+06;%��λ N*mm^1.5
K_oj =1.1356e+06;%��λ N*mm^1.5
A=(f(1)+f(2)-1)*D;%�������ľ��� mm
%ϵͳ����
dF=1;W_x=[];W_y=[];W_z=[];W_c=[];rate=[];
omega_i=4000*pi/30;%��Ȧת�� rad/s
T=2*pi/omega_i;%��Ȧת������ s
omega_cj=0.5*omega_i*(1-D*cosd(alfa)/dm);
alfa_j=atan(sind(alfa)/(cosd(alfa)+D/dm));
omega_Rj=0.5*omega_i*dm/D*(1-(D*cosd(alfa)/dm)^2)*sin(alfa_j);
W1=2*[200 800 -100 82*0.9 2*pi/Zn 82*0.9 0];%���ٶȳ�ֵ������7*1
W=W1;
output=[];le=0;dnum=1;WC=[];
for F=50;
for i=1:50
deta_a=0.0107;
X=[0.1129 0.2328 0.023 0.023 deta_a];%��ֵ������5*1
options =optimoptions('fsolve','Algorithm','Levenberg-Marquardt');
[Y,fval,exitflag]=fsolve(@fun,X,options);
if exitflag==1
    disp('�������');
else
    disp('�����ֵ��Ч');
end
Q_1j=K_oj*Y(3)^1.5;%��λ/N
Q_2j=K_ij*Y(4)^1.5;%��λ/N
A_2j=A*cosd(alfa);
cos_1j=Y(2)/((f(1)-0.5)*D+Y(3));
cos_2j=(A_2j-Y(2))/((f(2)-0.5)*D+Y(4));
alfa_i=acos(cos_2j);
alfa_o=acos(cos_1j);
tspan = 0:1e-5:0.3*T;
[t,y] = ode15s(@fun1,tspan,W,options);
num=length(t);
W=y(num,:);%���ٶȳ�ֵ������7*1
omega_cj=y(num,4);
omega_Rj=y(num,1);
output(1+le:le+num,:)=y(1:num,:);
le=le+num;
end
temp=length(output(:,4));
WC(dnum)=output(temp,4);
dnum=dnum+1;
end
TT=WC'/418.88
% num1=length(output(1,:));
% Y_fin=output(num1,:)
% suboutput=fun4(Y_fin)