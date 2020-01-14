function [F] = fun3( x,y,W_j,a,b,n )%m
global  Q_1j Q_2j alfa_i alfa_o dm D omega_i 
%根据摩擦系数来求解椭圆区域磨擦力
Po=3*Q_1j*sqrt(1-(x/b(1)).^2-(y/a(1)).^2)/(2*pi*a(1)*b(1));%赫兹压力
Pi=3*Q_2j*sqrt(1-(x/b(2)).^2-(y/a(2)).^2)/(2*pi*a(2)*b(2));%赫兹压力
f=[0.519 0.519];%沟道曲率
Ro=2*f(1)*D*1e-3/(2*f(1)+1);
ro=(Ro^2-y.^2).^0.5-(Ro^2-a(1)^2)^0.5+((0.5*1e-3*D)^2-a(1)^2)^0.5;
Ri=2*f(2)*D*1e-3/(2*f(2)+1);
ri=(Ri^2-y.^2).^0.5-(Ri^2-a(2)^2)^0.5+((0.5*1e-3*D)^2-a(2)^2)^0.5;
%求速度差
wso=W_j(3)*sin(alfa_o)+W_j(2)*cos(alfa_o)+W_j(4)*sin(alfa_o);
uox=0.5*W_j(4)*1e-3*(dm)+ro*(W_j(3)*cos(alfa_o)-W_j(2)*sin(alfa_o)...
    +W_j(4)*cos(alfa_o));
uoy=0.5*D*1e-3*W_j(1);
wsi=W_j(3)*sin(alfa_i)+W_j(2)*cos(alfa_i)+(W_j(4)-omega_i)*sin(alfa_i);
uix=0.5*1e-3*(-W_j(4)+omega_i)*(dm)+ri*(W_j(3)*cos(alfa_i)-...
    W_j(2)*sin(alfa_i)-(-W_j(4)+omega_i)*cos(alfa_i));
uiy=-0.5*D*1e-3*W_j(1);
switch n
    case 1
        F=f.*Po;
    case 2
        F=-yita_o.*uoy;
    case 3
        F=0.006.*Pi;
    case 4
        F=-yita_i.*uiy;
    case 5
        F=1*yita_o.*(wso.*(x.^2+y.^2));
    otherwise
        F=1*yita_i.*(wsi.*(x.^2+y.^2));
end

