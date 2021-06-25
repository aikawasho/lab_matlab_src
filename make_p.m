function P = make_p(k,a,X,Y,Z,sp_x,sp_y,sp_z,f_p)%,a_func, a_func_x)
% 位相がゼロ, 振幅1の複素制限波の音圧を計算
% k 波数
% a スピーカ半径
% XYZ 音場グリッド(3次元テンソル)
% sp_x,y,x スピーカ位置
%f_p アレイ焦点
    Dim = sqrt((X-sp_x).^2+(Y-sp_y).^2+(Z-sp_z).^2);
    ip1 =(X-sp_x)*(-sp_x)+(Y-sp_y)*(-sp_y)+(Z-sp_z)*(f_p-sp_z);
    ip2 = Dim*sqrt((0-sp_x).^2+(0-sp_y).^2+(f_p-sp_z).^2);
    theta = acos(ip1./ip2);
    dire = bessel_func2(k,a,theta);
    at = load('/Users/shota/Documents/klo_lab/matlab/archive/attenuation.mat');
    attenuation = at.a_func(at.a_func_x,Dim); 
    P =2*dire.*attenuation.*exp(1j*k*Dim);
end