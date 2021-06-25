function P = theory_p_2dflat(k,a,X,Y,Z,sp_x,sp_y,sp_z,muki)
% 位相がゼロ, 振幅1の複素制限波の音圧を計算
% k 波数
% a スピーカ半径
% XYZ 音場グリッド(3次元テンソル)
% sp_x,y,x スピーカ位置

    Dim = sqrt((X-sp_x).^2+(Y-sp_y).^2+(Z-sp_z).^2);
    %内積
    ip1 =(Z-sp_z)*(sp_z*muki);
    %ベクトル大きさ
    ip2 = Dim.*sqrt((sp_z).^2);
    theta = acos(ip1./ip2);
    
    dire = bessel_func2(k,a,theta);
    
    P =dire./Dim.*exp(1j*k*Dim);

end