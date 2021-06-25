

global wall_z
global c0
global f
global omega
global k
global a
global im_z
global sp_x
global sp_y
global sp_z
global tp
global wx
global wy
global wz
global wp
global reflect_on
global reverse


global theta_sp_num
wp = 1;
wx = 1;
wy = 1;
wz = 1;

close all
z_dis = (-10:10);

for tp_z = 80:10:80

    %トラップ位置
    tp = [0,0,tp_z];
    %力表示するかどうか
    force_on = 1; 
    %反射有りかどうか
    reflect_on = 1;
    %位相反転するかどうか
    reverse = 1;
    %アレイ位置表示するかどうk
    pos_mark = 0;
    %壁の位置
    wall_z = 85;
    %位相を読み込むかどうか
    load_on = 1;
    %位相保存するかどうか
    save_phi = 0;
    
    %グラフ保存するか
    save_graph = 0;
    %読み込みファイルパス
    file_name = sprintf('./phase/210120/W-25phased_%.1f.mat', tp(3));
    %振幅読み込むかどうか
    load_A = 0;
    
    delta_x = 1;
    delta_y = 1;
    delta_z = 1;
    x = (-10:delta_x:10);
    y = (-10:delta_y:10);
    z = (tp_z-10:delta_z:tp_z+10);
    [Y1,Z1] = meshgrid(y,z);
    [Y2,X2] = meshgrid(y,x);
    [Y00,Z_dis] = meshgrid(y,z_dis);

    c0 = 346*1000;
    f = 40000;
    omega = 2*pi*f;
    k = omega/c0;
    P00 = 0.17;
    A = 15;
    %ピストンの半径
    a =4.5;

    %スピーカ位置
    % 四角
    sp_x0 = -75:10:75;
    sp_y0 = -75:10:75;
    [sp_x,sp_y] = meshgrid(sp_x0,sp_y0);

    sp_x = reshape(sp_x,length(sp_y0)^2,1);
    sp_y = reshape(sp_y,length(sp_y0)^2,1);
    sp_z = ones(length(sp_y0)^2,1);

    %縦に並ぶトランデューサの数
    theta_sp_num = 8;
    %鏡z座標
    im_z = wall_z-abs(sp_z-wall_z);

    if load_on == 1
        phix0 = load(file_name);
        phix = phix0.phix;
        if load_A == 1
            sinA = phix0.sin_A;
        else
            sinA = ones(length(sp_z),1);
        end
    else
        phix = zeros(theta_sp_num,1);
    end



    ph_n = 1;
    P_yz = zeros(size(Y1));
    P_xy = zeros(size(X2));
    qc = zeros(length(sp_x),3);


    for n = 1:length(sp_x)

            P_im_1 = 0;
            P_im_2 = 0;

            P0_1 = theory_p_2dflat(k,a,0,Y1,Z1,sp_x(n),sp_y(n),sp_z(n),1);
            P0_2 = theory_p_2dflat(k,a,X2,Y2,tp_z,sp_x(n),sp_y(n),sp_z(n),1);


            if reflect_on ==1
                P_im_1 = theory_p_2dflat(k,a,0,Y1,Z1,sp_x(n),sp_y(n),im_z(n),-1);
                P_im_2 = theory_p_2dflat(k,a,X2,Y2,tp_z,sp_x(n),sp_y(n),im_z(n),-1);
            end

            ISO = phix(n);
            P_yz = P_yz+sinA(n)*P00*A*(P0_1+P_im_1)*exp(1j*ISO);
            P_xy = P_xy+sinA(n)*P00*A*(P0_2+P_im_2)*exp(1j*ISO);


    end


    U_yz = poten_cal_2d(P_yz,delta_y,delta_z,c0,omega);
    U_xy = poten_cal_2d(P_xy,delta_x,delta_z,c0,omega);
    %蚊にかかる重力24.5 mN
    %半径1.5mmのEPSにかかる重力 0.004 mN

    [F_y1, F_z1] = gradient(U_yz,delta_y,delta_z);
    [F_x2, F_z2] = gradient(U_xy,delta_y,delta_z);
    
    F_y1(isnan(F_y1)) = 0;
    F_z1(isnan(F_z1)) = 0;
    F_y1(isinf(F_y1)) = 0;
    F_z1(isinf(F_z1)) = 0;
    abs_f = sqrt(F_y1.^2+F_z1.^2);


    figure(1)
        quiver(Y00,Z_dis,-F_y1,-F_z1,0.8,'Color',[0,0,0])
        hold on 
        %point1 = scatter(tp(2),tp(3),250,'MarkerEdgeColor',[0.8,0.8,0.8],'MarkerFaceColor',[0.8,0.8,0.8],'MarkerFaceAlpha',0.5);
        %quiver(Y,Z,F_y,F_z)
        xlim([min(y),max(y)])
        ylim([min(z_dis),max(z_dis)])
        ax = gca;
        ax.FontSize = 15;
        xlabel('y-axis displacement of desired trap position (mm)');
        ylabel('z-axis displacement of desired trap position (mm)');
        title('Acoustic Radiation Force y-z plane')
        hold off

    if save_graph == 1
        saveas(gcf,sprintf('./20201209/y_z_LSw-25_%.1f.png', tp(3)))
    end

    figure(2)
        quiver(Y00,Z_dis,-F_x2,-F_z2,0.8,'Color',[0,0,0])
        hold on 
        %point2 = scatter(0,0,250,'MarkerEdgeColor',[0.8,0.8,0.8],'MarkerFaceColor',[0.8,0.8,0.8],'MarkerFaceAlpha',0.5);
        %quiver(Y,Z,F_y,F_z)
        xlim([min(y),max(y)])
        ylim([min(x),max(x)])
        ax = gca;
        ax.FontSize = 15;
        xlabel('x-axis displacement of desired trap position (mm)');
        ylabel('y-axis displacement of desired trap position (mm)');
        title('Acoustic Radiation Force x-y plane')
        hold off
    if save_graph == 1
        saveas(gcf,sprintf('./20201209/x-y_LSw-25_%.1f.png', tp(3)))
    end
end