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

%力表示するか
force_on = 1; 
%反射有りか
reflect_on = 1;
%位相反転するか
reverse = 1;
%アレイ位置表示するか
pos_mark = 1;
%壁の位置
wall_z = -25;
%位相を読み込むか
load_on = 1;


%グラフ保存するか
save_graph = 0;


delta_y = 1;
delta_z = 1;
x = 0;
y = (-20:delta_y:20);

c0 = 346*1000;
f = 40000;
omega = 2*pi*f;
k = omega/c0;
P00 = 0.17;
A = 15;
%ピストンの半径
a =4.5;
%アレイ半径
zahyo = load('./zahyo/20200720_180.mat');

%スピーカ位置
sp_x = zahyo.X;
sp_y = zahyo.Y;
sp_z = zahyo.Z;
%縦に並ぶトランデューサの数
theta_sp_num = 8;
%鏡z座標
im_z = wall_z-abs(sp_z-wall_z);



tp_z = -26:5:-26;
st_z_upper = zeros(length(y),length(tp_z));
st_z_under = zeros(size(st_z_upper));
st_y = zeros(size(st_z_upper));
st_Z = zeros(41,length(tp_z));
st_Zy = zeros(41,length(tp_z));
z_dis = (-20:20);
for num = 1:length(tp_z)
    %トラップ位置
    tp = [0,0,tp_z(num)];

    %読み込みファイルパス
    file_name = sprintf('./phase/210602/LSGw-30%.1f.mat', tp(3));
    %振幅読み込むかどうか
    load_A = 1;
    z = (tp_z(num)-20:delta_z:tp_z(num)+20);
    [Y,Z] = meshgrid(y,z);

    if load_on == 1
        phix0 = load(file_name);
        phix = phix0.phix(:);
        if load_A == 1
            sinA = phix0.sin_A(:);
        else
            sinA = ones(theta_sp_num,1);
        end
    else
        phix = zeros(theta_sp_num,1);
        sinA = ones(theta_sp_num,1);
    end



    ph_n = 1;
    P = zeros(size(Y));



        for n = 1:length(sp_x)

                xx = 0;
                if reverse  == 1
                    if sp_y(n) < 0 
                        xx =1;
                    end
                end

                P_im = 0;

                %P0 = theory_p_2d(k,a,0,Y,Z,sp_x(n),sp_y(n),sp_z(n),0);
                P0 = make_p(k,a,0,Y,Z,sp_x(n),sp_y(n),sp_z(n),0);


                if reflect_on ==1
                    %P_im = theory_p_2d(k,a,0,Y,Z,sp_x(n),sp_y(n),im_z(n),wall_z*2);
                    P0 = make_p(k,a,0,Y,Z,sp_x(n),sp_y(n),sp_z(n),0);
                end

                ISO = phix(ph_n)+pi*xx;
                P = P+sinA(ph_n)*P00*A*(P0+P_im)*exp(1j*ISO);

                if n < length(sp_x) && (sp_z(n) ~= sp_z(n+1))
                    ph_n = ph_n +1;
                    if ph_n == theta_sp_num+1
                        ph_n = 1;
                    end
                end


        end



        Power = 20*log10(abs(P));


        U = poten_cal_2d(P,delta_y,delta_z,c0,omega);
        %蚊にかかる重力24.5 mN
        %半径1.5mmのEPSにかかる重力 0.004 mN

        [F_y, F_z] = gradient(U,delta_y,delta_z);

%         st_z_under(:,num) =reshape(-F_z(6,:),size(y));
        st_z_upper(:,num) =reshape(-F_z(23,:),size(y));
        st_y(:,num) = reshape(-F_y(6,:),size(y));
        st_Z(:,num) = reshape(-F_z(:,21),size(z));
        st_Zy(:,num) = reshape(-F_y(:,20),size(y));


end

figure (1)

% subplot(1,2,1)
% hold on
% for num = 1:length(tp_z)
%     name = sprintf('%.1f mm ', tp_z(num));
%     plot(y,st_z_under(:,num),'LineWidth',2,'DisplayName',name);
%     ax = gca;
%     ax.FontSize = 15;
%     xlabel('y (mm)');
%     ylabel('Acoustic Radiation Force (mN)');
%     title('under')
% %     ylim([-0.2,0.05])
% end
% hold off
% legend show
% 
% subplot(1,2,2)
hold on 
for num = 1:length(tp_z)
    name = sprintf('%.1f mm ', tp_z(num));
    plot(y,st_z_upper(:,num),'LineWidth',2,'DisplayName',name);
    ax = gca;
    ax.FontSize = 15;
    xlabel('y-axis displacement of desired trap position (mm)');
    ylabel('Acoustic Radiation Force (mN)');
    title('z-axis force')
end
hold off
legend show
if save_graph == 1
    saveas(gcf,'./210416/2010y-z_2.png')
end

figure(2)
hold on
for num = 1:length(tp_z)    
    name = sprintf('%.1f mm ', tp_z(num));
    plot(y,st_y(:,num),'LineWidth',2,'DisplayName',name);
    ax = gca;
    ax.FontSize = 15;
    xlabel('y-axis displacement of desired trap position (mm)');
    ylabel('Acoustic Radiation Force (mN)');
    title('y-axis force')
end
hold off
legend show
if save_graph == 1
    saveas(gcf,'./210416/2010y_y_2.png')
end

figure(3)
hold on
bi = 0;
for num = 1:length(tp_z)    
    name = sprintf('%.1f mm ', tp_z(num));
    plot(z_dis,st_Z(:,num),'LineWidth',2,'DisplayName',name);
    ax = gca;
    ax.FontSize = 15;
    %ylim([-0.3 0.3])
    xlabel('z-axis displacement of desired trap position (mm)');
    ylabel('Acoustic Radiation Force (mN)');
    title('z-axis force')
    bi = bi+1;
end
hold off
legend show
if save_graph == 1
    saveas(gcf,'./210416/2010z_z_2.png')
end

figure(4)
hold on
bi = 0;
for num = 1:length(tp_z)    
    name = sprintf('%.1f mm ', tp_z(num));
    plot(z_dis,st_Zy(:,num)+0.15*(num-1),'LineWidth',2,'DisplayName',name);
    ax = gca;
    ax.FontSize = 15;
    %ylim([0 0.3])
    xlabel('z-xis displacement of desired trap position (mm)');
    ylabel('Acoustic Radiation Force (mN)');
    title('z-axis force')
    bi = bi+1;
end
hold off
legend show
if save_graph == 1
    saveas(gcf,'./210416/2010z_y_2.png')
end
