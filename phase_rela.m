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
save_graph = 1;


delta_y = 1;
delta_z = 1;
x = 0;
y = (-20:delta_y:20);

c0 = 346*1000;
f = 40000;
omega = 2*pi*f;
k = omega/c0;
P00 = 0.17;
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



tp_z = -20:0.1:-5;
PH = zeros(theta_sp_num,length(tp_z));
A = zeros(theta_sp_num,length(tp_z));

for num = 1:length(tp_z)
    %トラップ位置
    
    file_name = sprintf('./phase/210222/LSw-25%.1f.mat', tp_z(num));

    if load_on == 1
        phix0 = load(file_name);

        PH(:,num) = phix0.phix(:);
        if load_A == 1
            A(:,num)= phix0.sin_A(:);
        else
            A(:,num)= ones(theta_sp_num,1);
        end
    else
        PH(:,num) = zeros(theta_sp_num,1);
        A(:,num)= ones(theta_sp_num,1);
    end


end
 
for num = 2:length(tp_z)
        for l = 1:theta_sp_num
            if PH(l,num)-PH(l,num-1)<-2.5
                PH(l,num:end) = PH(l,num:end)+2*pi;
            end
            
            if PH(l,num)-PH(l,num-1)>2.5
                PH(l,num:end) = PH(l,num:end)-2*pi;
            end
        end
end

fi = figure (1);
fi.InnerPosition =[0,0,800,800];
grid on
hold on 
for num = 1:theta_sp_num/2
    name = sprintf('%d ch ', num);
    
    subplot(length(theta_sp_num/2+1:theta_sp_num),1,num)
        plot(tp_z,PH(num,:),'LineWidth',2,'DisplayName',name)
 
%     title(name)
%     ylim([-pi 2*pi])
    ax = gca;
    ax.FontSize = 15;
    xlabel('position');
    ylabel('phase');
    legend show
%     
end
hold off

if save_graph == 1
    saveas(gcf,'./210427/iso1.png')
end


fi = figure (2);
fi.InnerPosition =[0,0,800,800];
grid on
hold on 
for num = theta_sp_num/2+1:theta_sp_num
    name = sprintf('%d ch ', num);
    
    subplot(length(theta_sp_num/2+1:theta_sp_num),1,num-4)
   
        plot(tp_z,PH(num,:),'LineWidth',2,'DisplayName',name)

%     title(name)
%     ylim([-pi 2*pi])
    ax = gca;
    ax.FontSize = 15;
    xlabel('position');
    ylabel('phase');
    legend show
%     
end
hold off

if save_graph == 1
    saveas(gcf,'./210427/iso2.png')
end

fi = figure (3);
fi.InnerPosition =[0,0,800,800];

hold on 
for num = 1:theta_sp_num
    name = sprintf('%d ', num);
    plot(tp_z,A(num,:)+num-1,'LineWidth',2,'DisplayName',name);
    ax = gca;
    ax.FontSize = 15;
    xlabel('position');
    ylabel('Amplitude');
    
end
hold off
legend show
if save_graph == 1
    saveas(gcf,'./210427/amplitude.png')
end

% figure (4)
% 
% hold on 
% for num = theta_sp_num/2+1:theta_sp_num
%     name = sprintf('%d CH ', num);
%     plot(tp_z,A(num,:),'LineWidth',2,'DisplayName',name);
%     ax = gca;
%     ax.FontSize = 15;
%     xlabel('position');
%     ylabel('Amplitude');
%     
% end
% hold off
% legend show
% if save_graph == 1
%     saveas(gcf,'./210416/A2.png')
% end