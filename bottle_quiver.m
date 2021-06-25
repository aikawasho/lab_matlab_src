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

for tp_z = -16:4:-16

    %�g���b�v�ʒu
    tp = [0,0,tp_z];
    %�͕\�����邩�ǂ���
    force_on = 1; 
    %���˗L�肩�ǂ���
    reflect_on = 1;
    %�ʑ����]���邩�ǂ���
    reverse = 1;
    %�A���C�ʒu�\�����邩�ǂ�k
    pos_mark = 0;
    %�ǂ̈ʒu
    wall_z = -25;
    %�ʑ���ǂݍ��ނ��ǂ���
    load_on = 1;
    %�ʑ��ۑ����邩�ǂ���
    save_phi = 0;
     %�U���ǂݍ��ނ��ǂ���
    load_A = 0;
    %�O���t�ۑ����邩
    save_graph = 0;
    %�ǂݍ��݃t�@�C���p�X
    file_name = sprintf('./phase/210501/BoW-25%.1f.mat', tp(3));

    
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
    %�s�X�g���̔��a
    a =4.5;
    %�A���C���a
    zahyo = load('./zahyo/20200720_180.mat');

    %�X�s�[�J�ʒu
    sp_x = zahyo.X;
    sp_y = zahyo.Y;
    sp_z = zahyo.Z;
    %�c�ɕ��ԃg�����f���[�T�̐�
    theta_sp_num = 8;
    %��z���W
    im_z = wall_z-abs(sp_z-wall_z);

    if load_on == 1
        phix0 = load(file_name);
        phix = phix0.phix(:,1);
%           phix = phix0.phix;
          
        if load_A == 1
            sinA = phix0.phix(:,2);
%             sinA = phix0.sin_A;
        else
            sinA = ones(8,1);
        end
    else
        phix = zeros(theta_sp_num,1);
        sinA = ones(8,1);
    end



    ph_n = 1;
    P_yz = zeros(size(Y1));
    P_xy = zeros(size(X2));
    qc = zeros(length(sp_x),3);


    for n = 1:length(sp_x)

            xx = 0;
            if reverse  == 1
                if sp_z(n) > 40 
                    xx =1;
                end
            end

            P_im_1 = 0;
            P_im_2 = 0;

            P0_1 = theory_p_2d(k,a,0,Y1,Z1,sp_x(n),sp_y(n),sp_z(n),0);
            P0_2 = theory_p_2d(k,a,X2,Y2,tp_z,sp_x(n),sp_y(n),sp_z(n),0);


            if reflect_on ==1
                P_im_1 = theory_p_2d(k,a,0,Y1,Z1,sp_x(n),sp_y(n),im_z(n),wall_z*2);
                P_im_2 = theory_p_2d(k,a,X2,Y2,tp_z,sp_x(n),sp_y(n),im_z(n),wall_z*2);
            end

            ISO = phix(ph_n)+pi*xx;
%             P_yz = P_yz+sinA(ph_n)*P00*A*(P0_1+P_im_1)*exp(1j*ISO);
%             P_xy = P_xy+sinA(ph_n)*P00*A*(P0_2+P_im_2)*exp(1j*ISO);

            P_yz = P_yz+abs(1.5*sin(sinA(ph_n)))*P00*A*(P0_1+P_im_1)*exp(1j*ISO);
            P_xy = P_xy+abs(1.5*sin(sinA(ph_n)))*P00*A*(P0_2+P_im_2)*exp(1j*ISO);
            if n < length(sp_x) && (sp_z(n) ~= sp_z(n+1))
                ph_n = ph_n +1;
                if ph_n == theta_sp_num+1
                    ph_n = 1;
                end
            end


    end


    U_yz = poten_cal_2d(P_yz,delta_y,delta_z,c0,omega);
    U_xy = poten_cal_2d(P_xy,delta_x,delta_z,c0,omega);
    %��ɂ�����d��24.5 mN
    %���a1.5mm��EPS�ɂ�����d�� 0.004 mN

    [F_y1, F_z1] = gradient(U_yz,delta_y,delta_z);
    [F_x2, F_z2] = gradient(U_xy,delta_y,delta_z);

%     abs_f_yz = sqrt(F_y1.^2+F_z1.^2);
%     abs_f_xy = sqrt(F_x2.^2+F_z2.^2);
%     if save_graph == 1
%         save(sprintf('./210121/abs_f_yz_%.1f.mat', tp(3)),'abs_f_yz')
%         save(sprintf('./210121/abs_f_xy_%.1f.mat', tp(3)),'abs_f_xy')
%     end

    figure(1)
        surf(Y2,Z_dis,U_yz,'FaceAlpha',0.5)
        shading interp
        hold on 
        quiver(Y00,Z_dis,-F_y1,-F_z1,1.0,'Color',[0,0,0])
        hold on 
        %point1 = scatter(tp(2),tp(3),250,'MarkerEdgeColor',[0.8,0.8,0.8],'MarkerFaceColor',[0.8,0.8,0.8],'MarkerFaceAlpha',0.5);
        %quiver(Y,Z,F_y,F_z)
        xlim([min(y),max(y)])
        ylim([min(z_dis),max(z_dis)])
        ax = gca;
        ax.FontSize = 15;
        xlabel('x-axis displacement of desired trap position (mm)');
        ylabel('z-axis displacement of desired trap position (mm)');
        title('Acoustic Radiation Force x-z plane')
        c = colorbar;
        c.Label.String = 'The Gor�fkov potential';
        caxis([0 1])
        view(0,90)
        hold off

    if save_graph == 1
        saveas(gcf,sprintf('./210416/nut_xz%.1f.png', tp(3)))
    end

    figure(2)
        surf(Y1,X2,U_xy,'FaceAlpha',0.5)
        shading interp
        hold on 
        quiver(Y00,Z_dis,-F_x2,-F_z2,1.0,'Color',[0,0,0])
        
        %point2 = scatter(0,0,250,'MarkerEdgeColor',[0.8,0.8,0.8],'MarkerFaceColor',[0.8,0.8,0.8],'MarkerFaceAlpha',0.5);
        %quiver(Y,Z,F_y,F_z)
        xlim([min(y),max(y)])
        ylim([min(x),max(x)])
        ax = gca;
        ax.FontSize = 15;
        xlabel('x-axis displacement of desired trap position (mm)');
        ylabel('y-axis displacement of desired trap position (mm)');
        title('Acoustic Radiation Force x-y plane')
        c = colorbar;
        c.Label.String = 'The Gor�fkov potential';
        caxis([0 1])
        view(0,90)
        hold off
    if save_graph == 1
        saveas(gcf,sprintf('./210416/nut_xy%.1f.png', tp(3)))
    end
end