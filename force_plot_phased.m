


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

%�͕\�����邩�ǂ���
force_on = 1; 
%���˗L�肩�ǂ���
reflect_on = 1;
%�ʑ����]���邩�ǂ���
reverse = 1;
%�A���C�ʒu�\�����邩�ǂ�k
pos_mark = 1;
%�ǂ̈ʒu
wall_z = 85;
%�ʑ���ǂݍ��ނ��ǂ���
load_on = 1;
%�ʑ��ۑ����邩�ǂ���
save_phi = 0;

%�O���t�ۑ����邩
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
A = 15;
%�s�X�g���̔��a
a =4.5;
 %�X�s�[�J�ʒu
% �l�p
sp_x0 = -75:10:75;
sp_y0 = -75:10:75;
[sp_x,sp_y] = meshgrid(sp_x0,sp_y0);

sp_x = reshape(sp_x,length(sp_y0)^2,1);
sp_y = reshape(sp_y,length(sp_y0)^2,1);
sp_z = ones(length(sp_y0)^2,1);%�c�ɕ��ԃg�����f���[�T�̐�
theta_sp_num = 8;
%��z���W
im_z = wall_z-abs(sp_z-wall_z);



tp_z = 60:5:85;
st_z_upper = zeros(length(y),length(tp_z));
st_z_under = zeros(size(st_z_upper));
st_y = zeros(size(st_z_upper));
st_Z = zeros(41,length(tp_z));
z_dis = (-20:20);
for num = 1:length(tp_z)
    %�g���b�v�ʒu
    tp = [0,0,tp_z(num)];

    %�ǂݍ��݃t�@�C���p�X
    file_name = sprintf('./phase/20201211/W-25phased_%.1f.mat', tp(3));
    %�U���ǂݍ��ނ��ǂ���
    load_A = 0;
    z = (tp_z(num)-20:delta_z:tp_z(num)+20);
    [Y,Z] = meshgrid(y,z);

    if load_on == 1
        phix0 = load(file_name);
        phix = phix0.phix;
        if load_A == 1
            sinA = phix0.sin_A;
        else
            sinA = ones(length(sp_y0)^2,1);
        end
    else
        phix = zeros(length(sp_y0)^2,1);
    end



    ph_n = 1;
    P = zeros(size(Y));



        for n = 1:length(sp_x)


                P_im = 0;

                P0 = theory_p_2dflat(k,a,0,Y,Z,sp_x(n),sp_y(n),sp_z(n),1);


                if reflect_on ==1
                    P_im = theory_p_2dflat(k,a,0,Y,Z,sp_x(n),sp_y(n),im_z(n),-1);
                end

                ISO = phix(n);
                P = P+sinA(n)*P00*A*(P0+P_im)*exp(1j*ISO);


        end



        Power = 20*log10(abs(P));


        U = poten_cal_2d(P,delta_y,delta_z,c0,omega);
        %��ɂ�����d��24.5 mN
        %���a1.5mm��EPS�ɂ�����d�� 0.004 mN

        [F_y, F_z] = gradient(U,delta_y,delta_z);

        st_z_under(:,num) =reshape(-F_z(17,:),size(y));
        st_z_upper(:,num) =reshape(-F_z(23,:),size(y));
        st_y(:,num) = reshape(-F_y(6,:),size(y));
        st_Z(:,num) = reshape(-F_z(:,21),size(z));
        

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
    plot(y,st_z_under(:,num),'LineWidth',2,'DisplayName',name);
    ax = gca;
    ax.FontSize = 15;
    xlabel('y-axis displacement of desired trap position (mm)');
    ylabel('Acoustic Radiation Force (mN)');
    title('z-axis force')
end
hold off
legend show
if save_graph == 1
    saveas(gcf,'./20201211/phased_y-z_no.png')
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
    saveas(gcf,'./20201211/phased_y_y_no.png')
end

figure(3)
hold on
for num = 1:length(tp_z)    
    name = sprintf('%.1f mm ', tp_z(num));
    plot(z_dis,st_Z(:,num),'LineWidth',2,'DisplayName',name);
    ax = gca;
    ax.FontSize = 15;
    xlabel('z-axis displacement of desired trap position (mm)');
    ylabel('Acoustic Radiation Force (mN)');
    title('z-axis force')
end
hold off
legend show
if save_graph == 1
    saveas(gcf,'./20201211/phased_z_z_no.png')
end
