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
for tp_z = -25:0.5:-0
    %�g���b�v�ʒu
    tp = [0,0,tp_z];
    %�͕\�����邩�ǂ���
    force_on = 0; 
    %���˗L�肩�ǂ���
    reflect_on = 1;
    %�ʑ����]���邩�ǂ���
    reverse = 1;
    %�A���C�ʒu�\�����邩�ǂ�k
    pos_mark = 0;
    %�ǂ̈ʒu
    wall_z = -25;
    %�ʑ���ǂݍ��ނ��ǂ���
    load_on = 0;
    %�ʑ��ۑ����邩�ǂ���
    save_phi = 1;
    %�O���t�\�����邩
    graph = 0;
    %�O���t�ۑ����邩
    save_graph = 0;
        %�ǂݍ��݃t�@�C���p�X
    file_name = './phase/20201112/LSw-25-20.0.mat';
    delta_x = 0.5;
    delta_y = 0.5;
    delta_z = 0.5;
    x = (-20:delta_x:20);
    y = (-20:delta_y:20);
    z = (wall_z:delta_z:wall_z+40);
    [X,Y,Z] = meshgrid(x,y,z);

    c0 = 346*1000;
    f = 40000;
    omega = 2*pi*f;
    k = omega/c0;
    P00 = 1;

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
    phi0 = ones(theta_sp_num,2);

    options = optimoptions(@fmincon,'Display','iter','FunvalCheck','on','MaxFunctionEvaluations',3000,'PlotFcn','optimplotfval','Algorithm','sqp');
%     options = optimoptions(@fmincon,'Display','iter','FunvalCheck','on','MaxFunctionEvaluations',3000,'PlotFcn','optimplotfval');

    % filename = 'PlotFcns0818_1'+string(tp)+'.fig';
    fun = @object_fun_3;
    A2 = [];
    b = [];
    Aeq = [];
    beq = [];
    lb = -Inf;
    ub = Inf;
    nonlcon = @func3;
    [phix,fval] = fmincon(fun,phi0,A2,b,Aeq,beq,lb,ub,nonlcon,options);
%     [phix,fval] = fminunc(fun,phi0,options);
    % saveas(gcf,filename)


    if save_phi == 1
        save(sprintf('./phase/210520/nonre_w-25%.1f.mat', tp(3)),'phix');
    %     save('./phase/20201013/gradation.mat','phix');
    end

    % % �O���f�[�V����
    % phix = zeros(8,1);
    % for m = 1:8
    %     
    %     phix(m) = -pi/8*(8-m);
    % end
    %     
    if graph == 1
        ph_n = 1;
        P = zeros(size(X));
        qc = zeros(length(sp_x),3);

        figure (2)
        hold on
        for n = 1:length(sp_x)

                xx = 0;
                if reverse  == 1
                    if sp_y(n) < 0 
                        xx =1;
                    end
                end

                P_im = 0;

                P0 = theory_p(k,a,X,Y,Z,sp_x(n),sp_y(n),sp_z(n),0);


                if reflect_on ==1
                    P_im = theory_p(k,a,X,Y,Z,sp_x(n),sp_y(n),im_z(n),wall_z*2);
                end

                ISO = phix(ph_n,1)+pi*xx;
                P = P+P00*abs(1.5*sin(phix(ph_n,2)))*(P0+P_im)*exp(1j*ISO);


                 if pos_mark ==1

                    tmp = mod(ISO,2*pi);

                     if tmp < pi
                         qc = [1-tmp/pi,0,0]; %��
                     else
                         qc = [(tmp-pi)/pi,(tmp-pi)/pi,1];  %��
                     end

                     q = quiver3(sp_x(n),sp_y(n),sp_z(n),-sp_x(n),-sp_y(n),-sp_z(n),0.05);
                     q.LineWidth = 5;
                     q.Color = qc;
                 end

                if n < length(sp_x) && (sp_z(n) ~= sp_z(n+1))
                    ph_n = ph_n +1;
                    if ph_n == theta_sp_num+1
                        ph_n = 1;
                    end
                end


        end



        Power = 20*log10(abs(P));

        xslice =0;
        yslice =[];
        zslice = wall_z;

        U = poten_cal(P,delta_x,delta_y,delta_z,c0,omega);
        slice(X,Y,Z,U,xslice,yslice,zslice)

        view(90,0)
        xlabel('x (mm)','FontSize',30);
        ylabel('y (mm)','FontSize',30);
        zlabel('z (mm)','FontSize',30);
        ax = gca;
        ax.FontSize = 15;
        title("Potential field")
        shading interp
        c = colorbar;
        c.Label.String = 'The Gor�fkov potential';
        caxis([0 0.25])
        axis equal
        point = plot3(tp(1),tp(2),tp(3),'o','Color','w','MarkerSize',8,'MarkerFaceColor',[0.8,0.8,0.8]);
        colormap jet
        if force_on == 1
            [F_x0, F_y0, F_z0] =  gradient(U,delta_x,delta_y,delta_z);
            span = 2;
            F_x0 = -F_x0(1:span:end,1:span:end,1:span:end); F_y0 = -F_y0(1:span:end,1:span:end,1:span:end); F_z0 = -F_z0(1:span:end,1:span:end,1:span:end);

            F_x0 = F_x0/max(max(max(F_x0)));
            F_x0(abs(F_x0)<0.5) = 0;

            F_y0 = F_y0/max(max(max(F_y0)));
            F_y0(abs(F_y0)<0.5) = 0;


            F_z0 = F_z0/max(max(max(F_z0)));
            F_z0(abs(F_z0)<0.5) = 0;

            quiver3(X(1:span:end,1:span:end,1:span:end),Y(1:span:end,1:span:end,1:span:end),Z(1:span:end,1:span:end,1:span:end),real(F_x0),real(F_y0),real(F_z0),1,"red")
        end
        

%         delete(point)
        hold off

        figure(3)
        slice(X,Y,Z,Power,xslice,yslice,zslice)
        view(90,0)
        ax = gca;
        ax.FontSize = 15;
        xlabel('x2');
        ylabel('x1');
        zlabel('z (mm)');

        c = colorbar;
        c.Label.String = 'Sound pressure level(dB)';
        caxis([-10,10])
        shading interp
        axis equal
        if save_graph == 1
            saveas(gcf,sprintf('./210121/conw-25_%.1f.png', tp(3)))
        end
    %      L = 6*del2(U);
    %      L_cent = reshape(L(21,21,:),length(z),1);
    %      figure(4)
    % 
    %      plot(L_cent)
    %csv�t�@�C���̕ۑ�
%         %�U���ʑ��̕ۑ�
%         filename = 'opw-25_tp-10.csv';
%         table1 = table(sp_x,sp_y,sp_z,x_r,y_r,z_r,ang,pow);
%         writetable(table1,filename);

        %�����̕ۑ�
        %filename = 'opdataw-25_tp-10.csv';
        %table2 = table(reshape(X,length(x)*length(y)*length(z),1),reshape(Y,length(x)*length(y)*length(z),1),reshape(Z,length(x)*length(y)*length(z),1),reshape(Power,length(x)*length(y)*length(z),1));
        %writetable(table2,filename);
    end
end