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
global A
global P00
global two_ch_list_rev
global two_ch_list
global theta_sp_num
wp = 1;
wx = 1;
wy = 1;
wz = 1;
index = 1;
close all
two_ch_list = [8,9,10,4,11,6,12];
two_ch_list_rev = [1,2,3,4,5,6,7,1,2,3,5,7];
v = VideoWriter('./opt12ch.mp4','MPEG-4');
v.Quality = 100;
open(v);
loops = length( -25:0.05:0);
F(loops) = struct('cdata',[],'colormap',[]);
for tp_z = -25:0.1:--25
    %トラップ位置
    tp = [0,0,tp_z];
    %力表示するかどうか
    force_on = 0; 
    %反射有りかどうか
    reflect_on =1;
    %位相反転するかどうか
    reverse = 1;
    %アレイ位置表示するかどうk
    pos_mark = 1;
    %壁の位置
    wall_z = -25;
   
    %位相を読み込むかどうか
    load_on = 0;
        %読み込みファイルパス
        file_name = sprintf('./phase/210706/W-25%.1f.mat', tp_z);
    %位相保存するかどうか
    save_phi = 0;
    %グラフ表示するか
    graph = 1;
    %グラフ保存するか
    save_graph =0;
    save_csv = 0;
    delta_x = 1;
    delta_y = 1;
    delta_z = 1;
    x = (-20:delta_x:20);
    y = (-20:delta_y:20);
    z = (-25:delta_z:15);
    [X,Y,Z] = meshgrid(x,y,z);

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
    x_r = -sp_x;
    y_r = -sp_y;
    z_r = -sp_z;
    %縦に並ぶトランデューサの数
    theta_sp_num = 7;
    %鏡z座標
    im_z = wall_z-abs(sp_z-wall_z);
   
    if load_on == 0
        phi0 = zeros(12,1);
        options = optimoptions(@fmincon,'Display','iter','FunvalCheck','on','MaxFunctionEvaluations',3000,'PlotFcn','optimplotfval');
        % filename = 'PlotFcns0818_1'+string(tp)+'.fig';
        fun = @object_fun_12;
        [phix,fval] = fminunc(fun,phi0,options);
        % saveas(gcf,filename)
        
        elseif load_on == 1
            phix = load(file_name);
            phix = phix.phix;
        else
            phix = zeros(theta_sp_num,1);
    end

    if save_phi == 1
        save(sprintf('./phase/210726/Op12ch-25%.1f.mat', tp(3)),'phix');
    %     save('./phase/20201013/gradation.mat','phix');
    end

    % % グラデーション
    % phix = zeros(8,1);
    % for m = 1:8
    %     
    %     phix(m) = -pi/8*(8-m);
    % end
    %     
    if graph == 1
        p_n = 1;
        P = zeros(size(X));
        qc = zeros(length(sp_x),3);

        figure (2)
        hold on
        for n = 1:length(sp_x)

                xx = 0;
                if reverse  == 1
                    if sp_y(n) > 0
                        xx =1;
                        p_n = two_ch_list(p_n);
                    end
                end

                P_im = 0;

                P0 = theory_p(k,a,X,Y,Z,sp_x(n),sp_y(n),sp_z(n),0);


                if reflect_on ==1
                    P_im = theory_p(k,a,X,Y,Z,sp_x(n),sp_y(n),im_z(n),wall_z*2);
                end

                ISO = phix(p_n)+pi*xx;
                P = P+P00*A*(P0+P_im)*exp(1j*ISO);


                 if pos_mark ==1

                    tmp = mod(ISO,2*pi);

                     if tmp < pi
                         qc = [1-tmp/pi,0,0]; %赤
                         %qc = [1-p_n/12,0,0]; %赤
                     else
                         qc = [(tmp-pi)/pi,(tmp-pi)/pi,1];  %青
                         %qc = [1-p_n/12,0,0]; %赤
                         
                     end

                     q = quiver3(sp_x(n),sp_y(n),sp_z(n),-sp_x(n),-sp_y(n),-sp_z(n),0.05);
                     q.LineWidth = 5;
                     q.Color = qc;
                 end
                 p_n = two_ch_list_rev(p_n);
                if n < length(sp_x) && (sp_z(n) ~= sp_z(n+1)) && (  sp_z(n+1) > 69 || sp_z(n+1) < 65 )
                    p_n = p_n +1;
                    if p_n == theta_sp_num+1
                        p_n = 1;
                    end
                end


        end



        Power = 20*log10(abs(P));

        xslice =0;
        yslice =[];
        zslice = tp_z;

        U = poten_cal(P,delta_x,delta_y,delta_z,c0,omega);
        slice(X,Y,Z,U,xslice,yslice,zslice)
        
        %発散
        L = 6*del2(U);
        view(90,0)
        xlabel('x (mm)','FontSize',30);
        ylabel('y (mm)','FontSize',30);
        zlabel('z (mm)','FontSize',30);
        ax = gca;
        ax.FontSize = 15;
        title("Potential field")
        shading interp
        c = colorbar;
        c.Label.String = 'The Gor’kov potential';
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
        

        delete(point)
        hold off

        figure(3)
        
      
        slice(X,Y,Z,Power,xslice,yslice,zslice)
        caxis([-10,10])
        colormap jet

        % caxis([-0.04,0.04])
        view(90,0)
        xlabel('y (mm)');
        ylabel('x (mm)');
        zlabel('z (mm)');
        ax = gca;
        ax.FontSize = 15;
%         title("Potential field")
        %title("Amplitude field")
        shading interp
        %caxis([-0.01 0.01])
%         c.Label.String = 'The Gor’kov potential';
        %c.Label.String = 'Sound pressure level(dB)';
        hold on
        %目標位置
        point = plot3(tp(1),tp(2),tp(3),'o','Color','w','MarkerSize',4,'MarkerFaceColor',[0.8,0.8,0.8]);
        %colormap hot
        axis equal
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
        F(index) = getframe(gca);
        for fff = 1:10
            writeVideo(v,F(index));
        end
        delete(point)
        hold off
        if save_graph == 1
            saveas(gcf,sprintf('./210621/LSGw-25_%.1f.png', tp))
        end
        
        figure(4)
        slice(X,Y,Z,L,xslice,yslice,zslice)
                view(90,0)
        xlabel('y (mm)');
        ylabel('x (mm)');
        zlabel('z (mm)');
        xlim([min(x),max(x)])
        ylim([min(y),max(y)])
        zlim([min(z),max(z)])
        ax = gca;
        ax.FontSize = 15;
%         title("Potential field")
        %title("Amplitude field")
        shading interp
        c = colorbar;
        %caxis([-10 10])
        c.Label.String = 'the divergence of the gor’kov potential';
        %c.Label.String = 'Sound pressure level(dB)';
        hold on
        %point = plot3(tp(1),tp(2),tp(3),'o','Color','w','MarkerSize',8,'MarkerFaceColor',[0.8,0.8,0.8]);
        colormap hot
        axis equal
        if save_graph == 1
            saveas(gcf,sprintf('./210621/noref_%.1f.png', tp(3)))
        end
    %      L = 6*del2(U);
    %      L_cent = reshape(L(21,21,:),length(z),1);
    %      figure(4)
    % 
    %      plot(L_cent)

    end
    %csvファイルの保存
    %振幅位相の保存
%     filename = './csv/210625/w-25_phase_0.csv';
%     table1 = table(sp_x/delta_x,sp_y/delta_z,sp_z/delta_z,x_r/delta_x,y_r/delta_y,z_r/delta_z);
%     writetable(table1,filename);
%     
    % %音圧の保存
    if save_csv == 1
        filename = sprintf('./csv/210706/OpW-25_%.0f.csv', index);
        table2 = table(reshape(X/delta_x,length(x)*length(y)*length(z),1),reshape(Y/delta_y,length(x)*length(y)*length(z),1),reshape(Z/delta_z,length(x)*length(y)*length(z),1),reshape(L,length(x)*length(y)*length(z),1));
        writetable(table2,filename);
    end
    index = index + 1;
end
close(v);