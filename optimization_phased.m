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
wx = 1000;
wy = 1;
wz = 1;

close all

for tp_z = 60:1:85

    %トラップ位置
    tp = [0,0,tp_z];
    %力表示するかどうか
    force_on = 0; 
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
    %グラフ表示するか
    graph = 1;
    %グラフ保存するか
    save_graph = 1;
    %読み込みファイルパス
    file_name = './phase/20201211/W-25phased_80.0.mat';
    delta_x = 1;
    delta_y = 1;
    delta_z = 1;
    x = (-10:delta_x:10);
    y = (-20:delta_y:20);
    z = (1:delta_z:wall_z);
    [X,Y,Z] = meshgrid(x,y,z);

    c0 = 346*1000;
    f = 40000;
    omega = 2*pi*f;
    k = omega/c0;
    P00 = 0.17;
    A = 15;
    %ピストンの半径
    a =4.5;
    % 四角
    sp_x0 = -75:10:75;
    sp_y0 = -75:10:75;
    [sp_x,sp_y] = meshgrid(sp_x0,sp_y0);

    sp_x = reshape(sp_x,length(sp_y0)^2,1);
    sp_y = reshape(sp_y,length(sp_y0)^2,1);
    sp_z = ones(length(sp_y0)^2,1);

    %鏡z座標
    im_z = wall_z-abs(sp_z-wall_z);

    if load_on == 0
        phi0 = zeros(size(sp_z));
        options = optimoptions(@fmincon,'Display','iter','FunvalCheck','on','MaxFunctionEvaluations',4000,'PlotFcn','optimplotfval');
        % filename = 'PlotFcns0818_1'+string(tp)+'.fig';
        fun = @object_fun2;
        [phix,fval] = fminunc(fun,phi0,options);
        % saveas(gcf,filename)
    elseif load_on == 1
        phix = load(file_name);
        phix = phix.phix;
    else
        phix = zeros(theta_sp_num,1);
    end

    if save_phi == 1
        save(sprintf('./phase/20201211/W-25phased_%.1f.mat', tp(3)),'phix');
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
        ph_n = 1;
        P = zeros(size(X));
        qc = zeros(length(sp_x),3);

        figure (2)
        hold on
        for n = 1:length(sp_x)


                P_im = 0;

                P0 = theory_p_flat(k,a,X,Y,Z,sp_x(n),sp_y(n),sp_z(n),1);


                if reflect_on ==1
                    P_im = theory_p_flat(k,a,X,Y,Z,sp_x(n),sp_y(n),im_z(n),-1);
                end

                ISO = phix(n);
                P = P+P00*A*(P0+P_im)*exp(1j*ISO);


                 if pos_mark ==1

                    tmp = mod(ISO,2*pi);

                     if tmp < pi
                         qc = [1-tmp/pi,0,0]; %赤
                     else
                         qc = [(tmp-pi)/pi,(tmp-pi)/pi,1];  %青
                     end

                     q = quiver3(sp_x(n),sp_y(n),sp_z(n),0,0,1,1);
                     q.LineWidth = 1;
                     q.Color = qc;
                 end


        end



        Power = 20*log10(abs(P));

        xslice =0;
        yslice =[];
        zslice = tp_z;

        U = poten_cal(P,delta_x,delta_y,delta_z,c0,omega);
        slice(X,Y,Z,U,xslice,yslice,zslice)

        view(90,0)
        xlabel('x (mm)');
        ylabel('y (mm)');
        zlabel('z (mm)');
        ax = gca;
        ax.FontSize = 15;
        title("Potential field")
        shading interp
        c = colorbar;
        c.Label.String = 'The Gor’kov potential';
        caxis([0 0.40])
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
        zlim([min(z),max(z)])

        hold off

        figure(3)
        hold on
        slice(X,Y,Z,Power,xslice,yslice,zslice)
        view(90,0)
        ax = gca;
        ax.FontSize = 15;
        xlabel('x (mm)');
        ylabel('y (mm)');
        zlabel('z (mm)');
        title("Amplitude field")
        c = colorbar;
        c.Label.String = 'Sound pressure level(dB)';
        caxis([-10,10])
        shading interp
        axis equal
        point = plot3(tp(1),tp(2),tp(3),'o','Color','w','MarkerSize',8,'MarkerFaceColor',[0.8,0.8,0.8]);
        if save_graph == 1
            saveas(gcf,sprintf('./20201211/W85phased_%.1f.png', tp(3)))
        end
%         delete(point)
        zlim([min(z),max(z)])
        colormap hot
        hold off
    %      L = 6*del2(U);
    %      L_cent = reshape(L(21,21,:),length(z),1);
    %      figure(4)
    % 
    %      plot(L_cent)
    end
end