close all
delta_x = 1;
delta_y = 1;
delta_z = 1;
for tp_z = 0:0.1:10
    %トラップ位置Z
    tp = tp_z;
    %力表示するかどうか
    force_on = 0;
    %アレイ位置表示
    pos_mark = 0;
    %反射有りかどうか
    reflect_on = 1;
    %位相反転するかどうか
    reverse = 1;
    %壁の位置
    wall_z = -25;
    %振幅, 位相保存するかどうか
    save_w = 1;
    %グラフ保存するかどうか
    save_graph = 0;
    %グラフ表示するかどうか
    write_graph = 0;

    x = (-20:delta_x:20);
    y = (-20:delta_x:20);
    z = (wall_z :delta_x:wall_z +40);
    [X,Y,Z] = meshgrid(x,y,z);

    c0 = 346*1000;
    f = 40000;
    omega = 2*pi*f;
    k = omega/c0;

    %ピストンの半径
    a =4.5;
    %アレイ半径
    zahyo = load('./zahyo/20200720_180.mat');
    R = zahyo.R;

    %スピーカ位置
    sp_x = zahyo.X;
    sp_y = zahyo.Y;
    sp_z = zahyo.Z;
    %鏡z座標
    im_z = wall_z-abs(sp_z-wall_z);
    %縦に並ぶトランデューサの数
    theta_sp_num = 8;

    w = zeros(theta_sp_num,1);

    CP = load("CP3.mat");
    CP_p = reshape(CP.P,len,1);
    CP_X = reshape(CP.cx,len,1);
    CP_Y = reshape(CP.cy,len,1);
    CP_Z = reshape(CP.cz,len,1)+tp;
    G = zeros(length(CP_Z(CP_Z>=wall_z)),8);
    D = zeros(length(CP_Z(CP_Z>=wall_z)),1);
    wet = eye(length(CP_Z(CP_Z>=wall_z)),length(CP_Z(CP_Z>=wall_z)));

    ind = 0;
    for c_n = 1:len
        if CP_Z(c_n)>=wall_z
            ind = ind + 1;
            g = zeros(8,1);
            p_n = 1;

            if (abs(tp-CP_Z(c_n))< 10) && (abs(CP_Y(c_n))< 10) && (abs(CP_X(c_n))< 10)
                wet(ind,ind) =wet(ind,ind) + 5;
                if (abs(tp-CP_Z(c_n))< 5) && (abs(CP_Y(c_n))< 5) && (abs(CP_X(c_n))< 5)
                    wet(ind,ind) =wet(ind,ind) + 5;
                end
            end


            for n = 1:length(sp_x)

                    xx = 0;
                    if reverse  == 1
                        if sp_y(n) > 0 
                            xx =1;
                        end
                    end

                    g(p_n) = g(p_n)+theory_p_one(k,a,CP_X(c_n),CP_Y(c_n),CP_Z(c_n),sp_x(n),sp_y(n),sp_z(n),0)*exp(1j*(pi*xx));

                    if reflect_on == 1
                        g(p_n) = g(p_n)+theory_p_one(k,a,CP_X(c_n),CP_Y(c_n),CP_Z(c_n),sp_x(n),sp_y(n),im_z(n),wall_z*2)*exp(1j*(pi*xx));
                    end


                    if n < length(sp_x) && (sp_z(n) ~= sp_z(n+1))
                        p_n = p_n +1;
                        if p_n == theta_sp_num+1
                            p_n = 1;
                        end
                    end
            end
            D(ind) = real(CP_p(c_n))+imag(CP_p(c_n))*1j;
    %         D(ind) = CP_p(c_n);
            G(ind,:) = reshape(g,1,8);
        end


    end

    A = G'*wet*G;
    B =G'*wet*D;

    beta2 =eigs(A,1);
    beta2 = beta2*10^(-10);
    I=eye(size(A));

    w=(A+0)\B;
    %振幅正規化
    w = w/max(abs(w))*1.5;
%     for l = 1:theta_sp_num
%         if abs(w(l))>1.5
%             w(l) = w(l)/(abs(w(l)));
%         end
%     end

    if save_w == 1
        sin_A = abs(w);
        phix = angle(w);

        save(sprintf('./phase/210222/LSw-25%.1f.mat', tp),'sin_A','phix');
    end

    
    if write_graph == 1
        figure(1)
        hold on
        p_n = 1;
        P = zeros(size(X));
        %xmlファイル作成のための
        xms = zeros(length(sp_x),8);
        ab = sqrt(sp_x.^2+sp_y.^2+sp_z.^2);

        x_r = -sp_x;
        y_r = -sp_y;
        z_r = -sp_z;
        ang = zeros(size(sp_x));
        pow = zeros(size(sp_x));
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


                P = P+(P0+P_im)*exp(1j*(pi*xx))*w(p_n);
                ISO = mod(angle(w(p_n))+pi*xx,2*pi);
                ang(n) = ISO;
                pow(n) = abs(w(p_n));
                 if pos_mark ==1

                     tmp = mod(ISO,2*pi);

                     if tmp < pi
                         qc = [1-tmp/pi,0,0]; %赤
                     else
                         qc = [(tmp-pi)/pi,(tmp-pi)/pi,1];  %青
                     end

                     q = quiver3(sp_x(n),sp_y(n),sp_z(n),-sp_x(n),-sp_y(n),-sp_z(n),0.05);

                     q.LineWidth = 5+abs(w(p_n))*5;
                     q.Color = qc;

                 end

                if n < length(sp_x) && (sp_z(n) ~= sp_z(n+1))
                    p_n = p_n +1;
                    if p_n == theta_sp_num+1
                        p_n = 1;
                    end
                end
        end


        Power = 20*log10(abs(P));
        xslice =0;
        yslice =[];
        zslice = wall_z;

        U = poten_cal(P,delta_x,delta_y,delta_z,c0,omega);

        % % 最大値1に正規化
        % U = U/max(max(max(abs(U))));

        slice(X,Y,Z,Power,xslice,yslice,zslice)
        % caxis([-0.04,0.04])
        view(90,0)
        xlabel('x (mm)');
        ylabel('y (mm)');
        zlabel('z (mm)');
        ax = gca;
        ax.FontSize = 24;
%         title("Potential field")
        %title("Amplitude field")
        shading interp
        c = colorbar;
        caxis([-10 10])
%         c.Label.String = 'The Gor’kov potential';
        c.Label.String = 'Sound pressure level(dB)';
        point = plot3(0,0,tp,'o','Color','w','MarkerSize',8,'MarkerFaceColor',[0.8,0.8,0.8]);
        colormap hot
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


        hold off
        if save_graph == 1
            saveas(gcf,sprintf('.//LSw-25%.1f.png', tp))
        end
        
        delete(point)
    end
 
end
% figure(2)
% slice(X,Y,Z,Power,xslice,yslice,zslice)
% view(90,0)
% ax = gca;
% ax.FontSize = 30;
% xlabel('x (mm)','FontSize',30);
% ylabel('y (mm)','FontSize',30);
% zlabel('z (mm)','FontSize',30);
% title("Amplitude field")
% c = colorbar;
% c.Label.String = 'Sound pressure level(dB)';
% shading interp
% axis equal
% caxis([-10 10])
% %  L = 6*del2(U);
% %  L_cent = reshape(L(21,21,:),length(z),1);
% %  figure(3)
% %  plot(L_cent)
%csvファイルの保存
%振幅位相の保存
% filename = 'w-25_tp-25.csv';
% table1 = table(sp_x,sp_y,sp_z,x_r,y_r,z_r,ang,pow);
% writetable(table1,filename);
% 
% %音圧の保存
% filename = 'dataw-25_tp-25.csv';
% table2 = table(reshape(X,length(x)*length(y)*length(z),1),reshape(Y,length(x)*length(y)*length(z),1),reshape(Z,length(x)*length(y)*length(z),1),reshape(Power,length(x)*length(y)*length(z),1));
% writetable(table2,filename);