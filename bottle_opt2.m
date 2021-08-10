close all
global wall_z
global omega
global c0
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
global theta_sp_num
global reflect_on
global P00
global use_channels
wp = 1;
wx = 1;
wy = 1;
wz = 1;
delta_x = 1;
delta_y = 1;
delta_z = 1;
index = 1;
%  v = VideoWriter('./optbottle123468.mp4','MPEG-4');
%  v.Quality = 100;
%  open(v);
%  loops = length( -24:0.05:-20);
%  F(loops) = struct('cdata',[],'colormap',[]);
for tp_z =  -24:0.1:-24
    %3D計算
    yes_3D = 1;
    %トラップ位置Z
    tp = [0,0,tp_z];
    %力表示するかどうか
    force_on = 0;
    %アレイ位置表示
    pos_mark = 0;
    %反射有りかどうか
    reflect_on = 1;
    %2次反射有りかどうか(仮)
    reflect2_on =0;
    %位相反転するかどうか
    reverse = 0;
    %壁の位置
    wall_z = -25;
    %位相読み込む
    load_on = 0;
    %振幅, 位相保存するかどうか
    save_phi = 1;
    %グラフ保存するかどうか
    save_graph = 0;
    %グラフ表示するかどうか
    write_graph = 1;
    %音圧保存するか
    save_csv = 0;
    
    %重り分散
    dd = 20;
    delta = [dd,0,0;0,dd,0;0,0,dd];

    x = (-20:delta_x:20);
    y = (-20:delta_x:20);
    %z = (wall_z :delta_x:wall_z +40);
    z = (-25 :delta_x:15);

    [X,Y,Z] = meshgrid(x,y,z);
%     x = (-10:delta_x:10);
%     y = (-10:delta_y:10);
%     z = (tp_z-10:delta_z:tp_z+10);
%     z_dis = (-10:10);
    [Y1,Z1] = meshgrid(y,z);
%     [Y2,X2] = meshgrid(y,x);
    %[Y00,Z_dis] = meshgrid(y,z_dis);

    c0 = 346*1000;
    f = 40000;
    omega = 2*pi*f;
    k = omega/c0;

    %ピストンの半径
    a =4.5;
    %アレイ半径
    zahyo = load('./zahyo/20200720_180.mat');
    R = 64.1;
    P00 = 0.17;
    %スピーカ位置
    sp_x = zahyo.X;
    sp_y = zahyo.Y;
    sp_z = zahyo.Z;
    %鏡z座標
    im_z = wall_z-abs(sp_z-wall_z);
    im_z2 = (-1)*im_z + 30;
    %縦に並ぶトランデューサの数
    theta_sp_num = 7;
     
    use_channels = [1,2,3,4,5,6,7];
    if load_on == 0
        %phi0 = zeros(theta_sp_num,1);
        phi0 = ones(theta_sp_num,2);

        options = optimoptions(@fminunc,'Display','iter','FunvalCheck','on','MaxFunctionEvaluations',3000,'PlotFcn','optimplotfval');
        % filename = 'PlotFcns0818_1'+string(tp)+'.fig';
        fun = @bottle_ob_fun;
        [phix,fval] = fminunc(fun,phi0,options);

        % saveas(gcf,filename)
    elseif load_on == 1
        phix = load(file_name);
        phix = phix.phix;
    else
        phix = zeros(theta_sp_num,1);
    end

    if save_phi == 1
        save(sprintf('./phase/210726/op2tbottle123468-25%.1f.mat', tp(3)),'phix');
    %     save('./phase/20201013/gradation.mat','phix');
    end

    
    if write_graph == 1
        

        figure(1)
        hold on
        p_n = 1;
        if yes_3D  == 1 
            P = zeros(size(X));
        end
        %xmlファイル作成のための
        xms = zeros(length(sp_x),8);
        ab = sqrt(sp_x.^2+sp_y.^2+sp_z.^2);

        x_r = -sp_x;
        y_r = -sp_y;
        z_r = -sp_z;
        ang = zeros(size(sp_x));
        pow = zeros(size(sp_x));
        
        P_yz = zeros(size(Y1));
        %P_xy = zeros(size(X2));
        qc = zeros(length(sp_x),3);
        for n = 1:length(sp_x)
                A = 0;

                xx = 0;
                if reverse  == 1
                    if sp_y(n) < 0 
                        xx =1;
                    end
                end
                    
                if  ismember(p_n,use_channels) 
                    A = 15;
                end
                
                P_im = 0;
                P_im_1 = 0;
%                 P_im_2 = 0;
                
                if yes_3D  == 1 
                    P0 = theory_p(k,a,X,Y,Z,sp_x(n),sp_y(n),sp_z(n),0);
                end
                P0_1 = theory_p_2d(k,a,0,Y1,Z1,sp_x(n),sp_y(n),sp_z(n),0);
%                 P0_2 = theory_p_2d(k,a,X2,Y2,tp_z,sp_x(n),sp_y(n),sp_z(n),0);
                
                if reflect_on ==1
                    if yes_3D  == 1 
                        P_im = theory_p(k,a,X,Y,Z,sp_x(n),sp_y(n),im_z(n),wall_z*2);
                    end
                    P_im_1 = theory_p_2d(k,a,0,Y1,Z1,sp_x(n),sp_y(n),im_z(n),wall_z*2);
%                     P_im_2 = theory_p_2d(k,a,X2,Y2,tp_z,sp_x(n),sp_y(n),im_z(n),wall_z*2);
                end
                ISO = phix(p_n,1)+pi*xx;
                if yes_3D  == 1 
                 P = P+P00*A*abs(1.5*sin(phix(p_n,2)))*(P0+P_im)*exp(1j*(pi*xx))*exp(1j*ISO);
                end
%                 P_yz = P_yz+P00*A*(P0_1+P_im_1)*exp(1j*(pi*xx))*w(Gw_index(p_n));

                P_yz = P_yz+A*P00*abs(1.5*sin(phix(p_n,2)))*(P0_1+P_im_1)*exp(1j*ISO);
                %P_yz = P_yz+A*P00*(P0_1+P_im_1)*exp(1j*ISO);
%                 P_xy = P_xy+P00*A*(P0_2+P_im_2)*exp(1j*(pi*xx))*w(Gw_index(p_n));
                
                ang(n) = ISO;
                if yes_3D  == 1 

                 if pos_mark ==1

                     tmp = mod(ISO,2*pi);

                     %if tmp < pi
%                          qc = [1-tmp/pi,0,0]; %赤
%                      else
%                          qc = [(tmp-pi)/pi,(tmp-pi)/pi,1];  %青
%                      end

                     q = quiver3(sp_x(n),sp_y(n),sp_z(n),-sp_x(n),-sp_y(n),-sp_z(n),0.05);

                     %q.LineWidth = 5+abs(w(p_n))*5;
                     q.LineWidth = 5;
                     q.Color = [0,0,0];
                     
                     q1 = quiver3(sp_x(n),sp_y(n),im_z(n),-sp_x(n),-sp_y(n),-im_z(n),0.05);

                     q1.LineWidth = 5;
                     q1.Color = [0.6,0.6,0.6];
                     
                     q2 = quiver3(sp_x(n),sp_y(n),im_z2(n),-sp_x(n),-sp_y(n),-im_z2(n),0.05);

                     q2.LineWidth = 5;
                     q2.Color = [0.6,0.6,0.6];

                 end
                end

                if n < length(sp_x) && (sp_z(n) ~= sp_z(n+1)) && (  sp_z(n+1) > 69 || sp_z(n+1) < 65 )                        
                    p_n = p_n +1;
                    if p_n == theta_sp_num+1
                        p_n = 1;
                    end
                end
        end

        if yes_3D  == 1 

            Power = 20*log10(abs(P));
            U = poten_cal(P,delta_x,delta_y,delta_z,c0,omega);
                    %発散
            L = 6*del2(U);
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
            %point = plot3(tp,'o','Color','w','MarkerSize',3,'MarkerFaceColor',[0.8,0.8,0.8]);
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
%             F(index) = getframe(gca);
%             for fff = 1:5
%                 writeVideo(v,F(index));
%             end
            delete(point)
            hold off
            if save_graph == 1
                saveas(gcf,sprintf('./210621/LSGw-25_%.1f.png', tp(3)))
            end

            %delete(point)
            weight_graph = zeros(size(X));    
            for i = 1:length(x)
                for j = 1:length(y)
                    for l = 1:length(z)
                        tmp = [x(i);y(j);z(l)];
                        myu = [0;0;tp_z];
                        weight_graph(j,i,l) = 1/sqrt(2*pi)^(3/2)/sqrt(det(delta))*exp(-1/2*(tmp-myu)'*delta^(-1)*(tmp-myu));
                    end
                end
            end

            figure(2)
            slice(X,Y,Z,weight_graph,xslice,yslice,zslice)
            % caxis([-0.04,0.04])
            view(90,0)
            xlabel('x (mm)');
            ylabel('y (mm)');
            zlabel('z (mm)');
            ax = gca;
    %         title("Potential field")
            %title("Amplitude field")
            shading interp
            % %音圧の保存
            if save_csv == 1
                filename = sprintf('./csv/210706/LSGW-25_%.0f.csv', index);
                table2 = table(reshape(X/delta_x,length(x)*length(y)*length(z),1),reshape(Y/delta_y,length(x)*length(y)*length(z),1),reshape(Z/delta_z,length(x)*length(y)*length(z),1),reshape(L,length(x)*length(y)*length(z),1));
                writetable(table2,filename);
            end
        end

        xslice =0;
        yslice =[];
        zslice = tp_z;

        U_yz = poten_cal_2d(P_yz,delta_y,delta_z,c0,omega);
%         U_xy = poten_cal_2d(P_xy,delta_x,delta_z,c0,omega);
        P_yzower = 20*log10(abs(P_yz));
        %蚊にかかる重力24.5 mN
        %半径1.5mmのEPSにかかる重力 0.004 mN
        L_yz = 6*del2(U_yz);
%         L_xy = 6*del2(U_xy);
        [F_y1, F_z1] = gradient(U_yz,delta_y,delta_z);
        F_y1((Z1<= wall_z)) = 0;
        F_z1((Z1<= wall_z)) = 0;
%         U_yz((Z1<= wall_z)) = 0;
%         [F_x2, F_z2] = gradient(U_xy,delta_y,delta_z);

        % % 最大値1に正規化
        % U = U/max(max(max(abs(U))));


       
        
        
            figure(3)
            surf(Y1,Z1,P_yzower)%,'FaceAlpha',0.5)
            shading interp
            %quiver(Y00,Z_dis,-F_y1,-F_z1,1.1,'Color',[0,0,0])
            hold on 
            point1 = scatter(tp(2),tp(3),30,'MarkerFaceColor',[0.8,0.8,0.8]);%,'MarkerFaceAlpha',0.5);
            %quiver(Y,Z,F_y,F_z)
%             xlim([min(y),max(y)])
%             ylim([min(z_dis),max(z_dis)])
            ax = gca;
            ax.FontSize = 15;
           % xlabel('x-axis displacement of desired trap position (mm)');
           % ylabel('z-axis displacement of desired trap position (mm)');
            xlabel('x (mm)');
            ylabel('z (mm)');
            %title('Acoustic Radiation Force x-z plane')
            c = colorbar;
            c.Label.String = 'Sound pressure level(dB)';
        %         c.Label.String  = 'divergence of potential field';
            caxis([-10 10])
            view(0,90)
            colormap jet
            axis equal
%             F(index) = getframe(gca);
%             for fff = 1:10
%                 writeVideo(v,F(index));
%             end
            delete(point1)
            hold off
    index = index + 1;
        if save_graph == 1
            filename = strcat(name,'xz%.1f.png');
            saveas(gcf,sprintf(filename, tp(3)))
        end

%         figure(4)
%             surf(Y1,X2,U_xy,'FaceAlpha',0.5)
%             shading interp
%             hold on 
%             %quiver(Y00,Z_dis,-F_x2,-F_z2,1.1,'Color',[0,0,0])
% 
%             point2 = scatter(0,0,250,'MarkerEdgeColor',[0.8,0.8,0.8],'MarkerFaceColor',[0.8,0.8,0.8],'MarkerFaceAlpha',0.5);
%             %quiver(Y,Z,F_y,F_z)
%             xlim([min(y),max(y)])
%             ylim([min(x),max(x)])
%             ax = gca;
%             ax.FontSize = 15;
%             xlabel('x-axis displacement of desired trap position (mm)');
%             ylabel('y-axis displacement of desired trap position (mm)');
%             title('Acoustic Radiation Force x-y plane')
%             c = colorbar;
%             c.Label.String = 'The Gor’kov potential';
%         %         c.Label.String = 'divergence of potential field';
%             caxis([0 1])
%             colormap jet
%             view(0,90)
%             axis equal
%             hold off
%         if save_graph == 1
%             filename = strcat(name,'xy%.1f.png');
%             saveas(gcf,sprintf(filename, tp))
%         end
        
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
%  table1 = table(sp_x,sp_y,sp_z,x_r,y_r,z_r,ang,pow);
% writetable(table1,filename);
% 
close(v);