close all
delta_x = 0.5;
delta_y = 0.5;
delta_z = 0.5;
index = 0;
for tp_z = -21:1:-21
    %トラップ位置Z
    tp = tp_z;
    %力表示するかどうか
    force_on = 0;
    %アレイ位置表示
    pos_mark = 0;
    %反射有りかどうか
    reflect_on = 1;
    %2次反射有りかどうか(仮)
    reflect2_on =0;
    %位相反転するかどうか
    reverse = 1;
    %壁の位置
    wall_z = -25;
    %振幅, 位相保存するかどうか
    save_w = 0;
    %グラフ保存するかどうか
    save_graph = 0;
    %グラフ表示するかどうか
    write_graph = 1;
    %音圧保存するか
    save_csv = 0;
    %重り分散
    dd = 100;
    delta = [dd,0,0;0,dd,0;0,0,dd];

    x = (-20:delta_x:20);
    y = (-20:delta_x:20);
    %z = (wall_z :delta_x:wall_z +40);
    z = (-25 :delta_x:15);
    [X,Y,Z] = meshgrid(x,y,z);

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
    A = 15;
    %スピーカ位置
    sp_x = zahyo.X;
    sp_y = zahyo.Y;
    sp_z = zahyo.Z;
    %鏡z座標
    im_z = wall_z-abs(sp_z-wall_z);
    im_z2 = (-1)*im_z + 30;
    %縦に並ぶトランデューサの数
    theta_sp_num = 7;

    w = zeros(theta_sp_num,1);

    CP = load("CP3.mat");
    CP_p = reshape(CP.P,len,1);
    CP_X = reshape(CP.cx,len,1);
    CP_Y = reshape(CP.cy,len,1);
    CP_Z = reshape(CP.cz,len,1)+tp;
    G = zeros(length(CP_Z(CP_Z>=wall_z)),theta_sp_num);
    D = zeros(length(CP_Z(CP_Z>=wall_z)),1);
    weight = eye(length(CP_Z(CP_Z>=wall_z)),length(CP_Z(CP_Z>=wall_z)));

    ind = 0;
    for c_n = 1:len
        if CP_Z(c_n)>=wall_z
            ind = ind + 1;
            g = zeros(theta_sp_num,1);
            p_n = 1;

            tmp = [CP_X(c_n);CP_Y(c_n);CP_Z(c_n)];
            myu = [0;0;tp_z];
            weight(ind,ind) =1/sqrt(2*pi)^(3/2)/sqrt(det(delta))*exp(-1/2*(tmp-myu)'*delta^(-1)*(tmp-myu));


            for n = 1:length(sp_x)
                    A = 15;
                    xx = 0;
                    if reverse  == 1
                        if sp_y(n) > 0 
                            xx =1;
                        end
                    end

                    g(p_n) = g(p_n)+P00*A*theory_p_one(k,a,CP_X(c_n),CP_Y(c_n),CP_Z(c_n),sp_x(n),sp_y(n),sp_z(n),0)*exp(1j*(pi*xx));

                    if reflect_on == 1
                        g(p_n) = g(p_n)+P00*A*theory_p_one(k,a,CP_X(c_n),CP_Y(c_n),CP_Z(c_n),sp_x(n),sp_y(n),im_z(n),wall_z*2)*exp(1j*(pi*xx));
                    end


                    if n < length(sp_x) && (sp_z(n) ~= sp_z(n+1))
                        p_n = p_n +1;
                        if p_n == theta_sp_num+1
                            p_n = 1;
                        end
                    end
            end
           % D(ind) = real(CP_p(c_n))+imag(CP_p(c_n))*1j;
            D(ind) = CP_p(c_n);
            G(ind,:) = reshape(g,1,theta_sp_num);
        end


    end

    A = G'*weight*G;
    B =G'*weight*D;

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

        save(sprintf('./phase/210629/LSGw-25_min%.1f.mat', tp),'sin_A','phix');
    end

    
    if write_graph == 1
        figure(1)
        hold on
        p_n = 1;
        P = zeros(size(X));
        %xmlファイル作成のための
        xms = zeros(length(sp_x),theta_sp_num);
        ab = sqrt(sp_x.^2+sp_y.^2+sp_z.^2);

        x_r = -sp_x;
        y_r = -sp_y;
        z_r = -sp_z;
        ang = zeros(size(sp_x));
        pow = zeros(size(sp_x));
        for n = 1:length(sp_x)
                A = 15;

                xx = 0;
                if reverse  == 1
                    if sp_y(n) < 0 
                        xx =1;
                    end
                end
                
                P_im = 0;
                P_im2 = 0;
                P0 = theory_p(k,a,X,Y,Z,sp_x(n),sp_y(n),sp_z(n),0);

                if reflect_on ==1
                    P_im = theory_p(k,a,X,Y,Z,sp_x(n),sp_y(n),im_z(n),wall_z*2);
                end
                if reflect2_on == 1
                    P_im2 = theory_p(k,a,X,Y,Z,sp_x(n),sp_y(n),im_z2(n),30);
                end

                P = P+P00*A*(P0+P_im+P_im2)*exp(1j*(pi*xx))*w(p_n);
                ISO = mod(angle(w(p_n))+pi*xx,2*pi);
                ang(n) = ISO;
                pow(n) = abs(w(p_n));
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
        zslice = tp_z;

        U = poten_cal(P,delta_x,delta_y,delta_z,c0,omega);

        % % 最大値1に正規化
        % U = U/max(max(max(abs(U))));
        %発散
        L = 6*del2(U);

        slice(X,Y,Z,L,xslice,yslice,zslice)
        
   
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
        c = colorbar;
        caxis([-0.01 0.01])
%         c.Label.String = 'The Gor’kov potential';
        %c.Label.String = 'Sound pressure level(dB)';
        hold on
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
        axis off

        hold off
        if save_graph == 1
            saveas(gcf,sprintf('./210621/LSGw-25_%.1f.png', tp))
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
        index = index + 1;
    end
    

 
end