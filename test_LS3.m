delta_x = 1;
delta_y = 1;
delta_z = 1;
%トラップ位置Z
tp = -0;
%力表示するかどうか
force_on = 0;
%反射有りかどうか
reflect_on = 1;
%位相反転するかどうか
reverse = 1;
%壁の位置
wall_z = -15;
x = (-20:delta_x:20);
y = (-20:delta_x:20);
z = (wall_z :delta_x:20);
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
len = 8;
CP_p = reshape(CP.P,len,1);
CP_X = reshape(CP.cx,len,1);
CP_Y = reshape(CP.cy,len,1);
CP_Z = reshape(CP.cz,len,1)+tp;
G = zeros(length(CP_Z(CP_Z>=wall_z)),8);
D = zeros(length(CP_Z(CP_Z>=wall_z)),1);
ind = 0;
for c_n = 1:len
    if CP_Z(c_n)>=wall_z
        ind = ind + 1;
        g = zeros(8,1);
        p_n = 1;
        
        if (abs(tp-CP_Z(c_n))< 5) && (abs(CP_Y(c_n))< 5) && (abs(CP_X(c_n))< 5)
            wet(ind,ind) = 10;
            
        end
            
        for n = 1:length(sp_x)

                xx = 0;
                if reverse  == 1
                    if sp_y(n) < 0 
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
        D(ind) = CP_p(c_n);
        G(ind,:) = reshape(g,1,8);
    end
        

end
beta2 =eigs(G,1);
beta2 = beta2*10^(-10);
I=eye(size(G));
w=(G+beta2*I)\D;
%振幅正規化
w = w/max(abs(w));

p_n = 1;
 
P = zeros(size(X));
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

        if n < length(sp_x) && (sp_z(n) ~= sp_z(n+1))
            p_n = p_n +1;
            if p_n == theta_sp_num+1
                p_n = 1;
            end
        end
   
end
Power = 20*log10(abs(P));
figure(5)
xslice =0;
yslice =[];
zslice = 0;

U = poten_cal(P,delta_x,delta_y,delta_z,c0,omega);

% % 最大値1に正規化
% U = U/max(max(max(abs(U))));

slice(X,Y,Z,U,xslice,yslice,zslice)
hold on
caxis([-0.04,0.04])
view(90,0)
xlabel('x (mm)','FontSize',15);
ylabel('y (mm)','FontSize',15);
zlabel('z (mm)','FontSize',15);
ax = gca;
ax.FontSize = 15;
title("Potential field")
shading interp
c = colorbar;
c.Label.String = 'The Gor’kov potential';
if force_on == 1
    [F_x0, F_y0, F_z0] =  gradient(U,delta_x,delta_y,delta_z);
    span = 2;
    F_x0 = -F_x0(1:span:end,1:span:end,1:span:end); F_y0 = -F_y0(1:span:end,1:span:end,1:span:end); F_z0 = -F_z0(1:span:end,1:span:end,1:span:end);

    F_x0 = F_x0/max(max(max(F_x0)));
    F_x0(abs(F_x0)<0.4) = 0;

    F_y0 = F_y0/max(max(max(F_y0)));
    F_y0(abs(F_y0)<0.4) = 0;


    F_z0 = F_z0/max(max(max(F_z0)));
    F_z0(abs(F_z0)<0.4) = 0;

    quiver3(X(1:span:end,1:span:end,1:span:end),Y(1:span:end,1:span:end,1:span:end),Z(1:span:end,1:span:end,1:span:end),real(F_x0),real(F_y0),real(F_z0),1,"red")
end

hold off

figure(6)
slice(X,Y,Z,Power,xslice,yslice,zslice)
view(90,0)
ax = gca;
ax.FontSize = 30;
xlabel('x (mm)','FontSize',30);
ylabel('y (mm)','FontSize',30);
zlabel('z (mm)','FontSize',30);
title("Amplitude field")
c = colorbar;
c.Label.String = 'Sound pressure level(dB)';
shading interp
axis equal

 L = 6*del2(U);
 L_cent = reshape(L(21,21,:),length(z),1);
 figure(7)
 plot(L_cent)