delta_x = 1;
delta_y = 1;
delta_z = 1;
%トラップ位置Z
tp = -8;
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
%スピーカ位置
sp_x = zahyo.X;
sp_y = zahyo.Y;
sp_z = zahyo.Z;
%鏡z座標
im_z = wall_z-abs(sp_z-wall_z);
%縦に並ぶトランデューサの数
theta_sp_num = 8;

w = zeros(theta_sp_num,1);

CP = load("CP_test2.mat");

CP_P = zeros(7,7,7);
CP_X = zeros(7,7,7);
CP_Y = zeros(7,7,7);
CP_Z = zeros(7,7,7);

% CP_p(1,1,1) = CP.P(21,21,21);
% CP_X(1,1,1) = CP.X(21,21,21);
% CP_Y(1,1,1) = CP.Y(21,21,21);
% CP_Z(1,1,1) = CP.Z(21,21,21);

for i = -6:6
   for j = -6:6
       for l = -6:6
        CP_P(i+7,j+7,l+7) = CP.P(21+i+j+l,21+i+j+l,21+i+j+l);
        CP_X(i+7,j+7,l+7) = CP.X(21+i+j+l,21+i+j+l,21+i+j+l);
        CP_Y(i+7,j+7,l+7) = CP.Y(21+i+j+l,21+i+j+l,21+i+j+l);
        CP_Z(i+7,j+7,l+7) = CP.Z(21+i+j+l,21+i+j+l,21+i+j+l);
       end
   end
end

len = 2197;
CP_X = reshape(CP_X,len,1);
CP_Y = reshape(CP_Y,len,1);
CP_Z = reshape(CP_Z,len,1);
CP_P = reshape(CP_P,len,1);
Power = 20*log10(abs(CP.P));
figure(1)
xslice =0;
yslice =[];
zslice = 0;

slice(CP.X(),CP.Y,CP.Z,Power,xslice,yslice,zslice)
hold on 
plot3(CP_X,CP_P,CP_Z,'o')
shading interp
axis equal
hold off
    

ind = 0;
for c_n = 1:len
    if CP_Z(c_n)>=wall_z
        ind = ind + 1;
        g = zeros(8,1);
        p_n = 1;
        
        if (abs(tp-CP_Z(c_n))< 5) && (abs(CP_Y(c_n))< 5) && (abs(CP_X(c_n))< 5)
            wet(ind,ind) = 20;
            
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
        D(ind) = CP_P(c_n);
        G(ind,:) = reshape(g,1,8);
    end
        

end

A = G'*wet*G;
B =G'*wet*D;

% beta2 =eigs(A,1);
% beta2 = beta2;%*10^(-10);
% I=eye(size(A));

w=A\B;
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