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

%トラップ位置
tp = [0,0,-10];
%力表示するかどうか
force_on = 0;
%反射有りかどうか
reflect_on = 0;
%位相反転するかどうか
reverse = 1;
%壁の位置
wall_z = -15;
%位相を読み込むかどうか
load_on = 5;
%読み込みファイルパス
file_name = './phase/20200824/now_phase-10.mat';
ph = 0:2*pi/10:2*pi-2*pi/10;
th = 0:pi/10:pi;
[phi,theta] = meshgrid(ph,th);
r0 = [0.5,1,1.5,2,2.5,3,3.5,4,4.5,5];
cx = zeros(length(r0),length(th),length(ph));
cy = zeros(size(cx));
cz = zeros(size(cy));
index= 0;
for r = r0
    index = index + 1;
    cx(index,:,:) = r^2*sin(theta).*cos(phi);
    cy(index,:,:) = r^2*sin(theta).*sin(phi);
    cz(index,:,:) = r^2*cos(theta);
end
% 
% [CP_X,CP_Y,CP_Z] = reshape(cx,60,1),reshape(cy,60,1),reshape(cz,60,1));
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

%縦に並ぶトランデューサの数
theta_sp_num = 8;
%鏡z座標
im_z = wall_z-abs(sp_z-wall_z);

if load_on == 0
    phi0 = zeros(theta_sp_num,1);
    options = optimoptions(@fmincon,'Display','iter','FunvalCheck','on','MaxFunctionEvaluations',3000,'PlotFcn','optimplotfval');
    % filename = 'PlotFcns0818_1'+string(tp)+'.fig';
    fun = @object_fun;
    [phix,fval] = fminunc(fun,phi0,options);
    % saveas(gcf,filename)
elseif load_on == 1
    phix = load(file_name);
    phix = phix.phix;
else
    phix = zeros(theta_sp_num,1);
end


ph_n = 1;
P = zeros(size(cx));
for n = 1:length(sp_x)
        xx = 0;
        if reverse  == 1
            if sp_y(n) < 0 
                xx =1;
            end
        end
        
        P_im = 0;
        P0 = theory_p(k,a,cx,cy,cz,sp_x(n),sp_y(n),sp_z(n),0);
        
        if reflect_on ==1
            P_im = theory_p(k,a,cx,cy,cz,sp_x(n),sp_y(n),im_z(n),wall_z*2);
        end
        
        ISO = phix(ph_n)+pi*xx;
        P = P+(P0+P_im)*exp(1j*ISO);

         
        if n < length(sp_x) && (sp_z(n) ~= sp_z(n+1))
            ph_n = ph_n +1;
            if ph_n == theta_sp_num+1
                ph_n = 1;
            end
        end
   
end
save("CP3.mat","P","cx","cy","cz");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%位置プロット
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delta_x = 1;
delta_y = 1;
delta_z = 1;
x = (-20:delta_x:20);
y = (-20:delta_y:20);
z = (wall_z :delta_z:20);
[X,Y,Z] = meshgrid(x,y,z);

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
        
        ISO = phix(ph_n)+pi*xx;
        P = P+(P0+P_im)*exp(1j*ISO);

         
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
zslice = 0;

U = poten_cal(P,delta_x,delta_y,delta_z,c0,omega);
slice(X,Y,Z,U,xslice,yslice,zslice)
len = length(r0)*length(th)*length(ph);
plot3(reshape(cx,len,1),reshape(cy,len,1),reshape(cz,len,1),'o')
view(90,0)
xlabel('x (mm)','FontSize',30);
ylabel('y (mm)','FontSize',30);
zlabel('z (mm)','FontSize',30);
caxis([-0.04,0.04])
ax = gca;
ax.FontSize = 15;
title("Potential field")
shading interp
c = colorbar;
c.Label.String = 'The Gor’kov potential';
axis equal

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

