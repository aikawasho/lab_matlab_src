delta_x = 1;
delta_y = 1;
delta_z = 1;
x = (-100:delta_x:100);
y = (-100:delta_x:100);
z = (-100:delta_x:100);
[X,Y,Z] = meshgrid(x,y,z);
c = 340*1000;
f = 40000;
omega = 2*pi*f;
k = omega/c;

%ピストンの半径
a =4.5;
%スピーカ位置
sp_pos = [0,0,-100];

P1 = make_p(k,a,X,Y,Z,sp_pos(1),sp_pos(2),sp_pos(3),10);
P2 = theory_p(k,a,X,Y,Z,sp_pos(1),sp_pos(2),sp_pos(3),10);

P1 = P1/max(max(max(abs(P1))))*100;
Power1 = 20*log10(abs(P1));
P2 = P2/max(max(max(abs(P2))))*100;
Power2 = 20*log10(abs(P2));

xslice =0;
yslice =0;
zslice = 0;

figure(1)
slice(X,Y,Z,Power1,xslice,yslice,zslice)
% colormap(jet);
view(90,0)
ax = gca;
ax.FontSize = 30;
xlabel('x (mm)','FontSize',30);
ylabel('y (mm)','FontSize',30);
zlabel('z (mm)','FontSize',30);
title("Amplitude field")
c = colorbar;
c.Label.String = 'Sound pressure level(dB)';
caxis([-40,40])
shading interp

axis equal

figure(2)
slice(X,Y,Z,Power2,xslice,yslice,zslice)
% colormap(jet);
view(90,0)
ax = gca;
ax.FontSize = 30;
xlabel('x (mm)','FontSize',30);
ylabel('y (mm)','FontSize',30);
zlabel('z (mm)','FontSize',30);
title("Amplitude field")
c = colorbar;
c.Label.String = 'Sound pressure level(dB)theory';
 caxis([-40,40])
shading interp

axis equal
