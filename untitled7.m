P2 = zeros(size(P));
P2(,1:20,1:4) = 0;
P2(:,22:end,1:4) = 0;
Power = 20*log10(abs(P2));

xslice =0;
yslice =[];
zslice = 0;

U = poten_cal(P2,delta_x,delta_y,delta_z,c0,omega);
figure(5)
slice(X,Y,Z,U,xslice,yslice,zslice)

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
c.Label.String = 'The GorÅfkov potential';
axis equal


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
caxis([-10 5])
shading interp
axis equal

 L = 6*del2(U);
 L_cent = reshape(L(21,21,:),length(z),1);
 figure(4)
 
 plot(L_cent)