Sum_P = reshape(sum(P,2),length(z),length(y));
figure(1)
surf(z,y,20*log10(abs(Sum_P)))
view(-90,90)
axis equal
shading interp
ax = gca;
ax.FontSize = 15;
xlabel('z (mm)');
ylabel('y (mm)');
title("Integral Amplitude Field")
%colormap jet
caxis([0,10])
hold on
point = plot(tp(3),tp(2),'o','Color','w','MarkerSize',8,'MarkerFaceColor',[0.8,0.8,0.8]);

hold off