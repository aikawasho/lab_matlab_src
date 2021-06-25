figure(2)
hold on
c1 = -22:4:0;
c2 = -22:4:0;
cLS = -22:4:0;
n = 0;
for tp_z = -2:-4:-22
    n = n+ 1;
%     con_yz1 = load(sprintf('./210120/abs_f_yz_%.1f.mat', tp_z));
    con_yz2 = load(sprintf('./210121/abs_f_yz_%.1f.mat', tp_z));
    con_yzLS = load(sprintf('./210120/LSabs_f_yz_%.1f.mat', tp_z));
    
%     con_xy1 = load(sprintf('./210120/abs_f_xy_%.1f.mat', tp_z));
    con_xy2 = load(sprintf('./210121/abs_f_xy_%.1f.mat', tp_z));
    con_xyLS = load(sprintf('./210120/LSabs_f_xy_%.1f.mat', tp_z));
%     c1(n) =  sum(sum(con_yz1.abs_f_yz+con_xy1.abs_f_xy));
    c2(n) =  sum(sum(con_yz2.abs_f_yz+con_xy2.abs_f_xy));
    cLS(n) =  sum(sum(con_yzLS.abs_f_yz+con_xyLS.abs_f_xy));
    
end


plot( 2:4:22,cLS,'LineWidth',2,'DisplayName','LS');
% plot(-22:4:0,c1,'LineWidth',2,'DisplayName','2');
plot( 2:4:22,c2,'LineWidth',2,'DisplayName','QP');
xlabel('distance from the array (mm)');
ylabel('Acoustic Radiation Force');
ax = gca;
ax.FontSize = 15;
hold off
legend show