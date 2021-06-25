

 %a ピストンの直径
% R アレイ半径
    a = 10;
    R = 65;
    %縦に何こスピーカが並ぶのか
    L = 2*pi*R/4;
    L2 = L-3;
    theta_num = round(L2/(a+3));
   theta = zeros(theta_num,1);
    theta(1) = pi/2/L*5+pi/2/L*L2/theta_num/2;
    for num = 2:theta_num
        theta(num) = theta(num-1)+pi/2/L*L2/theta_num;
    end
    
    phi = zeros(11,1);
    sp_pos0=zeros(11,2);
    
    ind = 1;
    for num =1:length(theta)
        R_tmp = R*sin(theta(num));
        circ = 2*pi*R_tmp;
        phi_sp_num = round(circ/4/(a+15));
        if phi_sp_num == 0
            phi_sp_num=1;
        end
            
 
        for num2 = 1:phi_sp_num
    
            step = circ/phi_sp_num;
            if num2 == 1
                phi(ind)=pi/2/circ*(step/2);
            
            else
                phi(ind) = phi(ind-1)+pi/2/circ*step;
            end
            sp_pos0(ind,1) = theta(num);
            sp_pos0(ind,2) = phi(ind);
            ind = ind +1;
            
        end
    end



    x = R*sin(sp_pos0(:,1)).*cos(sp_pos0(:,2));
    y = R*sin(sp_pos0(:,1)).*sin(sp_pos0(:,2));
    z = R*cos(sp_pos0(:,1));
    X = [x;-x;x;-x];
    Y = [y;y;-y;-y];
    Z = [z;z;z;z];
    figure(2)
    quiver3(X,Y,Z,-X,-Y,-Z,0.5)
    axis equal
    save('phase/20201027_88.mat','R','X','Y','Z');
