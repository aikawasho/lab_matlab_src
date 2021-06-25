function U = poten_cal(P, delta_x, delta_y, delta_z, c_0, omega) %c_0�͔}���̉���, f�͉��g�̎��g��, rho_p�͔������̖̂��x, rho_0�͔}���̖��x, V_p�͔������̂̑̐�
    
    V_p = 4/3*pi*(1.5)^3; % 3mm��
    c_p = 900*1000;%�|���X�`���� 
    rho_0 = 1.3/(1000^3); %��C
    %rho_p = 1700; % ���C�|���X�`�����͖��x�Ƒ��x���킩��Ȃ��̋����Ƃ����i���Ȃ�A�o�E�g�j   
    rho_p = 29/(1000^3); % ���C�|���X�`����

    [p_x, p_y, p_z] = gradient(P,delta_x,delta_y,delta_z);  
%     
    laplacianP = abs(p_x).^2+abs(p_y).^2+abs(p_z).^2;
    

    
    %�|�e���V�����̌v�Z

    K_1 = 1/4*V_p*(1/(c_0^2*rho_0) - 1/(c_p^2*rho_p));
    K_2 = 3/4*V_p*( (rho_p-rho_0)/(omega^2*rho_0*(rho_0+2*rho_p)));
    U = 2*K_1*(abs(P).^2) - 2*K_2*laplacianP;
    
    
%     hold on 
    
    %Force�̌v�Z
    % matlab��gradient���g�p���Ă݂�
%     [F_x, F_y, F_z] =  gradient(U,delta_x,delta_y,delta_z);
%     F_x = -F_x; F_y = -F_y; F_z = -F_z;

end