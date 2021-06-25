function P = theory_p_f3(k,a,X,Y,Z,sp_x,sp_y,sp_z,f_p)
% �ʑ����[��, �U��1�̕��f�����g�̉������v�Z
% k �g��
% a �X�s�[�J���a
% XYZ ����O���b�h(3�����e���\��)
% sp_x,y,x �X�s�[�J�ʒu
%f_p �A���C�œ_
    Dim = sqrt((X-sp_x).^2+(Y-sp_y).^2+(Z-sp_z).^2);
    ip1 =(X-sp_x)*(f_p(1)-sp_x)+(Y-sp_y)*(f_p(2)-sp_y)+(Z-sp_z)*(f_p(3)-sp_z);
    ip2 = Dim*sqrt((f_p(1)-sp_x).^2+(f_p(2)-sp_y).^2+(f_p(3)-sp_z).^2);
    theta = acos(ip1./ip2);
    
    dire = bessel_func2(k,a,theta);
    
    P =dire./Dim.*exp(1j*k*Dim);
    %inf��Ԗڂɑ傫���l
    %P_m2 = maxk(max(max(abs(P))),2);
    %P(P==inf) = P_m2(2);
end