function P = theory_p_2dflat(k,a,X,Y,Z,sp_x,sp_y,sp_z,muki)
% �ʑ����[��, �U��1�̕��f�����g�̉������v�Z
% k �g��
% a �X�s�[�J���a
% XYZ ����O���b�h(3�����e���\��)
% sp_x,y,x �X�s�[�J�ʒu

    Dim = sqrt((X-sp_x).^2+(Y-sp_y).^2+(Z-sp_z).^2);
    %����
    ip1 =(Z-sp_z)*(sp_z*muki);
    %�x�N�g���傫��
    ip2 = Dim.*sqrt((sp_z).^2);
    theta = acos(ip1./ip2);
    
    dire = bessel_func2(k,a,theta);
    
    P =dire./Dim.*exp(1j*k*Dim);

end