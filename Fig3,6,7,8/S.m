function Spin=S(theta_1,theta_2,phi_1,phi_2)
    %----------------------------------------------------------------------
    %Function: Change the sphere cordination into complex Cartisian array
    %Input: (theta_1,theta_2,phi_1,phi_2) is the Marsaglia of spin at one site
    %Output: An array 'Spin', three complex components are the Cartisian one.
    %----------------------------------------------------------------------
    d_xR=theta_2^(1/4)*theta_1^(1/2)*sin(phi_1);
    d_xI=theta_2^(1/4)*theta_1^(1/2)*cos(phi_1);
    d_yR=theta_2^(1/4)*sqrt(1-theta_1)*sin(phi_2);
    d_yI=theta_2^(1/4)*sqrt(1-theta_1)*cos(phi_2);
    d_z=sqrt(1-theta_2^(1/2))+0*1i;
    d_x=d_xR+1i*d_xI;
    d_y=d_yR+1i*d_yI;
    Spin=[d_x;d_y;d_z];
end