function E=Energy_per_site(i,j,CellSize,spin,Parameter)
    %----------------------------------------------------------------------
    %Function: Compute the Energy per each site under the Energy functional
    %Input: i,j is the index of sites.
%           CellSize: 2 times 1 array [i_site_num, j_site_num]
    %       spin: The array to storage spin variables in the form of site.
%           Parameter=[J_1,J_2,J_3]  
    %Output: Average energy of this site for specific 'spin' configuration
    %----------------------------------------------------------------------
    i_site_num=CellSize(1);
    j_site_num=CellSize(2);
    J_1=Parameter(1);
    J_2=Parameter(2);
    J_3=Parameter(3);
    i_prev=neighborindex(i,-1,i_site_num);
    j_prev=neighborindex(j,-1,j_site_num);
    i_next=neighborindex(i,1,i_site_num);
    j_next=neighborindex(j,1,j_site_num);

%     i_prev = int32(mod(i - 2, i_site_num))+int32(1);
%     j_prev = int32(mod(j - 2, j_site_num))+int32(1);
%     i_next = int32(mod(i, i_site_num))+int32(1);
%     j_next = int32(mod(j, j_site_num))+int32(1);
    s_x=spin(i,j,1);
    s_y=spin(i,j,2);
    s_z=spin(i,j,3);
    d=[s_x;s_y;s_z];
    %Interaction-term
    E=0;
    t_x=[spin(i_next,j,1),spin(i,j_next,1),spin(i_prev,j_next,1),spin(i_prev,j,1),spin(i,j_prev,1),spin(i_next,j_prev,1)];
    t_y=[spin(i_next,j,2),spin(i,j_next,2),spin(i_prev,j_next,2),spin(i_prev,j,2),spin(i,j_prev,2),spin(i_next,j_prev,2)];
    t_z=[spin(i_next,j,3),spin(i,j_next,3),spin(i_prev,j_next,3),spin(i_prev,j,3),spin(i,j_prev,3),spin(i_next,j_prev,3)];
    for l=1:6
        %Vector notation
        d_n=[t_x(l);t_y(l);t_z(l)];
        %SU(3) term
%         for m=1:8
%             E=E+J_3*(d'*lambda(:,:,m)*d)*(d_n'*lambda(:,:,m)*d_n);
%         end
        E=E+2*abs(d'*d_n)^2-2/3;
        %SU(2) term
%         for m=1:3
%             E=E+J_2*(d'*lambda(:,:,m)*d)*(d_n'*lambda(:,:,m)*d_n);
%         end
        d_xy=conj(d(1))*d(2);
        d_xyn=conj(d_n(1))*d_n(2);
        d_xy_R=2*real(d_xy);
        d_xyn_R=2*real(d_xyn);
        d_xy_I=2*imag(d_xy);
        d_xyn_I=2*imag(d_xyn);
        E=E+J_2*(d_xy_R*d_xyn_R+d_xy_I*d_xyn_I);
        E=E+J_2*(abs(d(1))^2-abs(d(2))^2)*(abs(d_n(1))^2-abs(d_n(2))^2);
        %U(1) term
        E=E+J_1*(1-3*abs(d(3))^2)*(1-3*abs(d_n(3))^2)/3;
    end
    E=E/2;
end
