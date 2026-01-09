function lambda123data=From_spin_to_lambda123(spin_data)
%     ----------------------------
%     Function: Compute average value of lambda^{su2} on each sublattice
%     Input: spin_data: Spin configuration in i_site_num \times
%            j_site_num*3 format
%     Output：a i_site_num*j_site_num*3
%     ----------------------------
    i_site_num=size(spin_data,1);
    j_site_num=size(spin_data,2);
    lambda123data=zeros(i_site_num,j_site_num,3);
    for i=1:i_site_num
        for j=1:j_site_num
            d=spin_data(i,j,:);
            d_xy=2*conj(d(1))*d(2);
            lambda123data(i,j,1)=real(d_xy);
            lambda123data(i,j,2)=imag(d_xy);
            lambda123data(i,j,3)=abs(d(1))^2-abs(d(2))^2;
        end
    end
end