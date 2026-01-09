function lambda8data=From_spin_to_lambda8(spin_data)
%     ----------------------------
%     Function: Compute average value of lambda^{s8} on each sublattice
%     Input: spin_data: Spin configuration in i_site_num \times
%            j_site_num*3 format
%     Output：a i_site_num*j_site_num*3
%     ----------------------------
    i_site_num=size(spin_data,1);
    j_site_num=size(spin_data,2);
    lambda8data=zeros(i_site_num,j_site_num);
    for i=1:i_site_num
        for j=1:j_site_num
            d_z=spin_data(i,j,3);
            lambda8data(i,j)=(1-3*abs(d_z)^2)/sqrt(3);
        end
    end
end