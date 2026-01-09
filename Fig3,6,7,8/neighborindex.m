function i_next=neighborindex(i,k,i_site_num)
%-------------------------
%      Function: Output the neighbor index of i with PBC
% 
    i_next=i+k;
    while(i_next>i_site_num)
       i_next=i_next-i_site_num;
    end
    while(i_next<=0)
        i_next=i_next+i_site_num;
    end
    i_next=int32(i_next);
end