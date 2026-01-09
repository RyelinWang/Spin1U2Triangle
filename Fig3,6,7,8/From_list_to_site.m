function spin=From_list_to_site(X,CellSize)
    %-----------------------------------------------------------------------------------------
    %Function: To change the list of site-variables into the mode of site,
    %within the form (d_+,d_0,d_-)
    %Input: -1) X:  The array, which storage the spin vector in (theta_1,theta_2,phi_1,phi_2) and are ready to be minimize.
    %       -2) i_site_num,j_site_num:  Test order length on x,y direction
    %Output: An array, the meaning of each components defined as spin(i-label,j-label,1/0/-1-label)
    %Notes!: psi=d_+|1>+d_0|0>+d_-|-1>
    %-----------------------------------------------------------------------------------------
    n=length(X);
    i_site_num=int32(CellSize(1));
    j_site_num=int32(CellSize(2));
    spin=zeros(i_site_num,j_site_num,3);
    if((i_site_num*j_site_num)~=n/4)
        % Check whether the number of X isn't consistent to the set of.
        error("Wrong! The number of variables is inconsistent!");
    end
    
    List_label=1;  %A step counter which is used to count the label of array X
    for i=1:i_site_num
        for j=1:j_site_num 
            s_tmp=S(X(List_label),X(List_label+1),X(List_label+2),X(List_label+3));
            spin(i,j,1)=s_tmp(1);
            spin(i,j,2)=s_tmp(2);
            spin(i,j,3)=s_tmp(3);
            List_label=List_label+4;
        end
    end
    if(List_label-1~=n)
        error("The distribution of variables got wrong!");
    end
end