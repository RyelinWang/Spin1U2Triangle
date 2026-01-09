function E=Energy(X)
    %----------------------------------------------------------------------
    %Function: The Energy functional that will be minized
    %Input: The spin variables list. The latter fourth are Constant set:(...,J_1,J_2,i_site_num,j_site_num)
    %Output: An Energy under a certain configuration described by X.
    %----------------------------------------------------------------------
    l=length(X);
    CellSize=X(l-1:l);
    Parameter=X(l-4:l-2);
    spin=From_list_to_site(X(1:l-5),CellSize);
    %------Calculate current Energy--------
    i_site_num=CellSize(1);
    j_site_num=CellSize(2);
    E=0.0;
    for i=1:i_site_num
        for j=1:j_site_num
            E=E+Energy_per_site(i,j,CellSize,spin,Parameter);
        end
    end
    E=E/(i_site_num*j_site_num);
    if(imag(E)<1e-8)
        E=real(E);
    else
        error("Complex energy output!");
    end
end