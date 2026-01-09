function E=SWH_eigen(k_1,k_2,V_AB,V_BC,V_CA,V_BA,V_CB,V_AC,T_AB,T_BC,T_CA,Delta_AB,Delta_BC,Delta_CA,h,lambda_A,lambda_B,lambda_C,delta,Z,G_A,G_B,G_C)
    k=[k_1,k_2];
    gamma=0.0; 
    for nn=1:3
        gamma=gamma+exp(-1i*2*pi*k*delta(:,nn));
    end
    %Hamiltonian construct
    H_1=[Z*V_AB+Z*V_AC+2*G_A,T_AB*gamma,T_CA'*conj(gamma);T_AB'*conj(gamma),Z*V_BC+Z*V_BA+2*G_B,T_BC*gamma;T_CA*gamma,T_BC'*conj(gamma),Z*V_CA+Z*V_CB+2*G_C];
    H_2=[Z*V_AB.'+Z*V_AC.'+2*G_A.',conj(T_AB)*gamma,T_CA.'*conj(gamma);T_AB.'*conj(gamma),Z*V_BC.'+Z*V_BA.'+2*G_B.',conj(T_BC)*gamma;conj(T_CA)*gamma,T_BC.'*conj(gamma),Z*V_CA.'+Z*V_CB.'+2*G_C.'];
    Delta=[zeros(2,2),Delta_AB*gamma,Delta_CA.'*conj(gamma);Delta_AB.'*conj(gamma),zeros(2,2),Delta_BC*gamma;Delta_CA*gamma,Delta_BC.'*conj(gamma),zeros(2,2)];
%        Delta=[zeros(2,2),2*Delta_AB*gamma,2*Delta_CA*gamma;zeros(2,2),zeros(2,2),2*Delta_BC*gamma;zeros(2,2),zeros(2,2),zeros(2,2)];
    H_1=H_1+h*[lambda_A(2:3,2:3,8),zeros(2,2),zeros(2,2);zeros(2,2),lambda_B(2:3,2:3,8),zeros(2,2);zeros(2,2),zeros(2,2),lambda_C(2:3,2:3,8)];
    H_2=H_2+h*[lambda_A(2:3,2:3,8),zeros(2,2),zeros(2,2);zeros(2,2),lambda_B(2:3,2:3,8),zeros(2,2);zeros(2,2),zeros(2,2),lambda_C(2:3,2:3,8)];
    H=[H_1,Delta;Delta',H_2]/2;
    Sigma=[eye(6,6),zeros(6,6);zeros(6,6),-eye(6,6)];
    [Vec,EigenAll]=eig(Sigma*H);
    E_eigen=[];
    for m=1:12
        if Vec(:,m)'*Sigma*Vec(:,m)>0
            E_eigen=[E_eigen;EigenAll(m,m)];
        end
    end
    if (sum(abs(imag(E_eigen)))>sum(abs(real(E_eigen)))*1e-3)
            fprintf("k=(%.3f,%.3f), energy instable\n",k_1,k_2);
            E=nan*ones(6,1);
    else if(numel(E_eigen)~=6)
            fprintf("k=(%.2f,%.2f),No solution\n",k_1,k_2);
            E=nan*ones(6,1);
         else
            E_eigen=sort(real(E_eigen));
            %E(i,j,:)=E_eigen(numel(E_eigen)/2+1:numel(E_eigen));
            %E(i,j,:)=E_eigen(1:numel(E_eigen)/2);
            E=E_eigen;
         end
    end
end