%%%%% Tri-U(2)-Model Ground-state SSF from Classical variational method %%%%%
%-----------------------------
% Final result will be series of 3-row array. 
%        J1| ? | ? | ...
%        ------------
%        J2| ? | ? | ...
%        ------------
%        one type of SSF| ? | ? | ...
%% Simulation Parameter Setup
Iterated=40  %Searching Time for each parameter points.
%Range of cell size
L_series=[3]    %Expansion factor of unit cell(n_{Lambda})
%Range of parameter
J_3=1
J2_lb=-5    %Lower limit of model parameter J_2
J2_ub=5     %Upper limit of model parameter J_2
J1_lb=0     %Lower limit of model parameter J_1
J1_ub=5     %Upper limit of model parameter J_1
J_step=0.05     %Search step for each parameter
% DataFile Name
FileName="PD_total_data.mat"     %The file name which would storage the calculation results
% SSF extend times
Ext_Time=4      %Times for expanding order unit cell when computing SSF.

%Data Storage and exhibiting setup
SSF123_Gamma=zeros(numel(L_series)+2,int32((abs(J2_ub-J2_lb)/abs(J_step)+1)*(abs(J1_ub-J1_lb)/abs(J_step)+1)));
SSF123_K=zeros(numel(L_series)+2,int32((abs(J2_ub-J2_lb)/abs(J_step)+1)*(abs(J1_ub-J1_lb)/abs(J_step)+1)));
SSF8_Gamma=zeros(numel(L_series)+2,int32((abs(J2_ub-J2_lb)/abs(J_step)+1)*(abs(J1_ub-J1_lb)/abs(J_step)+1)));
SSF8_K=zeros(numel(L_series)+2,int32((abs(J2_ub-J2_lb)/abs(J_step)+1)*(abs(J1_ub-J1_lb)/abs(J_step)+1)));
options=optimoptions('fmincon','Display','notify-detailed','MaxFunctionEvaluations',1e6,'MaxIterations',1e6);

%% Computing
L_index=1;
a_1=[1;0];
a_2=[0.5;sqrt(3)/2];
for L=1:numel(L_series)
    i_site_num=L_series(L);
    j_site_num=L_series(L);
    J_index=1;
    for J_2=J2_lb:J_step:J2_ub
        for J_1=J1_lb:J_step:J1_ub
            SSF123_Gamma(1,J_index)=J_2;
            SSF123_K(1,J_index)=J_2;
            SSF8_Gamma(1,J_index)=J_2;
            SSF8_K(1,J_index)=J_2;
            SSF123_Gamma(2,J_index)=J_1;
            SSF123_K(2,J_index)=J_1;
            SSF8_Gamma(2,J_index)=J_1;
            SSF8_K(2,J_index)=J_1;
            Spin=[];  %Storage the most excellent configuration so far.
            E_g=100000; %Ground State Energy memory
            n=i_site_num*j_site_num;
            Parameter=[J_1,J_2,J_3];
            CellSize=[i_site_num,j_site_num];
            lb=[zeros(1,n*4),Parameter,CellSize];
            ub=[kron(ones(1,n),[1,1,2*pi,2*pi]),Parameter,CellSize];
            for i=1:Iterated
                fprintf("L=%d,J_2=%f,J_1=%f,Times:%d\n",L_series(L),J_2,J_1,i);
                Theta_1=rand(1,n);
                Theta_2=rand(1,n);
                Phi_1=2*pi*rand(1,n);
                Phi_2=2*pi*rand(1,n);
                x_0=zeros(1,n*4);
                x_0(1:4:n*4)=Theta_1;
                x_0(2:4:n*4)=Theta_2;
                x_0(3:4:n*4)=Phi_1;
                x_0(4:4:n*4)=Phi_2;
                X_0=[x_0,Parameter,CellSize];
                [spin,e_g]=fmincon(@(X) Energy(X),X_0,[],[],[],[],lb,ub,[],options);
                if e_g<E_g
                    Spin=spin;
                    E_g=e_g;
                    spin_data=From_list_to_site(Spin(1:length(Spin)-length(Parameter)-length(CellSize)),CellSize);
                    disp("Correlation Function calculating...")
                    lambda123data=From_spin_to_lambda123(spin_data);
                    lambda8data=From_spin_to_lambda8(spin_data);
                    S_123=zeros(i_site_num,j_site_num);
                    S_8=zeros(i_site_num,j_site_num);
                    for i=1:i_site_num
                        for j=1:j_site_num
                            for I=1:i_site_num
                                for J=1:j_site_num
                                    S_123(i,j)=S_123(i,j)+sum(lambda123data(I,J,:).*lambda123data(neighborindex(I,i,i_site_num),neighborindex(J,j,j_site_num),:));
                                    S_8(i,j)=S_8(i,j)+lambda8data(I,J).*lambda8data(neighborindex(I,i,i_site_num),neighborindex(J,j,j_site_num));
                                end
                            end
                        end
                    end
                    for r=1:Ext_Time
                        S_123=[S_123,S_123;S_123,S_123];
                        S_8=[S_8,S_8;S_8,S_8];
                    end
                    Gamma=[0;0];
                    K=[4*pi/3;0];
                    for i=1:size(S_123,1)
                        for j=1:size(S_123,2)
                            R=i*a_1+j*a_2;
                            SSF123_Gamma(L_index+2,J_index)=SSF123_Gamma(L_index+2,J_index)+exp(1i*Gamma'*R)*S_123(i,j);
                            SSF123_K(L_index+2,J_index)=SSF123_K(L_index+2,J_index)+exp(1i*K'*R)*S_123(i,j);
                            SSF8_Gamma(L_index+2,J_index)=SSF8_Gamma(L_index+2,J_index)+exp(1i*Gamma'*R)*S_8(i,j);
                            SSF8_K(L_index+2,J_index)=SSF8_K(L_index+2,J_index)+exp(1i*K'*R)*S_8(i,j);
                        end
                    end
                    SSF123_Gamma(L_index+2,J_index)=real(SSF123_Gamma(L_index+2,J_index))/(size(S_123,1)*size(S_123,2));
                    SSF123_K(L_index+2,J_index)=real(SSF123_K(L_index+2,J_index))/(size(S_123,1)*size(S_123,2));
                    SSF8_Gamma(L_index+2,J_index)=real(SSF8_Gamma(L_index+2,J_index))/(size(S_8,1)*size(S_8,2));
                    SSF8_K(L_index+2,J_index)=real(SSF8_K(L_index+2,J_index))/(size(S_8,1)*size(S_8,2));
                    save(FileName,"SSF123_K","SSF123_Gamma","SSF8_K","SSF8_Gamma");
                    disp("Calcultating finished.");
                end
            end
            J_index=J_index+1;
        end
    end
    L_index=L_index+1;
end




