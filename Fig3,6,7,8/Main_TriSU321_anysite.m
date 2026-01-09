%%%%%% S=1 TriSU321 Any-site-sublattice %%%%%
%% Parameter setup
J_3=1
J_2=5
J_1=3
i_site_num=3;
j_site_num=3;
Spin=[]; %Spin storage
E_g=10000; %Ground state energy storage
Iterated=30; %Searching time
options=optimoptions('fmincon','Display','notify-detailed','MaxFunctionEvaluations',1e5,'MaxIterations',1e5);


%% Computing

n=i_site_num*j_site_num;   %Total Site number
Parameter=[J_1,J_2,J_3];
CellSize=[i_site_num,j_site_num];
lb=[zeros(1,n*4),Parameter,CellSize];
ub=[kron(ones(1,n),[1,1,2*pi,2*pi]),Parameter,CellSize];
for i=1:Iterated
    fprintf("Iterated:%d, E_g=%f\n",i,E_g);
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
        Spin=spin;E_g=e_g;
        disp("Accept");
    end
end
spin_data=From_list_to_site(Spin(1:length(Spin)-length(Parameter)-length(CellSize)),CellSize);
c_A=[spin_data(1,1,1);spin_data(1,1,2);spin_data(1,1,3)]
c_B=[spin_data(2,1,1);spin_data(2,1,2);spin_data(2,1,3)]
c_C=[spin_data(1,2,1);spin_data(1,2,2);spin_data(1,2,3)]

%% Following process
lambda8_show;
lambda123_show;
SW_3sub_loop;
Mode1
Mode+Mode1





