%% Model parameter setup (De-comment the following six row, one can do the spin-wave at a certain reference state)
% J_3=1
% J_2=5 
% J_1=3
% c_A=[1;0;0];
% c_B=[0;0;1];
% c_C=[0;1;0];
g=0     %Coupling strength with lambda^3.
%% Gell-mann matrix setup
h=0;
lambda=zeros(3,3,8);
lambda_A=zeros(3,3,8);
lambda_B=zeros(3,3,8);
lambda_C=zeros(3,3,8);
lambda(:,:,1)=[0,1,0;1,0,0;0,0,0];
lambda(:,:,2)=[0,-1i,0;1i,0,0;0,0,0];
lambda(:,:,3)=[1,0,0;0,-1,0;0,0,0];
% lambda(:,:,3)=[1,0,0;0,0,0;0,0,-1];
lambda(:,:,4)=[0,0,1;0,0,0;1,0,0];
lambda(:,:,5)=[0,0,-1i;0,0,0;1i,0,0];
lambda(:,:,6)=[0,0,0;0,0,1;0,1,0];
lambda(:,:,7)=[0,0,0;0,0,-1i;0,1i,0];
lambda(:,:,8)=[1,0,0;0,1,0;0,0,-2]/sqrt(3);
% lambda(:,:,8)=[1,0,0;0,-2,0;0,0,1]/sqrt(3);
R_A=eye(3,3);
R_B=[0,1,0;0,0,1;1,0,0];
R_C=[0,1,0;1,0,0;0,0,1];
% R_A=spin1_transform(c_A);
% R_B=spin1_transform(c_B);
% R_C=spin1_transform(c_C);
for m=1:8
    lambda_A(:,:,m)=R_A'*lambda(:,:,m)*R_A;
    lambda_B(:,:,m)=R_B'*lambda(:,:,m)*R_B;
    lambda_C(:,:,m)=R_C'*lambda(:,:,m)*R_C;
end
J=[J_3+J_2,J_3+J_2,J_3+J_2,J_3,J_3,J_3,J_3,J_3+J_1]; %Parameter data
% J=[J_3,J_3,J_3+J_2,J_3+J_2,J_3+J_2,J_3,J_3,J_3+J_1];
%% Hamiltonian cell setup
V_AB=zeros(2,2);
V_BC=zeros(2,2);
V_CA=zeros(2,2);
V_BA=zeros(2,2);
V_CB=zeros(2,2);
V_AC=zeros(2,2);
T_AB=zeros(2,2);
T_BC=zeros(2,2);
T_CA=zeros(2,2);
Delta_AB=zeros(2,2);
Delta_BC=zeros(2,2);
Delta_CA=zeros(2,2);
Z=3;
G_A=lambda_A(2:3,2:3,3)*g;
G_B=lambda_B(2:3,2:3,3)*g;
G_C=lambda_C(2:3,2:3,3)*g;
for m=1:8
    V_AB=V_AB+J(m)*(lambda_A(2:3,2:3,m)-lambda_A(1,1,m)*eye(2,2))*lambda_B(1,1,m);
    V_BA=V_BA+J(m)*(lambda_B(2:3,2:3,m)-lambda_B(1,1,m)*eye(2,2))*lambda_A(1,1,m);
    V_BC=V_BC+J(m)*(lambda_B(2:3,2:3,m)-lambda_B(1,1,m)*eye(2,2))*lambda_C(1,1,m);
    V_CB=V_CB+J(m)*(lambda_C(2:3,2:3,m)-lambda_C(1,1,m)*eye(2,2))*lambda_B(1,1,m);
    V_AC=V_AC+J(m)*(lambda_A(2:3,2:3,m)-lambda_A(1,1,m)*eye(2,2))*lambda_C(1,1,m);
    V_CA=V_CA+J(m)*(lambda_C(2:3,2:3,m)-lambda_C(1,1,m)*eye(2,2))*lambda_A(1,1,m);
    for alpha=1:2
        for beta=1:2
            T_AB(alpha,beta)=T_AB(alpha,beta)+J(m)*lambda_A(alpha+1,1,m)*lambda_B(1,beta+1,m);
            T_CA(alpha,beta)=T_CA(alpha,beta)+J(m)*lambda_C(alpha+1,1,m)*lambda_A(1,beta+1,m);
            T_BC(alpha,beta)=T_BC(alpha,beta)+J(m)*lambda_B(alpha+1,1,m)*lambda_C(1,beta+1,m);
            Delta_AB(alpha,beta)=Delta_AB(alpha,beta)+J(m)*lambda_A(alpha+1,1,m)*lambda_B(beta+1,1,m);
            Delta_BC(alpha,beta)=Delta_BC(alpha,beta)+J(m)*lambda_B(alpha+1,1,m)*lambda_C(beta+1,1,m);
            Delta_CA(alpha,beta)=Delta_CA(alpha,beta)+J(m)*lambda_C(alpha+1,1,m)*lambda_A(beta+1,1,m);
        end
    end
end
delta=[1,1;-2,1;1,-2]'/3;
%% Compute
K_1=[];%Storage k_1 for plot
K_2=[];%Storage k_2 for plot
E_eig_eff=[];%Storage energy for magnons
k_1=0;
k_2=0;
k_step=0.001;
breakpoint=[];
breakpoint(1)=1;
Energy_point_num=0;
while k_1<0.5
    E=SWH_eigen(k_1,k_2,V_AB,V_BC,V_CA,V_BA,V_CB,V_AC,T_AB,T_BC,T_CA,Delta_AB,Delta_BC,Delta_CA,h,lambda_A,lambda_B,lambda_C,delta,Z,G_A,G_B,G_C);
    E_eig_eff=[E_eig_eff,E];
    Energy_point_num=Energy_point_num+1;
    k_1=k_1+k_step;
end
breakpoint=[breakpoint,Energy_point_num+1];
k_1=1/2;
k_2=0;
while k_1<2/3 || k_2<1/3
    E=SWH_eigen(k_1,k_2,V_AB,V_BC,V_CA,V_BA,V_CB,V_AC,T_AB,T_BC,T_CA,Delta_AB,Delta_BC,Delta_CA,h,lambda_A,lambda_B,lambda_C,delta,Z,G_A,G_B,G_C);
    E_eig_eff=[E_eig_eff,E];
    Energy_point_num=Energy_point_num+1;
    k_1=k_1+k_step/sqrt(3);
    k_2=k_2+2*k_step/sqrt(3);
end
breakpoint=[breakpoint,Energy_point_num+1];
k_1=2/3;
k_2=1/3;
while k_1>0 || k_2>0
    E=SWH_eigen(k_1,k_2,V_AB,V_BC,V_CA,V_BA,V_CB,V_AC,T_AB,T_BC,T_CA,Delta_AB,Delta_BC,Delta_CA,h,lambda_A,lambda_B,lambda_C,delta,Z,G_A,G_B,G_C);
    E_eig_eff=[E_eig_eff,E];
    Energy_point_num=Energy_point_num+1;
    k_1=k_1-2*k_step/sqrt(3);
    k_2=k_2-k_step/sqrt(3);
end
breakpoint=[breakpoint,Energy_point_num];


%% Plot
figure
hold on
plot(E_eig_eff(1,:),"LineWidth",2);
plot(E_eig_eff(2,:),"LineWidth",2);
plot(E_eig_eff(3,:),"LineWidth",2);
plot(E_eig_eff(4,:),"LineWidth",2);
plot(E_eig_eff(5,:),"LineWidth",2);
plot(E_eig_eff(6,:),"LineWidth",2);
xticks(breakpoint)
xline(breakpoint,"LineStyle","--","LineWidth",1)
xticklabels(["\Gamma","M","K","\Gamma"])
title(["SU(2)-Supersolid",sprintf("(J_2=%.2f,J_1=%.2f,g=%.2f)",J_2,J_1,g)],"FontSize",20)
xlim([0,1370])
ylim([0,20])
fig=gcf;
ax=findall(fig,"Type","axes");
set(ax,"FontSize",18)
%ylim([-0.5,18])
