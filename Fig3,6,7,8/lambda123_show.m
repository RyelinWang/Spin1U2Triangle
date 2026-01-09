L=1; %Size of configuration show 
a=2; 
a_1=[1;0]*a;
a_2=[1/2;sqrt(3)/2]*a;
%Location in one order cell
X_A_incell=zeros(i_site_num,j_site_num);
Y_A_incell=zeros(i_site_num,j_site_num);
for i=1:i_site_num
    for j=1:j_site_num
        X_A_incell(i,j)=(i-1)*a_1(1)+(j-1)*a_2(1);
        Y_A_incell(i,j)=(i-1)*a_1(2)+(j-1)*a_2(2);
    end
end
%Location for whole sample
X_A=X_A_incell;
Y_A=Y_A_incell;
for l=1:L-1
    X_A=[X_A,X_A+j_site_num*a_2(1);X_A+i_site_num*a_1(1),X_A+i_site_num*a_1(1)+j_site_num*a_2(1)];
    Y_A=[Y_A,Y_A+j_site_num*a_2(2);Y_A+i_site_num*a_1(2),Y_A+i_site_num*a_1(2)+j_site_num*a_2(2)];
end
Z_A=zeros(size(X_A));

axis equal
lambda123Data=From_spin_to_lambda123(spin_data);
lambda_1=lambda123Data(:,:,1);
lambda_2=lambda123Data(:,:,2);
lambda_3=lambda123Data(:,:,3);
for l=1:L-1
    lambda_1=[lambda_1,lambda_1;lambda_1,lambda_1];
    lambda_2=[lambda_2,lambda_2;lambda_2,lambda_2];
    lambda_3=[lambda_3,lambda_3;lambda_3,lambda_3];
end


factor=0.7;
startX=X_A(:)-0.5*lambda_1(:);
startY=Y_A(:)-0.5*lambda_2(:);
startZ=Z_A(:)-0.5*lambda_3(:);
figure
quiver3(startX,startY,startZ,lambda_1(:),lambda_2(:),lambda_3(:),0,'Color', 'b',"Marker",'.','MarkerSize',0.5,'LineWidth',1.5);
axis off
title(sprintf("\\lambda_{123} configuration, J_2=%.1f,J_1=%.1f ",J_2,J_1),"FontSize",10);

hold on
% 获取所有子格点的坐标（A、B、C）
pointsX = [X_A(:)];
pointsY = [Y_A(:)];

% 计算所有点之间的距离矩阵
distMatrix = squareform(pdist([pointsX, pointsY]));
threshold = a+0.1; % 根据晶格常数设置阈值，a=2时连接距离为1的边

% 绘制所有相邻边
[i, j] = find(distMatrix < threshold & distMatrix > 0);
for k = 1:length(i)
    if i(k) < j(k) % 避免重复绘制
        plot([pointsX(i(k)), pointsX(j(k))], [pointsY(i(k)), pointsY(j(k))], ...
             'k-', 'LineWidth', 1,'LineStyle','--');
    end
end
axis equal;
hold off;
