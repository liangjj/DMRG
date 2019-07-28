function W=rho_env_eigmat(Psi,n,N) %从env约化密度矩阵中取前N个最大本征值对应的本征矢，Psi列矢中env提供n个基
size_Psi=size(Psi);
m=size_Psi(1)/4/n;
D=zeros(1,2*m);
i=1:2*m;
D(i)=2*n;
RHO=mat2cell(Psi*Psi',D,D);
Rho=cell(1,1,2*m);
for i=1:2*m
    Rho{i}=RHO{i,i};
end
Rho=reshape(cell2mat(Rho),4*n^2,2*m);
rho=zeros(4*n^2,1);
for i=1:4*n^2
    rho(i)=sum(Rho(i,:));
end
rho=reshape(rho,2*n,2*n);
[rho_V,rho_D]=eig(rho);
rho_d=diag(rho_D);
w=zeros(2*n,N);
for i=1:N
    rho_max_position=find(rho_d==max(rho_d));
    w(:,i)=rho_V(:,rho_max_position(1));
    rho_d(rho_max_position(1))=min(rho_d);
end
W=w; %此即描述ne-env这个整体的N个基，以列矢存于此