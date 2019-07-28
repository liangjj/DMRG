function W=rho_sys_eigmat(Psi,m,N) %从sys约化密度矩阵中取前N个最大本征值对应的本征矢，Psi列矢的前m行是以sys为基的系数
size_Psi=size(Psi);
n=size_Psi(1)/4/m;
D=zeros(1,2*m);
i=1:2*m;
D(i)=2*n;
RHO=Psi*Psi';
Rho=reshape(mat2cell(RHO,D,D),[],1);
rho=zeros(4*m^2,1);
for i=1:4*m^2
rho(i)=trace(Rho{i});
end
rho=reshape(rho,2*m,2*m);
[rho_V,rho_D]=eig(rho);
rho_d=diag(rho_D);
w=zeros(2*m,N);
for i=1:N
    rho_max_position=find(rho_d==max(rho_d));
    rho_max=rho_V(:,rho_max_position(1));
    w(:,i)=rho_max;
    rho_d(rho_max_position(1))=min(rho_d);
end
W=w; %此即描述sys-ns这个整体的N个基，以列矢存于此