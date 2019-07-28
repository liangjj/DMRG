function W=rho_sys_eigmat(Psi,m,N) %��sysԼ���ܶȾ�����ȡǰN�������ֵ��Ӧ�ı���ʸ��Psi��ʸ��ǰm������sysΪ����ϵ��
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
W=w; %�˼�����sys-ns��������N����������ʸ���ڴ�