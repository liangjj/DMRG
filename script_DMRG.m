tic;
clear;
save DMRG.mat -v7.3;
L=12; %����4���ܸ������ż����
N=20; %sweepѭ������
dmrg=matfile('DMRG.mat','Writable',true);
dmrg.H_sys=cast(zeros(2),'like',1i); %���м�����5�У���ΪDMRG.mat�ļ��е�6�ľ��������H_sys,H_sys_ns�ȣ���ʼ����һ��ȷ��ÿ���������������Ϊ��������֤�����ھ��������д�븴������������������Ҫʹ�ó�ʼ���ĸ�ֵ��
dmrg.H_sys_ns=cast(Hamiltonian(2),'like',1i);
dmrg.H_ne_env=cast(Hamiltonian(2),'like',1i);
dmrg.H_env=cast(zeros(2),'like',1i);
dmrg.rho_sys_W=cast(zeros(4,2),'like',1i);
dmrg.rho_env_W=cast(zeros(4,2),'like',1i);
H=Hamiltonian(4);
for i=1:L/2-2 %�ԳƳ���
    [H_V,H_D]=eig(H);
    H_d=real(diag(H_D));
    H_min_position=find(H_d==min(H_d));
    Psi_0=H_V(:,H_min_position(1));
    rho_sys_W=rho_sys_eigmat(Psi_0,2,2); %�˴�����rho_sys_eigmat��������������sys��size=i+1����2����
    dmrg.rho_sys_W(4*i+1:4*i+4,:)=rho_sys_W;
    rho_env_W=rho_env_eigmat(Psi_0,2,2); %�˴�����rho_env_eigmat��������������env��size=i+1����2����
    dmrg.rho_env_W(4*i+1:4*i+4,:)=rho_env_W;
    H_sys=rho_sys_W'*(kron(dmrg.H_sys(2*i-1:2*i,:),eye(2))+dmrg.H_sys_ns(4*i-3:4*i,:))*rho_sys_W; %��i����H_sys��H_sys_ns����i+1����H_sys
    dmrg.H_sys(2*i+1:2*i+2,:)=H_sys;
    H_sys_ns=Hamiltonian_sys_ns(rho_sys_W); %����Hamiltonian_sys_ns��������sys���ܶȾ�����ʸrho_sys_W����i+1����H_sys_ns
    dmrg.H_sys_ns(4*i+1:4*i+4,:)=H_sys_ns;
    H_ne_env=Hamiltonian_ne_env(rho_env_W); %����Hamiltonian_ne_env��������env���ܶȾ�����ʸrho_env_W����i+1����H_ne_env
    dmrg.H_ne_env(4*i+1:4*i+4,:)=H_ne_env;
    H_env=rho_env_W'*(kron(eye(2),dmrg.H_env(2*i-1:2*i,:))+dmrg.H_ne_env(4*i-3:4*i,:))*rho_env_W; %��i����H_env��H_ne_env����i+1����H_env
    dmrg.H_env(2*i+1:2*i+2,:)=H_env;
    H=kron(H_sys,eye(8))+kron(H_sys_ns,eye(4))+kron(eye(2),kron(Hamiltonian(2),eye(2)))+kron(eye(4),H_ne_env)+kron(eye(8),H_env); %16*16��size(sys)=i+1��Hamiltonian
end
for j=1:N %sweep
    for i=L/2-1:L-4 %���м�����sweep
        [H_V,H_D]=eig(H);
        H_d=real(diag(H_D));
        H_min_position=find(H_d==min(H_d));
        Psi_0=H_V(:,H_min_position(1));
        rho_sys_W=rho_sys_eigmat(Psi_0,2,2); %�˴�����rho_sys_eigmat��������������sys��size=i+1����2����
        dmrg.rho_sys_W(4*i+1:4*i+4,:)=rho_sys_W;
        H_sys=rho_sys_W'*(kron(dmrg.H_sys(2*i-1:2*i,:),eye(2))+dmrg.H_sys_ns(4*i-3:4*i,:))*rho_sys_W; %��i����H_sys��H_sys_ns����i+1����H_sys
        dmrg.H_sys(2*i+1:2*i+2,:)=H_sys;
        H_sys_ns=Hamiltonian_sys_ns(rho_sys_W); %����Hamiltonian_sys_ns��������sys���ܶȾ�����ʸrho_sys_W����i+1����H_sys_ns
        dmrg.H_sys_ns(4*i+1:4*i+4,:)=H_sys_ns;
        H_ne_env=dmrg.H_ne_env(4*(L-i-3)-3:4*(L-i-3),:); %��dmrg.mat����size=L-i-3��H_ne_env
        H_env=dmrg.H_env(2*(L-i-3)-1:2*(L-i-3),:); %��dmrg.mat����size=L-i-3��H_env
        H=kron(H_sys,eye(8))+kron(H_sys_ns,eye(4))+kron(eye(2),kron(Hamiltonian(2),eye(2)))+kron(eye(4),H_ne_env)+kron(eye(8),H_env); %16*16��size(sys)=i+1��size(env)=L-i-3��Hamiltonian
    end
    for i=1:L-4 %��������sweep
        [H_V,H_D]=eig(H);
        H_d=real(diag(H_D));
        H_min_position=find(H_d==min(H_d));
        Psi_0=H_V(:,H_min_position(1));
        rho_env_W=rho_env_eigmat(Psi_0,2,2); %�˴�����rho_env_eigmat��������������env��size=i+1����2����
        dmrg.rho_env_W(4*i+1:4*i+4,:)=rho_env_W;
        H_env=rho_env_W'*(kron(eye(2),dmrg.H_env(2*i-1:2*i,:))+dmrg.H_ne_env(4*i-3:4*i,:))*rho_env_W; %��i����H_env��H_ne_env����i+1����H_env
        dmrg.H_env(2*i+1:2*i+2,:)=H_env;
        H_ne_env=Hamiltonian_ne_env(rho_env_W); %����Hamiltonian_ne_env��������env���ܶȾ�����ʸrho_env_W����i+1����H_ne_env
        dmrg.H_ne_env(4*i+1:4*i+4,:)=H_ne_env;
        H_sys=dmrg.H_sys(2*(L-i-3)-1:2*(L-i-3),:); %��dmrg.mat����size=L-i-3��H_sys
        H_sys_ns=dmrg.H_sys_ns(4*(L-i-3)-3:4*(L-i-3),:); %��dmrg.mat����size=L-i-3��H_sys_ns
        H=kron(H_sys,eye(8))+kron(H_sys_ns,eye(4))+kron(eye(2),kron(Hamiltonian(2),eye(2)))+kron(eye(4),H_ne_env)+kron(eye(8),H_env); %16*16��size(env)=i+1��size(sys)=L-i-3��Hamiltonian
    end
    for i=1:L/2-2 %�������м�sweep
        [H_V,H_D]=eig(H);
        H_d=real(diag(H_D));
        H_min_position=find(H_d==min(H_d));
        Psi_0=H_V(:,H_min_position(1));
        rho_sys_W=rho_sys_eigmat(Psi_0,2,2); %�˴�����rho_sys_eigmat��������������sys��size=i+1����2����
        dmrg.rho_sys_W(4*i+1:4*i+4,:)=rho_sys_W;
        H_sys=rho_sys_W'*(kron(dmrg.H_sys(2*i-1:2*i,:),eye(2))+dmrg.H_sys_ns(4*i-3:4*i,:))*rho_sys_W; %��i����H_sys��H_sys_ns����i+1����H_sys
        dmrg.H_sys(2*i+1:2*i+2,:)=H_sys;
        H_sys_ns=Hamiltonian_sys_ns(rho_sys_W); %����Hamiltonian_sys_ns��������sys���ܶȾ�����ʸrho_sys_W����i+1����H_sys_ns
        dmrg.H_sys_ns(4*i+1:4*i+4,:)=H_sys_ns;
        H_ne_env=dmrg.H_ne_env(4*(L-i-3)-3:4*(L-i-3),:); %��dmrg.mat����size=L-i-3��H_ne_env
        H_env=dmrg.H_env(2*(L-i-3)-1:2*(L-i-3),:); %��dmrg.mat����size=L-i-3��H_env
        H=kron(H_sys,eye(8))+kron(H_sys_ns,eye(4))+kron(eye(2),kron(Hamiltonian(2),eye(2)))+kron(eye(4),H_ne_env)+kron(eye(8),H_env); %16*16��size(sys)=i+1��size(env)=L-i-3��Hamiltonian
    end
end
[H_V,H_D]=eig(H);
H_d=real(diag(H_D));
E_0=min(H_d);
dmrg.E_0=E_0;
H_min_position=find(H_d==E_0);
n=size(H_min_position);
dmrg.Psi_0=cast(zeros(16*n(1),1),'like',1i);
for i=1:n(1)
    dmrg.Psi_0((i-1)*16+1:i*16,1)=H_V(:,H_min_position(i));
end
fprintf('E_0=%f\n',E_0);
if n(1)==1
    fprintf('No degeneracy in groud-state.\n');
else
    fprintf('%d-fold degeneracy in groud-state.\n',n(1));
end
toc;