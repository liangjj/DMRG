tic;
clear;
save DMRG.mat -v7.3;
dmrg=matfile('DMRG.mat','Writable',true);
L=100; %����4���ܸ������ż����
M=10; %sys��env��ౣ��2<=M<2^(L/2-2)��̬
N=30; %sweepѭ������
m=floor(log(M)/log(2));
dmrg.H_sys=cast(zeros(2,M),'like',1i); %���м�����5�У���ΪDMRG.mat�ļ��е�6�����������H_sys,H_sys_ns�ȣ���ʼ����һ��ȷ��ÿ���������������Ϊ��������֤�����ھ��������д�븴������������������Ҫʹ�ó�ʼ���ĸ�ֵ��
dmrg.H_sys_ns=cast([Hamiltonian(2),zeros(4,2*M-4)],'like',1i);
dmrg.H_ne_env=cast([Hamiltonian(2),zeros(4,2*M-4)],'like',1i);
dmrg.H_env=cast(zeros(2,M),'like',1i);
dmrg.rho_sys_W=cast(zeros(2,M),'like',1i);
dmrg.rho_env_W=cast(zeros(2,M),'like',1i);
H=Hamiltonian(4);
for i=1:m-1 %�ԳƳ��㣬ֱ��size(sys)=m
    [H_V,H_D]=eig(H);
    H_d=real(diag(H_D));
    H_min_position=find(H_d==min(H_d));
    Psi_0=H_V(:,H_min_position(1));
    rho_sys_W=rho_sys_eigmat(Psi_0,2^i,2^(i+1)); %����rho_sys_eigmat��������������sys��size=i+1����2^(i+1)����
    dmrg.rho_sys_W(2^(i+1)-1:2^(i+2)-2,:)=[rho_sys_W,zeros(2^(i+1),M-2^(i+1))];
    rho_env_W=rho_env_eigmat(Psi_0,2^i,2^(i+1)); %����rho_env_eigmat��������������env��size=i+1����2^(i+1)����
    dmrg.rho_env_W(2^(i+1)-1:2^(i+2)-2,:)=[rho_env_W,zeros(2^(i+1),M-2^(i+1))];
    H_sys=rho_sys_W'*(kron(dmrg.H_sys(2^i-1:2^(i+1)-2,1:2^i),eye(2))+dmrg.H_sys_ns(2^(i+1)-3:2^(i+2)-4,1:2^(i+1)))*rho_sys_W; %��i����H_sys��H_sys_ns����i+1����H_sys
    dmrg.H_sys(2^(i+1)-1:2^(i+2)-2,:)=[H_sys,zeros(2^(i+1),M-2^(i+1))];
    H_sys_ns=Hamiltonian_sys_ns(rho_sys_W); %����Hamiltonian_sys_ns��������sys���ܶȾ�����ʸrho_sys_W����i+1����H_sys_ns
    dmrg.H_sys_ns(2^(i+2)-3:2^(i+3)-4,:)=[H_sys_ns,zeros(2^(i+2),2*M-2^(i+2))];
    H_ne_env=Hamiltonian_ne_env(rho_env_W); %����Hamiltonian_ne_env��������env���ܶȾ�����ʸrho_env_W����i+1����H_ne_env
    dmrg.H_ne_env(2^(i+2)-3:2^(i+3)-4,:)=[H_ne_env,zeros(2^(i+2),2*M-2^(i+2))];
    H_env=rho_env_W'*(kron(eye(2),dmrg.H_env(2^i-1:2^(i+1)-2,1:2^i))+dmrg.H_ne_env(2^(i+1)-3:2^(i+2)-4,1:2^(i+1)))*rho_env_W; %��i����H_env��H_ne_env����i+1����H_env
    dmrg.H_env(2^(i+1)-1:2^(i+2)-2,:)=[H_env,zeros(2^(i+1),M-2^(i+1))];
    H=kron(H_sys,eye(2^(i+3)))+kron(H_sys_ns,eye(2^(i+2)))+kron(kron(eye(2^(i+1)),Hamiltonian(2)),eye(2^(i+1)))+kron(eye(2^(i+2)),H_ne_env)+kron(eye(2^(i+3)),H_env); %size(sys)=size(env)=i+1��Hamiltonian
end
[H_V,H_D]=eig(H); %���м�����rho_..._W���⣬�ԳƳ�1���㣬�����¾���ͬ
H_d=real(diag(H_D));
H_min_position=find(H_d==min(H_d));
Psi_0=H_V(:,H_min_position(1));
rho_sys_W=rho_sys_eigmat(Psi_0,2^m,M); %����rho_sys_eigmat��������������sys��size=m+1����M����
dmrg.rho_sys_W(2^(m+1)-1:2^(m+2)-2,:)=rho_sys_W;
rho_env_W=rho_env_eigmat(Psi_0,2^m,M); %����rho_env_eigmat��������������env��size=m+1����M����
dmrg.rho_env_W(2^(m+1)-1:2^(m+2)-2,:)=rho_env_W;
H_sys=rho_sys_W'*(kron(dmrg.H_sys(2^m-1:2^(m+1)-2,1:2^m),eye(2))+dmrg.H_sys_ns(2^(m+1)-3:2^(m+2)-4,1:2^(m+1)))*rho_sys_W; %��m����H_sys��H_sys_ns����m+1����H_sys
dmrg.H_sys(2^(m+1)-1:2^(m+1)+M-2,:)=H_sys;
H_sys_ns=Hamiltonian_sys_ns(rho_sys_W); %����Hamiltonian_sys_ns��������sys���ܶȾ�����ʸrho_sys_W����m+1����H_sys_ns
dmrg.H_sys_ns(2^(m+2)-3:2^(m+2)-4+2*M,:)=H_sys_ns;
H_ne_env=Hamiltonian_ne_env(rho_env_W); %����Hamiltonian_ne_env��������env���ܶȾ�����ʸrho_env_W����m+1����H_ne_env
dmrg.H_ne_env(2^(m+2)-3:2^(m+2)-4+2*M,:)=H_ne_env;
H_env=rho_env_W'*(kron(eye(2),dmrg.H_env(2^m-1:2^(m+1)-2,1:2^m))+dmrg.H_ne_env(2^(m+1)-3:2^(m+2)-4,1:2^(m+1)))*rho_env_W; %��m����H_env��H_ne_env����m+1����H_env
dmrg.H_env(2^(m+1)-1:2^(m+1)+M-2,:)=H_env;
H=kron(H_sys,eye(4*M))+kron(H_sys_ns,eye(2*M))+kron(kron(eye(M),Hamiltonian(2)),eye(M))+kron(eye(2*M),H_ne_env)+kron(eye(4*M),H_env); %size(sys)=size(env)=m+1��Hamiltonian  %��������rho_..._W���⣬�ԳƳ�1���㣬�����¾���ͬ
for i=m+1:L/2-2 %�ԳƳ���
    [H_V,H_D]=eig(H);
    H_d=real(diag(H_D));
    H_min_position=find(H_d==min(H_d));
    Psi_0=H_V(:,H_min_position(1));
    rho_sys_W=rho_sys_eigmat(Psi_0,M,M); %����rho_sys_eigmat��������������sys��size=i+1����M����
    dmrg.rho_sys_W(2^(m+2)-2+2*M*(i-m-1)+1:2^(m+2)-2+2*M*(i-m),:)=rho_sys_W;
    rho_env_W=rho_env_eigmat(Psi_0,M,M); %����rho_env_eigmat��������������env��size=i+1����M����
    dmrg.rho_env_W(2^(m+2)-2+2*M*(i-m-1)+1:2^(m+2)-2+2*M*(i-m),:)=rho_env_W;
    H_sys=rho_sys_W'*(kron(dmrg.H_sys(2^(m+1)+(i-m-1)*M-1:2^(m+1)+(i-m)*M-2,:),eye(2))+dmrg.H_sys_ns(2^(m+2)-3+(i-m-1)*2*M:2^(m+2)-4+(i-m)*2*M,:))*rho_sys_W; %��i����H_sys��H_sys_ns����i+1����H_sys
    dmrg.H_sys(2^(m+1)-1+(i-m)*M:2^(m+1)-2+(i-m+1)*M,:)=H_sys;
    H_sys_ns=Hamiltonian_sys_ns(rho_sys_W); %����Hamiltonian_sys_ns��������sys���ܶȾ�����ʸrho_sys_W����i+1����H_sys_ns
    dmrg.H_sys_ns(2^(m+2)-3+(i-m)*2*M:2^(m+2)-4+(i-m+1)*2*M,:)=H_sys_ns;
    H_ne_env=Hamiltonian_ne_env(rho_env_W); %����Hamiltonian_ne_env��������env���ܶȾ�����ʸrho_env_W����i+1����H_ne_env
    dmrg.H_ne_env(2^(m+2)-3+(i-m)*2*M:2^(m+2)-4+(i-m+1)*2*M,:)=H_ne_env;
    H_env=rho_env_W'*(kron(eye(2),dmrg.H_env(2^(m+1)+(i-m-1)*M-1:2^(m+1)+(i-m)*M-2,:))+dmrg.H_ne_env(2^(m+2)-3+(i-m-1)*2*M:2^(m+2)-4+(i-m)*2*M,:))*rho_env_W; %��i����H_env��H_ne_env����i+1����H_env
    dmrg.H_env(2^(m+1)-1+(i-m)*M:2^(m+1)-2+(i-m+1)*M,:)=H_env;
    H=kron(H_sys,eye(4*M))+kron(H_sys_ns,eye(2*M))+kron(kron(eye(M),Hamiltonian(2)),eye(M))+kron(eye(2*M),H_ne_env)+kron(eye(4*M),H_env); %size(sys)=size(env)=i+1��Hamiltonian
end
for j=1:N
    for i=L/2-1:L-m-4 %���м�����sweep��size(sys)=L-m-3
        [H_V,H_D]=eig(H);
        H_d=real(diag(H_D));
        H_min_position=find(H_d==min(H_d));
        Psi_0=H_V(:,H_min_position(1));
        rho_sys_W=rho_sys_eigmat(Psi_0,M,M); %����rho_sys_eigmat��������������sys��size=i+1����M����
        dmrg.rho_sys_W(2^(m+2)-2+2*M*(i-m-1)+1:2^(m+2)-2+2*M*(i-m),:)=rho_sys_W;
        H_sys=rho_sys_W'*(kron(dmrg.H_sys(2^(m+1)+(i-m-1)*M-1:2^(m+1)+(i-m)*M-2,:),eye(2))+dmrg.H_sys_ns(2^(m+2)-3+(i-m-1)*2*M:2^(m+2)-4+(i-m)*2*M,:))*rho_sys_W; %��i����H_sys��H_sys_ns����i+1����H_sys
        dmrg.H_sys(2^(m+1)-1+(i-m)*M:2^(m+1)-2+(i-m+1)*M,:)=H_sys;
        H_sys_ns=Hamiltonian_sys_ns(rho_sys_W); %����Hamiltonian_sys_ns��������sys���ܶȾ�����ʸrho_sys_W����i+1����H_sys_ns
        dmrg.H_sys_ns(2^(m+2)-3+(i-m)*2*M:2^(m+2)-4+(i-m+1)*2*M,:)=H_sys_ns;
        H_ne_env=dmrg.H_ne_env(2^(m+2)-3+(L-i-4-m)*2*M:2^(m+2)-4+(L-i-4-m+1)*2*M,:);
        H_env=dmrg.H_env(2^(m+1)-1+(L-i-4-m)*M:2^(m+1)-2+(L-i-4-m+1)*M,:);
        H=kron(H_sys,eye(4*M))+kron(H_sys_ns,eye(2*M))+kron(kron(eye(M),Hamiltonian(2)),eye(M))+kron(eye(2*M),H_ne_env)+kron(eye(4*M),H_env); %size(sys)=i+1��size(env)=L-i-3��Hamiltonian
    end
    for i=L-m-3:L-4 %��������sweep
        [H_V,H_D]=eig(H);
        H_d=real(diag(H_D));
        H_min_position=find(H_d==min(H_d));
        Psi_0=H_V(:,H_min_position(1));
        rho_sys_W=rho_sys_eigmat(Psi_0,M,M); %����rho_sys_eigmat��������������sys��size=i+1����M����
        dmrg.rho_sys_W(2^(m+2)-2+2*M*(i-m-1)+1:2^(m+2)-2+2*M*(i-m),:)=rho_sys_W;
        H_sys=rho_sys_W'*(kron(dmrg.H_sys(2^(m+1)+(i-m-1)*M-1:2^(m+1)+(i-m)*M-2,:),eye(2))+dmrg.H_sys_ns(2^(m+2)-3+(i-m-1)*2*M:2^(m+2)-4+(i-m)*2*M,:))*rho_sys_W; %��i����H_sys��H_sys_ns����i+1����H_sys
        dmrg.H_sys(2^(m+1)-1+(i-m)*M:2^(m+1)-2+(i-m+1)*M,:)=H_sys;
        H_sys_ns=Hamiltonian_sys_ns(rho_sys_W); %����Hamiltonian_sys_ns��������sys���ܶȾ�����ʸrho_sys_W����i+1����H_sys_ns
        dmrg.H_sys_ns(2^(m+2)-3+(i-m)*2*M:2^(m+2)-4+(i-m+1)*2*M,:)=H_sys_ns;
        H_ne_env=dmrg.H_ne_env(2^(L-i-4+2)-3:2^(L-i-4+3)-4,1:2^(L-i-2));
        H_env=dmrg.H_env(2^(L-i-4+1)-1:2^(L-i-4+2)-2,1:2^(L-i-3));
        H=kron(H_sys,eye(2^(L-i-1)))+kron(H_sys_ns,eye(2^(L-i-2)))+kron(kron(eye(M),Hamiltonian(2)),eye(2^(L-i-3)))+kron(eye(M*2),H_ne_env)+kron(eye(M*4),H_env); %size(sys)=i+1��size(env)=L-i-3��Hamiltonian
    end
    for i=1:m-1 %��������sweep��size(env)=m
        [H_V,H_D]=eig(H);
        H_d=real(diag(H_D));
        H_min_position=find(H_d==min(H_d));
        Psi_0=H_V(:,H_min_position(1));
        rho_env_W=rho_env_eigmat(Psi_0,2^i,2^(i+1)); %����rho_env_eigmat��������������env��size=i+1����2^(i+1)����
        dmrg.rho_env_W(2^(i+1)-1:2^(i+2)-2,:)=[rho_env_W,zeros(2^(i+1),M-2^(i+1))];
        H_ne_env=Hamiltonian_ne_env(rho_env_W); %����Hamiltonian_ne_env��������env���ܶȾ�����ʸrho_env_W����i+1����H_ne_env
        dmrg.H_ne_env(2^(i+2)-3:2^(i+3)-4,:)=[H_ne_env,zeros(2^(i+2),2*M-2^(i+2))];
        H_env=rho_env_W'*(kron(eye(2),dmrg.H_env(2^i-1:2^(i+1)-2,1:2^i))+dmrg.H_ne_env(2^(i+1)-3:2^(i+2)-4,1:2^(i+1)))*rho_env_W; %��i����H_env��H_ne_env����i+1����H_env
        dmrg.H_env(2^(i+1)-1:2^(i+2)-2,:)=[H_env,zeros(2^(i+1),M-2^(i+1))];
        H_sys=dmrg.H_sys(2^(m+1)-1+(L-i-4-m)*M:2^(m+1)-2+(L-i-4-m+1)*M,:);
        H_sys_ns=dmrg.H_sys_ns(2^(m+2)-3+(L-i-4-m)*2*M:2^(m+2)-4+(L-i-4-m+1)*2*M,:);
        H=kron(H_sys,eye(2^(i+3)))+kron(H_sys_ns,eye(2^(i+2)))+kron(kron(eye(M),Hamiltonian(2)),eye(2^(i+1)))+kron(eye(2*M),H_ne_env)+kron(eye(4*M),H_env); %size(sys)=L-i-3��size(env)=i+1��Hamiltonian
    end
    [H_V,H_D]=eig(H);
    H_d=real(diag(H_D));
    H_min_position=find(H_d==min(H_d));
    Psi_0=H_V(:,H_min_position(1));
    rho_env_W=rho_env_eigmat(Psi_0,2^m,M); %����rho_env_eigmat��������������env��size=m+1����M����
    dmrg.rho_env_W(2^(m+1)-1:2^(m+2)-2,:)=rho_env_W;
    H_ne_env=Hamiltonian_ne_env(rho_env_W); %����Hamiltonian_ne_env��������env���ܶȾ�����ʸrho_env_W����m+1����H_ne_env
    dmrg.H_ne_env(2^(m+2)-3:2^(m+2)-4+2*M,:)=H_ne_env;
    H_env=rho_env_W'*(kron(eye(2),dmrg.H_env(2^m-1:2^(m+1)-2,1:2^m))+dmrg.H_ne_env(2^(m+1)-3:2^(m+2)-4,1:2^(m+1)))*rho_env_W; %��m����H_env��H_ne_env����m+1����H_env
    dmrg.H_env(2^(m+1)-1:2^(m+1)+M-2,:)=H_env;
    H_sys=dmrg.H_sys(2^(m+1)-1+(L-m-4-m)*M:2^(m+1)-2+(L-m-4-m+1)*M,:);
    H_sys_ns=dmrg.H_sys_ns(2^(m+2)-3+(L-m-4-m)*2*M:2^(m+2)-4+(L-m-4-m+1)*2*M,:);
    H=kron(H_sys,eye(4*M))+kron(H_sys_ns,eye(2*M))+kron(kron(eye(M),Hamiltonian(2)),eye(M))+kron(eye(2*M),H_ne_env)+kron(eye(4*M),H_env); %size(sys)=L-i-3��size(env)=m+1��Hamiltonian
    for i=m+1:L-m-4
        [H_V,H_D]=eig(H);
        H_d=real(diag(H_D));
        H_min_position=find(H_d==min(H_d));
        Psi_0=H_V(:,H_min_position(1));
        rho_env_W=rho_env_eigmat(Psi_0,M,M); %����rho_env_eigmat��������������env��size=i+1����M����
        dmrg.rho_env_W(2^(m+2)-2+2*M*(i-m-1)+1:2^(m+2)-2+2*M*(i-m),:)=rho_env_W;
        H_ne_env=Hamiltonian_ne_env(rho_env_W); %����Hamiltonian_ne_env��������env���ܶȾ�����ʸrho_env_W����i+1����H_ne_env
        dmrg.H_ne_env(2^(m+2)-3+(i-m)*2*M:2^(m+2)-4+(i-m+1)*2*M,:)=H_ne_env;
        H_env=rho_env_W'*(kron(eye(2),dmrg.H_env(2^(m+1)+(i-m-1)*M-1:2^(m+1)+(i-m)*M-2,:))+dmrg.H_ne_env(2^(m+2)-3+(i-m-1)*2*M:2^(m+2)-4+(i-m)*2*M,:))*rho_env_W; %��i����H_env��H_ne_env����i+1����H_env
        dmrg.H_env(2^(m+1)-1+(i-m)*M:2^(m+1)-2+(i-m+1)*M,:)=H_env;
        H_sys=dmrg.H_sys(2^(m+1)-1+(L-i-4-m)*M:2^(m+1)-2+(L-i-4-m+1)*M,:);
        H_sys_ns=dmrg.H_sys_ns(2^(m+2)-3+(L-i-4-m)*2*M:2^(m+2)-4+(L-i-4-m+1)*2*M,:);
        H=kron(H_sys,eye(4*M))+kron(H_sys_ns,eye(2*M))+kron(kron(eye(M),Hamiltonian(2)),eye(M))+kron(eye(2*M),H_ne_env)+kron(eye(4*M),H_env); %size(sys)=L-i-3��size(env)=i+1��Hamiltonian
    end
    for i=L-m-3:L-4
        [H_V,H_D]=eig(H);
        H_d=real(diag(H_D));
        H_min_position=find(H_d==min(H_d));
        Psi_0=H_V(:,H_min_position(1));
        rho_env_W=rho_env_eigmat(Psi_0,M,M); %����rho_env_eigmat��������������env��size=i+1����M����
        dmrg.rho_env_W(2^(m+2)-2+2*M*(i-m-1)+1:2^(m+2)-2+2*M*(i-m),:)=rho_env_W;
        H_ne_env=Hamiltonian_ne_env(rho_env_W); %����Hamiltonian_ne_env��������env���ܶȾ�����ʸrho_env_W����i+1����H_ne_env
        dmrg.H_ne_env(2^(m+2)-3+(i-m)*2*M:2^(m+2)-4+(i-m+1)*2*M,:)=H_ne_env;
        H_env=rho_env_W'*(kron(eye(2),dmrg.H_env(2^(m+1)+(i-m-1)*M-1:2^(m+1)+(i-m)*M-2,:))+dmrg.H_ne_env(2^(m+2)-3+(i-m-1)*2*M:2^(m+2)-4+(i-m)*2*M,:))*rho_env_W; %��i����H_env��H_ne_env����i+1����H_env
        dmrg.H_env(2^(m+1)-1+(i-m)*M:2^(m+1)-2+(i-m+1)*M,:)=H_env;
        H_sys=dmrg.H_sys(2^(L-i-4+1)-1:2^(L-i-4+2)-2,1:2^(L-i-3));
        H_sys_ns=dmrg.H_sys_ns(2^(L-i-4+2)-3:2^(L-i-4+3)-4,1:2^(L-i-2));
        H=kron(H_sys,eye(4*M))+kron(H_sys_ns,eye(2*M))+kron(kron(eye(2^(L-i-3)),Hamiltonian(2)),eye(M))+kron(eye(2^(L-i-2)),H_ne_env)+kron(eye(2^(L-i-1)),H_env); %size(sys)=L-i-3��size(env)=i+1��Hamiltonian
    end
    for i=1:m-1
        [H_V,H_D]=eig(H);
        H_d=real(diag(H_D));
        H_min_position=find(H_d==min(H_d));
        Psi_0=H_V(:,H_min_position(1));
        rho_sys_W=rho_sys_eigmat(Psi_0,2^i,2^(i+1)); %����rho_sys_eigmat��������������sys��size=i+1����2^(i+1)����
        dmrg.rho_sys_W(2^(i+1)-1:2^(i+2)-2,:)=[rho_sys_W,zeros(2^(i+1),M-2^(i+1))];
        H_sys=rho_sys_W'*(kron(dmrg.H_sys(2^i-1:2^(i+1)-2,1:2^i),eye(2))+dmrg.H_sys_ns(2^(i+1)-3:2^(i+2)-4,1:2^(i+1)))*rho_sys_W; %��i����H_sys��H_sys_ns����i+1����H_sys
        dmrg.H_sys(2^(i+1)-1:2^(i+2)-2,:)=[H_sys,zeros(2^(i+1),M-2^(i+1))];
        H_sys_ns=Hamiltonian_sys_ns(rho_sys_W); %����Hamiltonian_sys_ns��������sys���ܶȾ�����ʸrho_sys_W����i+1����H_sys_ns
        dmrg.H_sys_ns(2^(i+2)-3:2^(i+3)-4,:)=[H_sys_ns,zeros(2^(i+2),2*M-2^(i+2))];
        H_ne_env=dmrg.H_ne_env(2^(m+2)-3+(L-i-4-m)*2*M:2^(m+2)-4+(L-i-4-m+1)*2*M,:);
        H_env=dmrg.H_env(2^(m+1)-1+(L-i-4-m)*M:2^(m+1)-2+(L-i-4-m+1)*M,:);
        H=kron(H_sys,eye(4*M))+kron(H_sys_ns,eye(2*M))+kron(kron(eye(2^(i+1)),Hamiltonian(2)),eye(M))+kron(eye(2^(i+2)),H_ne_env)+kron(eye(2^(i+3)),H_env); %size(sys)=i+1��size(env)=L-i-3��Hamiltonian
    end
    [H_V,H_D]=eig(H);
    H_d=real(diag(H_D));
    H_min_position=find(H_d==min(H_d));
    Psi_0=H_V(:,H_min_position(1));
    rho_sys_W=rho_sys_eigmat(Psi_0,2^m,M); %����rho_sys_eigmat��������������sys��size=m+1����M����
    dmrg.rho_sys_W(2^(m+1)-1:2^(m+2)-2,:)=rho_sys_W;
    H_sys=rho_sys_W'*(kron(dmrg.H_sys(2^m-1:2^(m+1)-2,1:2^m),eye(2))+dmrg.H_sys_ns(2^(m+1)-3:2^(m+2)-4,1:2^(m+1)))*rho_sys_W; %��m����H_sys��H_sys_ns����m+1����H_sys
    dmrg.H_sys(2^(m+1)-1:2^(m+1)+M-2,:)=H_sys;
    H_sys_ns=Hamiltonian_sys_ns(rho_sys_W); %����Hamiltonian_sys_ns��������sys���ܶȾ�����ʸrho_sys_W����m+1����H_sys_ns
    dmrg.H_sys_ns(2^(m+2)-3:2^(m+2)-4+2*M,:)=H_sys_ns;
    H_ne_env=dmrg.H_ne_env(2^(m+2)-3+(L-m-4-m)*2*M:2^(m+2)-4+(L-m-4-m+1)*2*M,:);
    H_env=dmrg.H_env(2^(m+1)-1+(L-m-4-m)*M:2^(m+1)-2+(L-m-4-m+1)*M,:);
    H=kron(H_sys,eye(4*M))+kron(H_sys_ns,eye(2*M))+kron(kron(eye(M),Hamiltonian(2)),eye(M))+kron(eye(2*M),H_ne_env)+kron(eye(4*M),H_env); %size(sys)=m+1��size(env)=L-i-3��Hamiltonian
    for i=m+1:L/2-2
        [H_V,H_D]=eig(H);
        H_d=real(diag(H_D));
        H_min_position=find(H_d==min(H_d));
        Psi_0=H_V(:,H_min_position(1));
        rho_sys_W=rho_sys_eigmat(Psi_0,M,M); %����rho_sys_eigmat��������������sys��size=i+1����M����
        dmrg.rho_sys_W(2^(m+2)-2+2*M*(i-m-1)+1:2^(m+2)-2+2*M*(i-m),:)=rho_sys_W;
        H_sys=rho_sys_W'*(kron(dmrg.H_sys(2^(m+1)+(i-m-1)*M-1:2^(m+1)+(i-m)*M-2,:),eye(2))+dmrg.H_sys_ns(2^(m+2)-3+(i-m-1)*2*M:2^(m+2)-4+(i-m)*2*M,:))*rho_sys_W; %��i����H_sys��H_sys_ns����i+1����H_sys
        dmrg.H_sys(2^(m+1)-1+(i-m)*M:2^(m+1)-2+(i-m+1)*M,:)=H_sys;
        H_sys_ns=Hamiltonian_sys_ns(rho_sys_W); %����Hamiltonian_sys_ns��������sys���ܶȾ�����ʸrho_sys_W����i+1����H_sys_ns
        dmrg.H_sys_ns(2^(m+2)-3+(i-m)*2*M:2^(m+2)-4+(i-m+1)*2*M,:)=H_sys_ns;
        H_ne_env=dmrg.H_ne_env(2^(m+2)-3+(L-i-4-m)*2*M:2^(m+2)-4+(L-i-4-m+1)*2*M,:);
        H_env=dmrg.H_env(2^(m+1)-1+(L-i-4-m)*M:2^(m+1)-2+(L-i-4-m+1)*M,:);
        H=kron(H_sys,eye(4*M))+kron(H_sys_ns,eye(2*M))+kron(kron(eye(M),Hamiltonian(2)),eye(M))+kron(eye(2*M),H_ne_env)+kron(eye(4*M),H_env); %size(sys)=i+1��size(env)=L-i-3��Hamiltonian
    end
end
[H_V,H_D]=eig(H);
H_d=real(diag(H_D));
E_0=min(H_d);
dmrg.E_0=E_0;
H_min_position=find(H_d==E_0);
n=size(H_min_position);
dmrg.Psi_0=cast(zeros(4*M^2*n(1),1),'like',1i);
for i=1:n(1)
    dmrg.Psi_0((i-1)*4*M^2+1:i*4*M^2,1)=H_V(:,H_min_position(i));
end
fprintf('E_0=%f\n',E_0);
if n(1)==1
    fprintf('No degeneracy in groud-state.\n');
else
    fprintf('%d-fold degeneracy in groud-state.\n',n(1));
end
toc;