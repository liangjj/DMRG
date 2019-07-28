tic;
clear;
save DMRG.mat -v7.3;
dmrg=matfile('DMRG.mat','Writable',true);
L=100; %大于4的总格点数（偶数）
M=10; %sys或env最多保留2<=M<2^(L/2-2)个态
N=30; %sweep循环次数
m=floor(log(M)/log(2));
dmrg.H_sys=cast(zeros(2,M),'like',1i); %此行及以下5行，是为DMRG.mat文件中的6个矩阵变量（H_sys,H_sys_ns等）初始化：一是确定每个矩阵变量的类型为复数，保证可以在矩阵变量中写入复数；二是能量矩阵需要使用初始化的赋值。
dmrg.H_sys_ns=cast([Hamiltonian(2),zeros(4,2*M-4)],'like',1i);
dmrg.H_ne_env=cast([Hamiltonian(2),zeros(4,2*M-4)],'like',1i);
dmrg.H_env=cast(zeros(2,M),'like',1i);
dmrg.rho_sys_W=cast(zeros(2,M),'like',1i);
dmrg.rho_env_W=cast(zeros(2,M),'like',1i);
H=Hamiltonian(4);
for i=1:m-1 %对称长点，直到size(sys)=m
    [H_V,H_D]=eig(H);
    H_d=real(diag(H_D));
    H_min_position=find(H_d==min(H_d));
    Psi_0=H_V(:,H_min_position(1));
    rho_sys_W=rho_sys_eigmat(Psi_0,2^i,2^(i+1)); %调用rho_sys_eigmat函数，生成描述sys（size=i+1）的2^(i+1)个基
    dmrg.rho_sys_W(2^(i+1)-1:2^(i+2)-2,:)=[rho_sys_W,zeros(2^(i+1),M-2^(i+1))];
    rho_env_W=rho_env_eigmat(Psi_0,2^i,2^(i+1)); %调用rho_env_eigmat函数，生成描述env（size=i+1）的2^(i+1)个基
    dmrg.rho_env_W(2^(i+1)-1:2^(i+2)-2,:)=[rho_env_W,zeros(2^(i+1),M-2^(i+1))];
    H_sys=rho_sys_W'*(kron(dmrg.H_sys(2^i-1:2^(i+1)-2,1:2^i),eye(2))+dmrg.H_sys_ns(2^(i+1)-3:2^(i+2)-4,1:2^(i+1)))*rho_sys_W; %用i格点的H_sys和H_sys_ns生成i+1格点的H_sys
    dmrg.H_sys(2^(i+1)-1:2^(i+2)-2,:)=[H_sys,zeros(2^(i+1),M-2^(i+1))];
    H_sys_ns=Hamiltonian_sys_ns(rho_sys_W); %调用Hamiltonian_sys_ns函数，由sys的密度矩阵本征矢rho_sys_W生成i+1格点的H_sys_ns
    dmrg.H_sys_ns(2^(i+2)-3:2^(i+3)-4,:)=[H_sys_ns,zeros(2^(i+2),2*M-2^(i+2))];
    H_ne_env=Hamiltonian_ne_env(rho_env_W); %调用Hamiltonian_ne_env函数，由env的密度矩阵本征矢rho_env_W生成i+1格点的H_ne_env
    dmrg.H_ne_env(2^(i+2)-3:2^(i+3)-4,:)=[H_ne_env,zeros(2^(i+2),2*M-2^(i+2))];
    H_env=rho_env_W'*(kron(eye(2),dmrg.H_env(2^i-1:2^(i+1)-2,1:2^i))+dmrg.H_ne_env(2^(i+1)-3:2^(i+2)-4,1:2^(i+1)))*rho_env_W; %用i格点的H_env和H_ne_env生成i+1格点的H_env
    dmrg.H_env(2^(i+1)-1:2^(i+2)-2,:)=[H_env,zeros(2^(i+1),M-2^(i+1))];
    H=kron(H_sys,eye(2^(i+3)))+kron(H_sys_ns,eye(2^(i+2)))+kron(kron(eye(2^(i+1)),Hamiltonian(2)),eye(2^(i+1)))+kron(eye(2^(i+2)),H_ne_env)+kron(eye(2^(i+3)),H_env); %size(sys)=size(env)=i+1的Hamiltonian
end
[H_V,H_D]=eig(H); %此行及以下rho_..._W特殊，对称长1个点，与上下均不同
H_d=real(diag(H_D));
H_min_position=find(H_d==min(H_d));
Psi_0=H_V(:,H_min_position(1));
rho_sys_W=rho_sys_eigmat(Psi_0,2^m,M); %调用rho_sys_eigmat函数，生成描述sys（size=m+1）的M个基
dmrg.rho_sys_W(2^(m+1)-1:2^(m+2)-2,:)=rho_sys_W;
rho_env_W=rho_env_eigmat(Psi_0,2^m,M); %调用rho_env_eigmat函数，生成描述env（size=m+1）的M个基
dmrg.rho_env_W(2^(m+1)-1:2^(m+2)-2,:)=rho_env_W;
H_sys=rho_sys_W'*(kron(dmrg.H_sys(2^m-1:2^(m+1)-2,1:2^m),eye(2))+dmrg.H_sys_ns(2^(m+1)-3:2^(m+2)-4,1:2^(m+1)))*rho_sys_W; %用m格点的H_sys和H_sys_ns生成m+1格点的H_sys
dmrg.H_sys(2^(m+1)-1:2^(m+1)+M-2,:)=H_sys;
H_sys_ns=Hamiltonian_sys_ns(rho_sys_W); %调用Hamiltonian_sys_ns函数，由sys的密度矩阵本征矢rho_sys_W生成m+1格点的H_sys_ns
dmrg.H_sys_ns(2^(m+2)-3:2^(m+2)-4+2*M,:)=H_sys_ns;
H_ne_env=Hamiltonian_ne_env(rho_env_W); %调用Hamiltonian_ne_env函数，由env的密度矩阵本征矢rho_env_W生成m+1格点的H_ne_env
dmrg.H_ne_env(2^(m+2)-3:2^(m+2)-4+2*M,:)=H_ne_env;
H_env=rho_env_W'*(kron(eye(2),dmrg.H_env(2^m-1:2^(m+1)-2,1:2^m))+dmrg.H_ne_env(2^(m+1)-3:2^(m+2)-4,1:2^(m+1)))*rho_env_W; %用m格点的H_env和H_ne_env生成m+1格点的H_env
dmrg.H_env(2^(m+1)-1:2^(m+1)+M-2,:)=H_env;
H=kron(H_sys,eye(4*M))+kron(H_sys_ns,eye(2*M))+kron(kron(eye(M),Hamiltonian(2)),eye(M))+kron(eye(2*M),H_ne_env)+kron(eye(4*M),H_env); %size(sys)=size(env)=m+1的Hamiltonian  %此行以上rho_..._W特殊，对称长1个点，与上下均不同
for i=m+1:L/2-2 %对称长点
    [H_V,H_D]=eig(H);
    H_d=real(diag(H_D));
    H_min_position=find(H_d==min(H_d));
    Psi_0=H_V(:,H_min_position(1));
    rho_sys_W=rho_sys_eigmat(Psi_0,M,M); %调用rho_sys_eigmat函数，生成描述sys（size=i+1）的M个基
    dmrg.rho_sys_W(2^(m+2)-2+2*M*(i-m-1)+1:2^(m+2)-2+2*M*(i-m),:)=rho_sys_W;
    rho_env_W=rho_env_eigmat(Psi_0,M,M); %调用rho_env_eigmat函数，生成描述env（size=i+1）的M个基
    dmrg.rho_env_W(2^(m+2)-2+2*M*(i-m-1)+1:2^(m+2)-2+2*M*(i-m),:)=rho_env_W;
    H_sys=rho_sys_W'*(kron(dmrg.H_sys(2^(m+1)+(i-m-1)*M-1:2^(m+1)+(i-m)*M-2,:),eye(2))+dmrg.H_sys_ns(2^(m+2)-3+(i-m-1)*2*M:2^(m+2)-4+(i-m)*2*M,:))*rho_sys_W; %用i格点的H_sys和H_sys_ns生成i+1格点的H_sys
    dmrg.H_sys(2^(m+1)-1+(i-m)*M:2^(m+1)-2+(i-m+1)*M,:)=H_sys;
    H_sys_ns=Hamiltonian_sys_ns(rho_sys_W); %调用Hamiltonian_sys_ns函数，由sys的密度矩阵本征矢rho_sys_W生成i+1格点的H_sys_ns
    dmrg.H_sys_ns(2^(m+2)-3+(i-m)*2*M:2^(m+2)-4+(i-m+1)*2*M,:)=H_sys_ns;
    H_ne_env=Hamiltonian_ne_env(rho_env_W); %调用Hamiltonian_ne_env函数，由env的密度矩阵本征矢rho_env_W生成i+1格点的H_ne_env
    dmrg.H_ne_env(2^(m+2)-3+(i-m)*2*M:2^(m+2)-4+(i-m+1)*2*M,:)=H_ne_env;
    H_env=rho_env_W'*(kron(eye(2),dmrg.H_env(2^(m+1)+(i-m-1)*M-1:2^(m+1)+(i-m)*M-2,:))+dmrg.H_ne_env(2^(m+2)-3+(i-m-1)*2*M:2^(m+2)-4+(i-m)*2*M,:))*rho_env_W; %用i格点的H_env和H_ne_env生成i+1格点的H_env
    dmrg.H_env(2^(m+1)-1+(i-m)*M:2^(m+1)-2+(i-m+1)*M,:)=H_env;
    H=kron(H_sys,eye(4*M))+kron(H_sys_ns,eye(2*M))+kron(kron(eye(M),Hamiltonian(2)),eye(M))+kron(eye(2*M),H_ne_env)+kron(eye(4*M),H_env); %size(sys)=size(env)=i+1的Hamiltonian
end
for j=1:N
    for i=L/2-1:L-m-4 %从中间向右sweep到size(sys)=L-m-3
        [H_V,H_D]=eig(H);
        H_d=real(diag(H_D));
        H_min_position=find(H_d==min(H_d));
        Psi_0=H_V(:,H_min_position(1));
        rho_sys_W=rho_sys_eigmat(Psi_0,M,M); %调用rho_sys_eigmat函数，生成描述sys（size=i+1）的M个基
        dmrg.rho_sys_W(2^(m+2)-2+2*M*(i-m-1)+1:2^(m+2)-2+2*M*(i-m),:)=rho_sys_W;
        H_sys=rho_sys_W'*(kron(dmrg.H_sys(2^(m+1)+(i-m-1)*M-1:2^(m+1)+(i-m)*M-2,:),eye(2))+dmrg.H_sys_ns(2^(m+2)-3+(i-m-1)*2*M:2^(m+2)-4+(i-m)*2*M,:))*rho_sys_W; %用i格点的H_sys和H_sys_ns生成i+1格点的H_sys
        dmrg.H_sys(2^(m+1)-1+(i-m)*M:2^(m+1)-2+(i-m+1)*M,:)=H_sys;
        H_sys_ns=Hamiltonian_sys_ns(rho_sys_W); %调用Hamiltonian_sys_ns函数，由sys的密度矩阵本征矢rho_sys_W生成i+1格点的H_sys_ns
        dmrg.H_sys_ns(2^(m+2)-3+(i-m)*2*M:2^(m+2)-4+(i-m+1)*2*M,:)=H_sys_ns;
        H_ne_env=dmrg.H_ne_env(2^(m+2)-3+(L-i-4-m)*2*M:2^(m+2)-4+(L-i-4-m+1)*2*M,:);
        H_env=dmrg.H_env(2^(m+1)-1+(L-i-4-m)*M:2^(m+1)-2+(L-i-4-m+1)*M,:);
        H=kron(H_sys,eye(4*M))+kron(H_sys_ns,eye(2*M))+kron(kron(eye(M),Hamiltonian(2)),eye(M))+kron(eye(2*M),H_ne_env)+kron(eye(4*M),H_env); %size(sys)=i+1，size(env)=L-i-3的Hamiltonian
    end
    for i=L-m-3:L-4 %继续向右sweep
        [H_V,H_D]=eig(H);
        H_d=real(diag(H_D));
        H_min_position=find(H_d==min(H_d));
        Psi_0=H_V(:,H_min_position(1));
        rho_sys_W=rho_sys_eigmat(Psi_0,M,M); %调用rho_sys_eigmat函数，生成描述sys（size=i+1）的M个基
        dmrg.rho_sys_W(2^(m+2)-2+2*M*(i-m-1)+1:2^(m+2)-2+2*M*(i-m),:)=rho_sys_W;
        H_sys=rho_sys_W'*(kron(dmrg.H_sys(2^(m+1)+(i-m-1)*M-1:2^(m+1)+(i-m)*M-2,:),eye(2))+dmrg.H_sys_ns(2^(m+2)-3+(i-m-1)*2*M:2^(m+2)-4+(i-m)*2*M,:))*rho_sys_W; %用i格点的H_sys和H_sys_ns生成i+1格点的H_sys
        dmrg.H_sys(2^(m+1)-1+(i-m)*M:2^(m+1)-2+(i-m+1)*M,:)=H_sys;
        H_sys_ns=Hamiltonian_sys_ns(rho_sys_W); %调用Hamiltonian_sys_ns函数，由sys的密度矩阵本征矢rho_sys_W生成i+1格点的H_sys_ns
        dmrg.H_sys_ns(2^(m+2)-3+(i-m)*2*M:2^(m+2)-4+(i-m+1)*2*M,:)=H_sys_ns;
        H_ne_env=dmrg.H_ne_env(2^(L-i-4+2)-3:2^(L-i-4+3)-4,1:2^(L-i-2));
        H_env=dmrg.H_env(2^(L-i-4+1)-1:2^(L-i-4+2)-2,1:2^(L-i-3));
        H=kron(H_sys,eye(2^(L-i-1)))+kron(H_sys_ns,eye(2^(L-i-2)))+kron(kron(eye(M),Hamiltonian(2)),eye(2^(L-i-3)))+kron(eye(M*2),H_ne_env)+kron(eye(M*4),H_env); %size(sys)=i+1，size(env)=L-i-3的Hamiltonian
    end
    for i=1:m-1 %从右向左sweep到size(env)=m
        [H_V,H_D]=eig(H);
        H_d=real(diag(H_D));
        H_min_position=find(H_d==min(H_d));
        Psi_0=H_V(:,H_min_position(1));
        rho_env_W=rho_env_eigmat(Psi_0,2^i,2^(i+1)); %调用rho_env_eigmat函数，生成描述env（size=i+1）的2^(i+1)个基
        dmrg.rho_env_W(2^(i+1)-1:2^(i+2)-2,:)=[rho_env_W,zeros(2^(i+1),M-2^(i+1))];
        H_ne_env=Hamiltonian_ne_env(rho_env_W); %调用Hamiltonian_ne_env函数，由env的密度矩阵本征矢rho_env_W生成i+1格点的H_ne_env
        dmrg.H_ne_env(2^(i+2)-3:2^(i+3)-4,:)=[H_ne_env,zeros(2^(i+2),2*M-2^(i+2))];
        H_env=rho_env_W'*(kron(eye(2),dmrg.H_env(2^i-1:2^(i+1)-2,1:2^i))+dmrg.H_ne_env(2^(i+1)-3:2^(i+2)-4,1:2^(i+1)))*rho_env_W; %用i格点的H_env和H_ne_env生成i+1格点的H_env
        dmrg.H_env(2^(i+1)-1:2^(i+2)-2,:)=[H_env,zeros(2^(i+1),M-2^(i+1))];
        H_sys=dmrg.H_sys(2^(m+1)-1+(L-i-4-m)*M:2^(m+1)-2+(L-i-4-m+1)*M,:);
        H_sys_ns=dmrg.H_sys_ns(2^(m+2)-3+(L-i-4-m)*2*M:2^(m+2)-4+(L-i-4-m+1)*2*M,:);
        H=kron(H_sys,eye(2^(i+3)))+kron(H_sys_ns,eye(2^(i+2)))+kron(kron(eye(M),Hamiltonian(2)),eye(2^(i+1)))+kron(eye(2*M),H_ne_env)+kron(eye(4*M),H_env); %size(sys)=L-i-3，size(env)=i+1的Hamiltonian
    end
    [H_V,H_D]=eig(H);
    H_d=real(diag(H_D));
    H_min_position=find(H_d==min(H_d));
    Psi_0=H_V(:,H_min_position(1));
    rho_env_W=rho_env_eigmat(Psi_0,2^m,M); %调用rho_env_eigmat函数，生成描述env（size=m+1）的M个基
    dmrg.rho_env_W(2^(m+1)-1:2^(m+2)-2,:)=rho_env_W;
    H_ne_env=Hamiltonian_ne_env(rho_env_W); %调用Hamiltonian_ne_env函数，由env的密度矩阵本征矢rho_env_W生成m+1格点的H_ne_env
    dmrg.H_ne_env(2^(m+2)-3:2^(m+2)-4+2*M,:)=H_ne_env;
    H_env=rho_env_W'*(kron(eye(2),dmrg.H_env(2^m-1:2^(m+1)-2,1:2^m))+dmrg.H_ne_env(2^(m+1)-3:2^(m+2)-4,1:2^(m+1)))*rho_env_W; %用m格点的H_env和H_ne_env生成m+1格点的H_env
    dmrg.H_env(2^(m+1)-1:2^(m+1)+M-2,:)=H_env;
    H_sys=dmrg.H_sys(2^(m+1)-1+(L-m-4-m)*M:2^(m+1)-2+(L-m-4-m+1)*M,:);
    H_sys_ns=dmrg.H_sys_ns(2^(m+2)-3+(L-m-4-m)*2*M:2^(m+2)-4+(L-m-4-m+1)*2*M,:);
    H=kron(H_sys,eye(4*M))+kron(H_sys_ns,eye(2*M))+kron(kron(eye(M),Hamiltonian(2)),eye(M))+kron(eye(2*M),H_ne_env)+kron(eye(4*M),H_env); %size(sys)=L-i-3，size(env)=m+1的Hamiltonian
    for i=m+1:L-m-4
        [H_V,H_D]=eig(H);
        H_d=real(diag(H_D));
        H_min_position=find(H_d==min(H_d));
        Psi_0=H_V(:,H_min_position(1));
        rho_env_W=rho_env_eigmat(Psi_0,M,M); %调用rho_env_eigmat函数，生成描述env（size=i+1）的M个基
        dmrg.rho_env_W(2^(m+2)-2+2*M*(i-m-1)+1:2^(m+2)-2+2*M*(i-m),:)=rho_env_W;
        H_ne_env=Hamiltonian_ne_env(rho_env_W); %调用Hamiltonian_ne_env函数，由env的密度矩阵本征矢rho_env_W生成i+1格点的H_ne_env
        dmrg.H_ne_env(2^(m+2)-3+(i-m)*2*M:2^(m+2)-4+(i-m+1)*2*M,:)=H_ne_env;
        H_env=rho_env_W'*(kron(eye(2),dmrg.H_env(2^(m+1)+(i-m-1)*M-1:2^(m+1)+(i-m)*M-2,:))+dmrg.H_ne_env(2^(m+2)-3+(i-m-1)*2*M:2^(m+2)-4+(i-m)*2*M,:))*rho_env_W; %用i格点的H_env和H_ne_env生成i+1格点的H_env
        dmrg.H_env(2^(m+1)-1+(i-m)*M:2^(m+1)-2+(i-m+1)*M,:)=H_env;
        H_sys=dmrg.H_sys(2^(m+1)-1+(L-i-4-m)*M:2^(m+1)-2+(L-i-4-m+1)*M,:);
        H_sys_ns=dmrg.H_sys_ns(2^(m+2)-3+(L-i-4-m)*2*M:2^(m+2)-4+(L-i-4-m+1)*2*M,:);
        H=kron(H_sys,eye(4*M))+kron(H_sys_ns,eye(2*M))+kron(kron(eye(M),Hamiltonian(2)),eye(M))+kron(eye(2*M),H_ne_env)+kron(eye(4*M),H_env); %size(sys)=L-i-3，size(env)=i+1的Hamiltonian
    end
    for i=L-m-3:L-4
        [H_V,H_D]=eig(H);
        H_d=real(diag(H_D));
        H_min_position=find(H_d==min(H_d));
        Psi_0=H_V(:,H_min_position(1));
        rho_env_W=rho_env_eigmat(Psi_0,M,M); %调用rho_env_eigmat函数，生成描述env（size=i+1）的M个基
        dmrg.rho_env_W(2^(m+2)-2+2*M*(i-m-1)+1:2^(m+2)-2+2*M*(i-m),:)=rho_env_W;
        H_ne_env=Hamiltonian_ne_env(rho_env_W); %调用Hamiltonian_ne_env函数，由env的密度矩阵本征矢rho_env_W生成i+1格点的H_ne_env
        dmrg.H_ne_env(2^(m+2)-3+(i-m)*2*M:2^(m+2)-4+(i-m+1)*2*M,:)=H_ne_env;
        H_env=rho_env_W'*(kron(eye(2),dmrg.H_env(2^(m+1)+(i-m-1)*M-1:2^(m+1)+(i-m)*M-2,:))+dmrg.H_ne_env(2^(m+2)-3+(i-m-1)*2*M:2^(m+2)-4+(i-m)*2*M,:))*rho_env_W; %用i格点的H_env和H_ne_env生成i+1格点的H_env
        dmrg.H_env(2^(m+1)-1+(i-m)*M:2^(m+1)-2+(i-m+1)*M,:)=H_env;
        H_sys=dmrg.H_sys(2^(L-i-4+1)-1:2^(L-i-4+2)-2,1:2^(L-i-3));
        H_sys_ns=dmrg.H_sys_ns(2^(L-i-4+2)-3:2^(L-i-4+3)-4,1:2^(L-i-2));
        H=kron(H_sys,eye(4*M))+kron(H_sys_ns,eye(2*M))+kron(kron(eye(2^(L-i-3)),Hamiltonian(2)),eye(M))+kron(eye(2^(L-i-2)),H_ne_env)+kron(eye(2^(L-i-1)),H_env); %size(sys)=L-i-3，size(env)=i+1的Hamiltonian
    end
    for i=1:m-1
        [H_V,H_D]=eig(H);
        H_d=real(diag(H_D));
        H_min_position=find(H_d==min(H_d));
        Psi_0=H_V(:,H_min_position(1));
        rho_sys_W=rho_sys_eigmat(Psi_0,2^i,2^(i+1)); %调用rho_sys_eigmat函数，生成描述sys（size=i+1）的2^(i+1)个基
        dmrg.rho_sys_W(2^(i+1)-1:2^(i+2)-2,:)=[rho_sys_W,zeros(2^(i+1),M-2^(i+1))];
        H_sys=rho_sys_W'*(kron(dmrg.H_sys(2^i-1:2^(i+1)-2,1:2^i),eye(2))+dmrg.H_sys_ns(2^(i+1)-3:2^(i+2)-4,1:2^(i+1)))*rho_sys_W; %用i格点的H_sys和H_sys_ns生成i+1格点的H_sys
        dmrg.H_sys(2^(i+1)-1:2^(i+2)-2,:)=[H_sys,zeros(2^(i+1),M-2^(i+1))];
        H_sys_ns=Hamiltonian_sys_ns(rho_sys_W); %调用Hamiltonian_sys_ns函数，由sys的密度矩阵本征矢rho_sys_W生成i+1格点的H_sys_ns
        dmrg.H_sys_ns(2^(i+2)-3:2^(i+3)-4,:)=[H_sys_ns,zeros(2^(i+2),2*M-2^(i+2))];
        H_ne_env=dmrg.H_ne_env(2^(m+2)-3+(L-i-4-m)*2*M:2^(m+2)-4+(L-i-4-m+1)*2*M,:);
        H_env=dmrg.H_env(2^(m+1)-1+(L-i-4-m)*M:2^(m+1)-2+(L-i-4-m+1)*M,:);
        H=kron(H_sys,eye(4*M))+kron(H_sys_ns,eye(2*M))+kron(kron(eye(2^(i+1)),Hamiltonian(2)),eye(M))+kron(eye(2^(i+2)),H_ne_env)+kron(eye(2^(i+3)),H_env); %size(sys)=i+1，size(env)=L-i-3的Hamiltonian
    end
    [H_V,H_D]=eig(H);
    H_d=real(diag(H_D));
    H_min_position=find(H_d==min(H_d));
    Psi_0=H_V(:,H_min_position(1));
    rho_sys_W=rho_sys_eigmat(Psi_0,2^m,M); %调用rho_sys_eigmat函数，生成描述sys（size=m+1）的M个基
    dmrg.rho_sys_W(2^(m+1)-1:2^(m+2)-2,:)=rho_sys_W;
    H_sys=rho_sys_W'*(kron(dmrg.H_sys(2^m-1:2^(m+1)-2,1:2^m),eye(2))+dmrg.H_sys_ns(2^(m+1)-3:2^(m+2)-4,1:2^(m+1)))*rho_sys_W; %用m格点的H_sys和H_sys_ns生成m+1格点的H_sys
    dmrg.H_sys(2^(m+1)-1:2^(m+1)+M-2,:)=H_sys;
    H_sys_ns=Hamiltonian_sys_ns(rho_sys_W); %调用Hamiltonian_sys_ns函数，由sys的密度矩阵本征矢rho_sys_W生成m+1格点的H_sys_ns
    dmrg.H_sys_ns(2^(m+2)-3:2^(m+2)-4+2*M,:)=H_sys_ns;
    H_ne_env=dmrg.H_ne_env(2^(m+2)-3+(L-m-4-m)*2*M:2^(m+2)-4+(L-m-4-m+1)*2*M,:);
    H_env=dmrg.H_env(2^(m+1)-1+(L-m-4-m)*M:2^(m+1)-2+(L-m-4-m+1)*M,:);
    H=kron(H_sys,eye(4*M))+kron(H_sys_ns,eye(2*M))+kron(kron(eye(M),Hamiltonian(2)),eye(M))+kron(eye(2*M),H_ne_env)+kron(eye(4*M),H_env); %size(sys)=m+1，size(env)=L-i-3的Hamiltonian
    for i=m+1:L/2-2
        [H_V,H_D]=eig(H);
        H_d=real(diag(H_D));
        H_min_position=find(H_d==min(H_d));
        Psi_0=H_V(:,H_min_position(1));
        rho_sys_W=rho_sys_eigmat(Psi_0,M,M); %调用rho_sys_eigmat函数，生成描述sys（size=i+1）的M个基
        dmrg.rho_sys_W(2^(m+2)-2+2*M*(i-m-1)+1:2^(m+2)-2+2*M*(i-m),:)=rho_sys_W;
        H_sys=rho_sys_W'*(kron(dmrg.H_sys(2^(m+1)+(i-m-1)*M-1:2^(m+1)+(i-m)*M-2,:),eye(2))+dmrg.H_sys_ns(2^(m+2)-3+(i-m-1)*2*M:2^(m+2)-4+(i-m)*2*M,:))*rho_sys_W; %用i格点的H_sys和H_sys_ns生成i+1格点的H_sys
        dmrg.H_sys(2^(m+1)-1+(i-m)*M:2^(m+1)-2+(i-m+1)*M,:)=H_sys;
        H_sys_ns=Hamiltonian_sys_ns(rho_sys_W); %调用Hamiltonian_sys_ns函数，由sys的密度矩阵本征矢rho_sys_W生成i+1格点的H_sys_ns
        dmrg.H_sys_ns(2^(m+2)-3+(i-m)*2*M:2^(m+2)-4+(i-m+1)*2*M,:)=H_sys_ns;
        H_ne_env=dmrg.H_ne_env(2^(m+2)-3+(L-i-4-m)*2*M:2^(m+2)-4+(L-i-4-m+1)*2*M,:);
        H_env=dmrg.H_env(2^(m+1)-1+(L-i-4-m)*M:2^(m+1)-2+(L-i-4-m+1)*M,:);
        H=kron(H_sys,eye(4*M))+kron(H_sys_ns,eye(2*M))+kron(kron(eye(M),Hamiltonian(2)),eye(M))+kron(eye(2*M),H_ne_env)+kron(eye(4*M),H_env); %size(sys)=i+1，size(env)=L-i-3的Hamiltonian
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