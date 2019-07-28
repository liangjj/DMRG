function H=Hamiltonian(N) %一维N格点只考虑最近邻相互作用的Hamiltonian
    Jx=0.4;
    Jy=1;
    Jz=1;
    hbar=1;
    Pauli_x=[0, 1; 1, 0];
    Pauli_y=[0, -1i; 1i, 0];
    Pauli_z=[1, 0; 0, -1];
    hx=zeros(2^N);
    hy=zeros(2^N);
    hz=zeros(2^N);
    for i=1:N-1
        hx=hx+kron(kron(eye(2^(i-1)), kron(Pauli_x, Pauli_x)), eye(2^(N-i-1)));
        hy=hy+kron(kron(eye(2^(i-1)), kron(Pauli_y, Pauli_y)), eye(2^(N-i-1)));
        hz=hz+kron(kron(eye(2^(i-1)), kron(Pauli_z, Pauli_z)), eye(2^(N-i-1)));
    end
    H=(hx*Jx+hy*Jy+hz*Jz)*hbar^2/4;