function H=Hamiltonian_sys_ns(W) %根据sys密度矩阵的本征矢W（m*n），生成2n*2n的sys_ns哈密顿矩阵
size_W=size(W);
h=Hamiltonian(2);
h=reshape(h,8,2);
h([2 3 6 7],:)=h([3 2 7 6],:);
H0=zeros(2*size_W(2));
h0=zeros(4,1);
for i=1:size_W(2)
    for j=1:size_W(2)
        for k=1:4
            for p=1:size_W(1)/2
                h0(k)=h0(k)+W(2*p-1:2*p,i)'*h(2*k-1:2*k,:)*W(2*p-1:2*p,j);
            end
        end
        H0(2*i-1:2*i,2*j-1:2*j)=reshape(h0,2,2);
        h0=zeros(4,1);
    end
end
H=H0;