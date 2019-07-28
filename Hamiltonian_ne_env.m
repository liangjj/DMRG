function H=Hamiltonian_ne_env(W) %根据env密度矩阵的本征矢W（m*n），生成2n*2n的ne_env哈密顿矩阵
size_W=size(W);
h=Hamiltonian(2);
h=[h(:,1:2);h(:,3:4)];
H0=zeros(4*size_W(2),size_W(2));
for k=1:4
    for i=1:size_W(2)
        for j=1:size_W(2)
            for p=1:size_W(1)/2
                H0(size_W(2)*(k-1)+i,j)=H0(size_W(2)*(k-1)+i,j)+W([p p+size_W(1)/2],i)'*h(2*k-1:2*k,:)*W([p p+size_W(1)/2],j);
            end
        end
    end
end
H=[H0(1:2*size_W(2),:),H0(2*size_W(2)+1:4*size_W(2),:)];