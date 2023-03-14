function [s_hat,H_hat,idices] = OMP(y,H,n)
[M,N]=size(H);
maxIter=N;
iterflag=1;
epsilon=1e-3;
r_i=y;
y_i=zeros(size(y));
idices=zeros(1,N);
iterNum=0;
H_i=zeros(size(H));
while iterflag
    iterNum=iterNum+1;
    r=abs(r_i'*H);
    idx=find(r==max(r));
    idices(idx)=1;
    H_i(:,idx)=H(:,idx);
    P_i=H_i*pinv(H_i);
    r_i=(eye(M)-P_i)*r_i;
    if norm(r_i,2)<epsilon||iterNum>=maxIter
        iterflag=0;
    end
end
H_hat=H_i;
s_hat=pinv(H_hat)*y;
end