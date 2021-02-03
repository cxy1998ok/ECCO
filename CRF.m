close all;
N=288; %ttotal
k=36; %blocksize, need k<=(N+1)/2
taomax=3; %say max lag is 3 years
tao=[0:taomax*12/k:taomax*12]; %lags 
P=nan(N-k+1,1);
S=nan(N-k+1,k);
for p=1:N-k+1
    for q=1:k
        S(p,q)=SAM(p-q+k);
    end
    P(p)=psi2(k+p-1);
end
g=S\P;
figure;
plot(g);
SST=zeros(1,N);
for t=k+1:N%start with k+1 to avoid negative indexes
    for i=1:length(g)
        SST(t)=SST(t)+g(i)*SAM(t-tao(i)); %convolution
    end
end
figure;
plot(SST);
