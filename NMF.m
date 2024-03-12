function [W,H] = NMF(V,k)


% V: m*n
% W m*k
% H k*n
N_iter=1000;
rng(0,'twister');
[m,n]=size(V);
W=rand(m,k)+eps;
H=rand(k,n)+eps;

for i=1:N_iter
    % update H
    WV=W'*V;
    WWH=W'*W*H;
    H=H.*(WV./WWH);
    H=H./sum(H,1); % Additional normalization step

    % update W
    VH=V*H';
    WHH=W*(H*H');
    W=W.*(VH./WHH);
    


end