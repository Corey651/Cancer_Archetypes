function [H] = NMF_New_Weights(V,W,k)

% Computes the archetype scores from previously trained weights W

% V: m*n
% W m*k
% H k*n

N_iter=1000;
rng(0,'twister');
[m,n]=size(V);
H=rand(k,n)+eps;

for i=1:N_iter
    % update H
    WV=W'*V;
    WWH=W'*W*H;
    H=H.*(WV./WWH);
    H=H./sum(H,1); % Additional normalization step

end