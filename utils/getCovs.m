function [Rx,W,G] = getCovs(thetan,xn,M,L,N)
    Rx = zeros(M, M);
    W = zeros(L, L);
    G = zeros(L, M);

    avg = N - M + 1;
    for i = M:N
        wi = exp(1i * thetan(i,:)');
        xi = xn(i - M + 1:i,:);

        Rx = Rx + xi * xi';
        W = W + wi * wi';
        G = G + wi * xi';
    end
    Rx = Rx / avg + eps * eye(M,M);
    W = W / avg;
    G = G / avg;
end
