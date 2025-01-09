function out = optFilt(data,f0,hyper)
    % Initialization of sampling frequency and size parameters
    fs = hyper.fs;
    ts = 1 / fs;
    [S,~] = size(data);

    % Extract hyperparameters
    N = hyper.N;
    M = hyper.M;
    L = hyper.L;
    
    % Constants
    oneMat = ones(L,1); % column selector
    
    % Preallocate output
    hAll = nan(S,M);
    aAll = nan(S,L);
    yMod = nan(S,1); % initialize modelled output
    yFilt = nan(S,1); % initialize filtered output

    % Get frequency vectors
    nSources = size(f0,2);
    f = [];
    for i = 1:nSources
        f = [f,f0(:,i) .* (1:round(L / nSources))];
    end

    % Loop over each window
    for n = N + 1:S
        % Set signal window
        winn = (n - N + 1:n)';
        tn = (0:N - 1)' * ts;
        xn = data(winn,:);
        fn = f(winn,:);
        thetan = 2 * pi * cumtrapz(tn,fn);
    
        % Set filter window
        winm = (n - M + 1:n)';
        tm = (N - M:N - 1)' * ts;
        xm = data(winm,:);
        fm = f(winm,:); 
        thetam = thetan((N - M + 1:N),:);
    
        % Build harmonic model
        Zn   = exp(1i * thetam);
            
        % Build model for present
        wn   = exp(1i * thetam(end,:)');

        % Calculate covariance matrices over sliding window
        [Rx,W,G] = getCovs(thetan,xn,M,L,N);
    
        % Approximate noise covariance 
        Q = Rx - G' * (W \ G); 
    
        % Get filter h and model coefficients a
        lambda = (Zn' / Q * Zn) \ oneMat;
        h   = Q \ Zn * lambda;
        a   = W \ G * h;
    
        % Reconstruct signals
        hAll(n - N + 1,:) = real(h);
        aAll(n - N + 1,:) = real(a);
        yMod(n - N + 1) = real(2 * a' * wn);
        yFilt(n - N + 1) = real(2 * h' * xm);
    
        % Verbose
        fprintf("Window " + num2str(n) + " / " + num2str(S) + " \n");
    end

    % Clip signals
    [~,remIdx] = rmmissing(yFilt);
    keepIdx = find(~remIdx);
    out.yFilt = yFilt;
    out.yMod = yMod;
    out.aAll = aAll;
    out.hAll = hAll;
    out.startIdx = keepIdx(1);
    out.endIdx = keepIdx(end);
end

