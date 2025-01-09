function out = harmRegWin(data,f0,hyper)
    % Initialization of sampling frequency and size parameters
    fs = hyper.fs;
    [nSamples,nChannels] = size(data);

    % Preallocate result matrices with NaN to handle missing data
    VRShat = nan(nSamples,nChannels);
    CRShat = nan(nSamples,nChannels);
    L = hyper.nHarms; % Number of harmonics for ventilation and perfusion
    
    % Create frequency matrix and detect overlaps
    f = [zeros(nSamples,1),     ... % DC
         f0(:,1) .* (1:L(1)),   ... % CRS
         f0(:,2) .* (1:L(2))];      % VRS
    idxV = 1 + (1:L(1))';
    idxP = idxV(end) + (1:L(2))';

    % Time vector calculation
    t = (0:1/fs:(nSamples - 1/fs))';

    % Set up parameters for harmonic analysis
    winSize = hyper.nWin; % Initial size of the analysis window
    
    % Main processing loop over the data samples within the valid window range
    for n = winSize:nSamples
        nWin = n - winSize + 1:n;
        tWin = t(nWin);
        yWin = data(nWin,:);
        fWin = f(nWin,:);
        thetaWin = 2 * pi * cumtrapz(tWin,fWin);
        
        % Build model
        Z   = exp(1i * thetaWin);

        % Calculate coefficients
        a   = 2 * (Z \ yWin);

        % Calculate outputs
        VRShat(n,:) = real(Z(end,idxV) * a(idxV,:));
        CRShat(n,:) = real(Z(end,idxP) * a(idxP,:));

        % Verbose
        fprintf("Window " + num2str(n) + " / " + num2str(nSamples) + " \n");
    end
    
    % Collect results
    [~,remTF] = rmmissing(VRShat,1);
    keepIdx = find(~remTF);
    out.CRShat = CRShat;
    out.VRShat = VRShat;
    out.startIdx = keepIdx(1);
    out.endIdx = keepIdx(end);
end