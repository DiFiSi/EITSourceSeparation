function out = dPCA(X,hyper)
    if ~exist('hyper.transWindow','var')
       hyper.transWin = 0.2;
    end
    
    [nSamples,~] = size(X);

    %% Prepare time delays
    deltaV = (0:hyper.dfDeltaV * hyper.fs:hyper.maxDeltaV * hyper.fs)';
    if isempty(deltaV)
        deltaV = 0;
    else
        deltaV = unique(round([-deltaV(2:end);deltaV]));
    end
    
    deltaC = (0:hyper.dfDeltaC * hyper.fs:hyper.maxDeltaC * hyper.fs)';
    if isempty(deltaV)
        deltaV = 0;
    else
        deltaC = unique(round([-deltaC(2:end);deltaC]));
    end

    %% First Approximation
    % PCA   - keep first component
    %       - make delayed matrix
    S = PCA(X,1);
    Bv = S(:,1);

    nDels = length(deltaV);
    BvDel = zeros(size(Bv,1),nDels);
    for i = 1:nDels
        BvDel(:,i) = circshift(Bv,deltaV(i),1);
    end

    % LMS   - project first PCA component onto data
    Xvp = BvDel * (BvDel \ X);

    % Remove VRS projection from data
    Xcp = X - Xvp;
    
    %% Second Approximation
    % BPF   - bandpass filter data w/o VRS projection
    Xcpbp = mirroredFilt(Xcp,hyper.dHp);
    Xcpbp = mirroredFilt(Xcpbp,hyper.dLpNoise);

    % PCA   - keep first two components and 
    %       - set delayed matrix
    S = PCA(Xcpbp,2);
    Bc = S(:,1:2);

    nDels = length(deltaC);
    BcDel = [];
    for i = 1:nDels
        BcDel = [BcDel,circshift(Bc,deltaC(i),1)];
    end
    
    % LMS   - project first two components onto data w/o VRS projection
    %       - add projection to CRS estimation
    Xcpp1 = BcDel * (BcDel \ Xcp);

    % LMS   - project first two components onto VRS projection
    %       - add CRS projection onto CRS estimation
    %       - subtract CRS projection from VRS projection to get VRS
    %       estimation
    Xcpp2 = BcDel * (BcDel \ Xvp);

    Xc = Xcpp1 + Xcpp2;
    Xv = Xvp - Xcpp2;

    % Collect outputs
    out.CRShat = Xc;
    out.VRShat = Xv;
    out.INThat = nan(size(out.VRShat));
    out.MIXhat = out.CRShat + out.VRShat;
    out.startIdx = 1;
    out.endIdx = nSamples;
end

function [U,S,Vt,idx] = PCA(data,thresh)
    [U,S,V] = svd(data,'econ');
    if thresh < 1
        s = diag(S);
        cS = cumsum(s ./ sum(s));
        idx = 1:find(cS >= thresh,1,'first');
    else
       idx = 1:thresh; 
    end
    
    Vt = V';
end

function [filtData, resData] = filtMix(mixMat, idx, sources, latent)
    select = zeros(1, size(sources, 2));
    select(idx) = 1;
    filtData = sources * latent * diag(select) * mixMat;
    resData = sources * latent * diag(~select) * mixMat;
end
