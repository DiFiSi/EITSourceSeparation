function out = ceemdanFilt(data, hyper)
    [nSamples,nChannels] = size(data);

    % CEEMDAN Parameters
    NR = 0.5;
    maxReal = 40;
    maxIter = 1e6;
    SNRFlag = true;

    % Run EMD
    out.VRShat = zeros(nSamples,nChannels);
    out.CRShat = zeros(nSamples,nChannels);

    for i = 1:nChannels
        IMF = ceemdan(data(:,i),NR,maxReal,maxIter,SNRFlag);
        IMF = IMF';
        
        % Sort IMFs to VRS or CRS
        [out.VRShat(:,i),out.CRShat(:,i)] = sortIMF(IMF,hyper.fs);
    end

    out.MIXhat = out.VRShat + out.CRShat;
    out.startIdx = 1;
    out.endIdx = nSamples;
end

function [VRS,CRS,DC,e] = sortIMF(IMF,fs)
    [nSamples,nIMFs] = size(IMF);

    % Initialize outputs
    DC          = zeros(nSamples,1);
    VRS         = zeros(nSamples,1);
    CRS         = zeros(nSamples,1);
    e           = zeros(nSamples,1);

    % Perform FFT
    nBins       = 2 ^ nextpow2(min(12000,nSamples * 100));
    halfBins    = round(nBins / 2);
    
    Y           = fft(IMF, nBins) / 2;
    Y           = Y(1:halfBins + 1, :);
    
    amp         = abs(Y);
    f           = fs * (0:halfBins)' / nBins;

    % Define decision frequency bands
    bandDC      = [0;1] / 60;
    bandVRS     = [1;30] / 60;
    bandCRS     = [30;200] / 60;
    bande       = [200 / 60; fs / 2];

    % Assign dominant frequencies
    for i = 1:nIMFs
        [pk,loc] = findpeaks(amp(:,i),f,"NPeaks",5,"SortStr","descend");
        fd       = pk' * loc / sum(pk);

        if bandDC(1) <= fd && fd < bandDC(2)
            DC  = DC + IMF(:,i);
        elseif bandVRS(1) <= fd && fd < bandVRS(2)
            VRS = VRS + IMF(:,i);
        elseif bandCRS(1) <= fd && fd < bandCRS(2)
            CRS = CRS + IMF(:,i);
        elseif bande(1) <= fd && fd < bande(2)
            e   = e + IMF(:,i);
        end
    end
end
