function out = bandFilt(data, hyper)
    % Bandpass filter signal
    out.VRShat = mirroredFilt(mirroredFilt(data,hyper.dHpDrift),hyper.dLp);
    out.CRShat = mirroredFilt(mirroredFilt(data,hyper.dLpNoise),hyper.dHp);

    % Collect results
    out.INThat = nan(size(out.VRShat));
    out.MIXhat = out.VRShat + out.CRShat;
    out.startIdx = 1;
    out.endIdx = size(data,1);
end

