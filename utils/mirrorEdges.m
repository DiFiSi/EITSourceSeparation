function [sig, extN] = mirrorEdges(sig, extN)
    N = size(sig,1);
    if ~exist('extN','var')
        extN = round(0.25 * N);
    else
        extN = min(N,extN) - 1;
    end
    
    pref = flipud(sig(2:extN + 1, :));
    suff = flipud(sig(end - extN:end - 1,:));
    
    sig = [pref; sig; suff];
end