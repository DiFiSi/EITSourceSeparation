function [VRS,CRS,t] = getResult(out,maskIdxs,t)
    startIdx = out.startIdx; endIdx = out.endIdx;
    VRS = out.VRShat(startIdx:endIdx,maskIdxs);
    CRS = out.CRShat(startIdx:endIdx,maskIdxs);
    t = t(startIdx:endIdx);
end