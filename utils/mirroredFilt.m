function [imgFilt,transWin] = mirroredFilt(img,filtObj)
    sz = size(img);
    if length(sz) > 2
        xSize = sz(1);
        ySize = sz(2); 
    
        imgFlat = reshape(img, xSize * ySize, [])';
    else
        imgFlat = img;
    end

    [imgFlat,transWin] = mirrorEdges(imgFlat, filtord(filtObj) * 100);
    imgFlatFilt = filtfilt(filtObj, imgFlat);
    imgFlatFilt = imgFlatFilt(transWin + 1:end - transWin,:);

    
    if length(sz) > 2
        imgFilt = reshape(imgFlatFilt', xSize, ySize, []);
    else
        imgFilt = imgFlatFilt;
    end
end

