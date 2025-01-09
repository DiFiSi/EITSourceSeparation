function [temp,tempT] = mkTemps(tempFile,tempJitter,fs)
    temps = readtable(tempFile);
    tempF0 = 1 / (100 / fs);
    tempT = (0:100 - 1)' / fs;
    temp = zeros(100,2);
    
    % VRS
    aVRS = rmmissing(table2array(temps(1,2:2:end)))';
    bVRS = rmmissing(table2array(temps(1,3:2:end)))';
    
    LVRS = length(aVRS);
    for i = 1:LVRS
        temp(:,1) = temp(:,1) + aVRS(i) * (1 + randn * tempJitter(1)) * cos(2 * pi * tempF0 * i * tempT) + ...
                                bVRS(i) * (1 + randn * tempJitter(1)) * sin(2 * pi * tempF0 * i * tempT);
    end
    
    % CRS
    aCRS = rmmissing(table2array(temps(2,2:2:end)))';
    bCRS = rmmissing(table2array(temps(2,3:2:end)))';
    
    LCRS = length(aCRS);
    for i = 1:LCRS
        temp(:,2) = temp(:,2) + aCRS(i) * (1 + randn * tempJitter(2)) * cos(2 * pi * tempF0 * i * tempT) + ...
                                bCRS(i) * (1 + randn * tempJitter(2)) * sin(2 * pi * tempF0 * i * tempT);
    end
    temp(:,2) = circshift(temp(:,2),-5);

    temp = temp ./ range(temp,1);
end