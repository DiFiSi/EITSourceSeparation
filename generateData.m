%% Set target directory
addpath(genpath(".\generator"));
dataFolder = ".\data\generated\";
qaFolder = ".\data\generated\QA\";

if ~exist(qaFolder,'dir')
    mkdir(qaFolder);
end

%% Read simulations file
simsFile = ".\generator\simulations.xlsx";
sims = readtable(simsFile);

%% Load IBS signals
load(".\data\example\ibsHeart.mat");
ibsHeart = y;
load(".\data\example\ibsLung.mat");
ibsLung = y;

clear y;

%% Parameter space
% Constants
fs      = 50;   % Hz
seg0Time = 15;  % s
segTime = 60;   % s
t   = (0:(seg0Time + segTime) * fs - 1)' / fs;
t0	= t(1:1 * seg0Time * fs);
t1	= t(1:1 * segTime * fs);

% Define PRQ as a function of HR
prq = @(HR) 0.024 * HR + 3.92;

% Simulation parameter bounds
bolusTp = [0;1;2];
bHR1    = linspace(sims.bHR1(1),sims.bHR1(2),sims.bHR1(3))';
bHR2    = linspace(sims.bHR2(1),sims.bHR2(2),sims.bHR2(3))';
kCard   = linspace(sims.kCard(1),sims.kCard(2),sims.kCard(3))';
FM      = linspace(sims.FM(1),sims.FM(2),sims.FM(3))';
AM      = linspace(sims.AM(1),sims.AM(2),sims.AM(3))';
SNR     = linspace(sims.SNR(1),sims.SNR(2),sims.SNR(3))';

% Simulation parameter delta
dbHR1   = diff(bHR1(1:2)) / range(bHR1);
dbHR2   = diff(bHR2(1:2)) / range(bHR2);
dkCard  = diff(kCard(1:2)) / range(kCard);
dFM     = diff(FM(1:2)) / range(FM);
dAM     = diff(AM(1:2)) / range(AM);
dSNR    = diff(SNR(1:2)) / range(SNR);

% Set up parameter space
[gbolusTp,gbHR1,gbHR2,gkCard,gFM,gAM,gSNR] = ndgrid(bolusTp,bHR1,bHR2,kCard,FM,AM,SNR);
gbolusTp = gbolusTp(:);
gbHR1 = gbHR1(:);
gbHR2 = gbHR2(:);
gkCard = gkCard(:);
gFM = gFM(:);
gAM = gAM(:);
gSNR = gSNR(:);

%% Set template variables
tempJitter = [0.15;0.15];
tempFile = ".\generator\templates.xlsx";

%% Create signals
nSims   = numel(gbHR1);
fid = fopen('.\generator\simulationLog.txt','w');

% Save parameter table
params = [gbolusTp,gbHR1,gbHR2,gkCard,gFM,gAM,gSNR];
dataFile = ".\generator\summary.csv";
if ~isfile(dataFile)
    % writematrix(params,dataFile);

    str = "SUMMARY PARAMETERS SAVED!";
    disp(str);
    fprintf(fid, str + "\n");
else
    str = "SUMMARY FILE EXISTS ALREADY!";
    disp(str);
    fprintf(fid, str + "\n");
end

rng(0);
startSim = 1;
left = round((nSims - 211010) / 2); 
for i = 1:nSims
    try 
        % Gather params for iteration
        bltp    = gbolusTp(i);
        bhr1    = min(max(gbHR1(i)	* (1 + randn * dbHR1 / 2)   ,sims.bHR1(1)),sims.bHR1(2));
        bhr2    = min(max(gbHR2(i)	* (1 + randn * dbHR2 / 2)   ,sims.bHR2(1)),sims.bHR2(2));
        kcard   = min(max(gkCard(i)	* (1 + randn * dkCard / 2)  ,sims.kCard(1)),sims.kCard(2));
        fm      = min(max(gFM(i)	* (1 + randn * dFM / 2)     ,sims.FM(1)),sims.FM(2));
        am      = min(max(gAM(i)	* (1 + randn * dAM / 2)     ,sims.AM(1)),sims.AM(2));
        snr     = gSNR(i);
        
        % Make signal
        % Mark test segment
        initT = zeros(size(t));
        initT(1:length(t0)) = 1;

        % Heart rate
        HR0 = bhr1 * ones(size(t0));
        HR1 = linspace(bhr1,bhr2,length(t1))';
        HR = [HR0;HR1];

        % Respiratory rate
        RR = 1 ./ prq(HR) .* HR;

        % Synthesize templates
        sim.HR              = HR;
        sim.RR              = RR;
        sim.MWFreq          = 6;
        sim.RSAFM           = fm;
        sim.MWFM            = 0.0005;
        sim.RSAJitter       = 0.00008;
        sim.MWJitter        = 0.008;
        sim.SimulationTime  = t(end);

        [phi,valid] = synthesizer(sim,fs);
        [temp,tempT] = mkTemps(tempFile,tempJitter,fs);
        yVRS = (1 - kcard) * interp1(phi.Temp,temp(:,1),phi.Resp,'spline');
        yCRS = - kcard * interp1(phi.Temp,temp(:,2),phi.Card,'spline');

        yIBS = zeros(size(yCRS));
        switch(bltp)
            case 1
                yIBS(~initT) = ibsHeart;
            case 2
                yIBS(~initT) = ibsLung;
        end

        % Simulate yCRS, yVRS, and yAM
        maxCRS = max(yCRS);
        yAM = am * maxCRS *...
              (yCRS / max(yCRS)) .* yVRS / max(yVRS);

        % Build final signal with noise e
        yClean = (10 * maxCRS * yIBS) + yCRS + yVRS + yAM;
        e = (1 / snr) * max(yClean) * randn(size(yClean));
        yNoisy = yClean + e;
        
        % Quality assurance
        qaName = strrep("SIMULATED " + num2str(i),"_"," ");
                    
        fig = figure;
        subplot(4,1,1);
        plot(t,yClean);
        title("Clean");
        subplot(4,1,2);
        plot(t,yNoisy);
        title("Noisy");
        subplot(4,1,3);
        plot(t,yVRS);
        title("VRS");
        subplot(4,1,4);
        plot(t,yCRS + yAM);
        title("CRS + AM");

        saveas(fig,fullfile(qaFolder,qaName + ".png"));
        close all;
        
        % Save data
        clear data;
        data.t      = t;
        data.initT  = initT;
        data.yNoisy = yNoisy;
        data.yVRS   = yVRS;
        data.yCRS   = yCRS;
        data.yAM    = yAM;
        data.yClean = yClean;
        data.RR     = RR;
        data.RRmod  = valid.RR;
        data.HR     = HR;
        data.HRmod  = valid.HR;
        dataFile = fullfile(dataFolder,num2str(i));
        save(dataFile,'data');

        clear data;
        data = [t,yNoisy,yVRS,yCRS,yAM,yClean,valid.RR,valid.HR];
        dataFile = fullfile(dataFolder,num2str(i) + ".csv");
        writematrix(data,dataFile);

        str = "SIMULATION " + num2str(i) + " / " + num2str(nSims);
        disp(str);
        fprintf(fid, str + "\n");

    catch
        str = "SIMULATION " + num2str(i) + " FAILED / " + num2str(nSims);
        disp(str);
        fprintf(fid, str + "\n");
    end
end
