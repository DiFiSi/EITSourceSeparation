addpath(genpath(".\algorithms"));
addpath(genpath(".\utils"));
addpath(genpath(".\other"));

%% Load and plot two-dimensional data
load(".\data\example\sample2D.mat");
Y = data.yNoisy;
t = data.t;
rlung = data.rlung;
llung = data.llung;
heart = data.heart;

figure; hold on;
plot(t,Y(:,rlung));
plot(t,Y(:,llung));
plot(t,Y(:,heart));
legend("Right Lung", "Left Lung", "Heart");

%% Showcase dPCA (only works for 2D data)
hyperdPCA.fs = fs;
hyperdPCA.transWin  = [];
hyperdPCA.dLpNoise 	= designfilt('lowpassiir','PassbandFrequency', 10,     ...
                                 'StopbandFrequency', 12,'PassbandRipple', 1,         ...
                                 'StopbandAttenuation', 40,'DesignMethod', 'butter',   ...
                                 'SampleRate', hyperdPCA.fs);
hyperdPCA.dLp       = designfilt('lowpassiir','PassbandFrequency', 0.6,     ...
                                 'StopbandFrequency', 0.8,'PassbandRipple', 1,         ...
                                 'StopbandAttenuation', 40,'DesignMethod', 'butter',   ...
                                 'SampleRate', hyperdPCA.fs);
hyperdPCA.dHp       = designfilt('highpassiir','PassbandFrequency', 0.8,    ...
                                 'StopbandFrequency', 0.6,'PassbandRipple', 1,         ...
                                 'StopbandAttenuation', 40,'DesignMethod', 'butter',   ...
                                 'SampleRate', hyperdPCA.fs);
hyperdPCA.maxDeltaV     = 0;       % seconds
hyperdPCA.maxDeltaC     = 1 / 3;   % seconds
hyperdPCA.dfDeltaV      = 0;       % seconds
hyperdPCA.dfDeltaC      = 1 / 3;   % seconds

out = dPCA(Y,hyperdPCA);
[yVRShat,yCRShat,tOut] = getResult(out,[rlung;llung;heart],t);

figure;
subplot(3,1,1);
plot(tOut,yVRShat + yCRShat);
subplot(3,1,2);
plot(tOut,yVRShat);
subplot(3,1,3);
plot(tOut,yCRShat);
legend("Right Lung", "Left Lung", "Heart");

%% Load and plot one-dimensional data
load(".\data\example\sample1D.mat");
nHarms = [10,6];
fs = 50;
f0VRS = data.RRmod / 60;
fVRS = f0VRS .* (1:nHarms(1));
f0CRS = data.HRmod / 60;
fCRS = f0CRS .* (1:nHarms(2));
t = data.t;
y = data.yNoisy;
yClean = data.yClean;
yCRS = data.yCRS + data.yAM;
yVRS = data.yVRS;

figure; hold on;
plot(t,y);
plot(t,yVRS);
plot(t,yCRS);
ylabel("Amplitude [a.u.]");
xlabel("Time [s]");
legend("Total","VRS","CRS");

%% Bandpass filtering
hyperBP.fs = fs;
hyperBP.dHp       = designfilt('highpassiir','PassbandFrequency', 0.8,    ...
                  'StopbandFrequency', 0.6,'PassbandRipple', 1,         ...
                  'StopbandAttenuation', 40,'DesignMethod', 'butter',   ...
                  'SampleRate', hyperBP.fs);
hyperBP.dHpDrift	= designfilt('highpassiir','PassbandFrequency', 0.2,    ...
                  'StopbandFrequency', 0.01,'PassbandRipple', 1,         ...
                  'StopbandAttenuation', 40,'DesignMethod', 'butter',   ...
                  'SampleRate', hyperBP.fs);
hyperBP.dLp       = designfilt('lowpassiir','PassbandFrequency', 0.6,     ...
                  'StopbandFrequency', 0.8,'PassbandRipple', 1,         ...
                  'StopbandAttenuation', 40,'DesignMethod', 'butter',   ...
                  'SampleRate', hyperBP.fs);
hyperBP.dLpNoise	= designfilt('lowpassiir','PassbandFrequency', 10,     ...
                  'StopbandFrequency', 12,'PassbandRipple', 1,         ...
                  'StopbandAttenuation', 40,'DesignMethod', 'butter',   ...
                  'SampleRate', hyperBP.fs);

out = bandFilt(yClean, hyperBP);

[yVRShat,yCRShat,tOut] = getResult(out,1,t);

figure;
subplot(3,1,1); hold on;
plot(t,yClean);
plot(tOut,yVRShat + yCRShat);
subplot(3,1,2); hold on;
plot(t,yVRS);
plot(tOut,yVRShat);
subplot(3,1,3); hold on;
plot(t,yCRS);
plot(tOut,yCRShat);

%% Empirical Mode Decomposition
hyperEMD.fs = fs;
hyperEMD.maxIMF = 10;

out = ceemdanFilt(yClean, hyperEMD);

[yVRShat,yCRShat,tOut] = getResult(out,1,t);

figure;
subplot(3,1,1); hold on;
plot(t,yClean);
plot(tOut,yVRShat + yCRShat);
subplot(3,1,2); hold on;
plot(t,yVRS);
plot(tOut,yVRShat);
subplot(3,1,3); hold on;
plot(t,yCRS);
plot(tOut,yCRShat);

%% Rolling Harmonic Regression
hyperHR.fs = fs;
hyperHR.nHarms = [10,5];
hyperHR.nWin = 400;

out = harmRegWin(yClean,[f0VRS,f0CRS],hyperHR);

[yVRShat,yCRShat,tOut] = getResult(out,1,t);

figure;
subplot(3,1,1); hold on;
plot(t,yClean);
plot(tOut,yVRShat + yCRShat);
subplot(3,1,2); hold on;
plot(t,yVRS);
plot(tOut,yVRShat);
subplot(3,1,3); hold on;
plot(t,yCRS);
plot(tOut,yCRShat);

%% Optimal Harmonic Filtering
hyperOpt.fs = fs;
hyperOpt.M = 150;
hyperOpt.N = 450;
hyperOpt.L = 5;

outVRS = optFilt(yClean,f0VRS,hyperOpt);
outCRS = optFilt(yClean,f0CRS,hyperOpt);

out = outVRS;
out.VRShat = outVRS.yMod;
out.CRShat = outCRS.yMod;
[yVRShat,yCRShat,tOut] = getResult(out,1,t);

figure;
subplot(3,1,1); hold on;
plot(t,yClean);
plot(tOut,yVRShat + yCRShat);
subplot(3,1,2); hold on;
plot(t,yVRS);
plot(tOut,yVRShat);
subplot(3,1,3); hold on;
plot(t,yCRS);
plot(tOut,yCRShat);
