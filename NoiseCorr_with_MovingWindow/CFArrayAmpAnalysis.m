% This function is use to do the amplitude analysis of the array cross-correlation functions
% by Huajian Yao  March 5, 2008  @ MIT
% modified by Huajian Yao, July 10, 2016 @ USTC
% contact: huajianyao@gmail.com

function CFArrayAmpAnalysis

global gfcn
% Define the structure of cross-correlation functions or empirical Green's functions
GreenFcnInfo = struct('Lat1',0,...
    'Lon1',0,...
    'Lat2',0,...
    'Lon2',0,...
    'StaDist',0,...
    'GreenFcn',zeros(2,1),...
    'Time',zeros(1,1),...
    'PtNum',0);
gfcn = struct(GreenFcnInfo);

%%%%%%%%%%%%%% important input parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folder = './CFs-obs/V-V/ascii/';   % folder contains all the CFs
filename = dir(strcat(folder, 'VV*'));
fs = 40;  % sampling frequency of the CFs
Ts = 0.1; % Starting period for band-pass filtering the CFs
Te = 2; % Ending period for band-pass filtering the CFs
deggrid = 20; % azimthal degree interval for estimating the amplitude of CFs, depending on how many CFs you have             
SigWinV = [0.5 4];  % group velocity window of CFs for estimating the amplitude of the signal (e.g., surface waves)
% NoiseWinV = [1.2 1.4]; % group velocity window of CFs for estimating the amplitude of noise in CFs
NoiseWinLength = 100; % length of noise window (s), which is defined right after the signal window
                     % SigWinV and NoiseWinV are highly dependent on the site and period you are looking at
SNRmin = 5; % define the minimum signal to noise ratio of CFs to be accepted for the 
            % use of estimation of azimuthal dependence of ambient noise energy
AzimCenter = 0;  % center of azimuthal angle for plotting CFs
AzimHalfBand = 90; % half band of azimth angle with respect to AzimCenter for plotting CFs
% plot CFs for station pairs with azimuthal angle within the interval (1) AzimCenter + [-AzimHalfBand, AzimHalfBand] 
% or(2) AzimCenter + 180 + [-AzimHalfBand, AzimHalfBand], the positive time part correspond to interval (1)
% and the negative time part correspond to interval (2)
AmpScaleRatio = 20; % amplitude scaling ratio for plotting CFs, this value is very empirical
ylimData = [-0 2]; % y axis plotting range
xlimData = [-400 400]; % x axis plotting range
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nfile = length(filename);

% filename = strcat(folder, filename);
[b, a] = butter(2, [(1/Te)/(fs/2) (1/Ts)/(fs/2)]);

centerdeg = 0:deggrid:(360-deggrid);
degnum = length(centerdeg);
degCFamp = zeros(degnum,1);
degCFnum = zeros(degnum,1);
degStaDist = zeros(degnum,1);

figure
kk = 0;
for mm = 1:nfile
    
    % read CFs
    RdGreenFcn(strcat(folder,filename(mm).name));
    % band-pass filtering the CFs
    gfcn.GreenFcn(1,:) = filtfilt(b,a,gfcn.GreenFcn(1,:));
    gfcn.GreenFcn(2,:) = filtfilt(b,a,gfcn.GreenFcn(2,:));
    nptCF = length(gfcn.GreenFcn(1,:));
    CFtemp = [gfcn.GreenFcn(2,nptCF:-1:2) (gfcn.GreenFcn(2,1) + gfcn.GreenFcn(1,1))/2 gfcn.GreenFcn(1,2:nptCF)];
    CFbpfilt = filtfilt(b,a,CFtemp);
    gfcn.GreenFcn(1,:) = CFbpfilt(nptCF:end);
    gfcn.GreenFcn(2,:) = CFbpfilt(nptCF:-1:1);

    SampleT = 1/fs;
    
    % get the station pair distance and azimuth angles
    gfcn.StaDist = deg2km(distance([gfcn.Lat1 gfcn.Lon1], [gfcn.Lat2 gfcn.Lon2]));
    AzimSta1Sta2 = azimuth(gfcn.Lat1, gfcn.Lon1, gfcn.Lat2, gfcn.Lon2);    
    AzimSta2Sta1 = azimuth(gfcn.Lat2, gfcn.Lon2, gfcn.Lat1, gfcn.Lon1);
    
    if gfcn.StaDist > 0
        % hilbert transform to obtain the wave envelope
        groupCF12 = abs(hilbert(gfcn.GreenFcn(1,:)));
        groupCF21 = abs(hilbert(gfcn.GreenFcn(2,:)));

        % obtain the time window for signal and noise of CFs    
        SigWinT = gfcn.StaDist./SigWinV;
        % NoiseWinT = gfcn.StaDist./NoiseWinV;

        % calculate the signal to noise ratio of CFs from station 1 --> station 2
        nnSigWin = ceil((SigWinT(2)/SampleT)):floor((SigWinT(1)/SampleT));
        nnNoiseWin = nnSigWin(end):(nnSigWin(end)+round(NoiseWinLength/SampleT));
        % nnNoiseWin = ceil((NoiseWinT(2)/SampleT)):floor((NoiseWinT(1)/SampleT));
        test = groupCF12(nnSigWin);
        SigAmp12 = max(groupCF12(nnSigWin));
        SigAmpAve12 = mean(groupCF12(nnSigWin));
        NoiseAmpAve12 = mean(groupCF12(nnNoiseWin));
        if NoiseAmpAve12 == 0
            if SigAmp12 > 0
                SNR12 = 100;
            else
                SNR12 = 0;
            end
        else
            SNR12 = SigAmp12/NoiseAmpAve12; % max signal envelope amplitue / average noise envelope amplitude
        end
        % calculate the signal to noise ratio of CFs from station 2 --> station 1
        SigAmp21 = max(groupCF21(nnSigWin));
        SigAmpAve21 = mean(groupCF21(nnSigWin));
        NoiseAmpAve21 = mean(groupCF21(nnNoiseWin));
        if NoiseAmpAve21 == 0
            if SigAmp21 > 0
                SNR21 = 100;
            else
                SNR21 = 0;
            end
        else
            SNR21 = SigAmp21/NoiseAmpAve21;  % max signal envelope amplitue / average noise envelope amplitude
        end

        % add the amplitude of CFs (two sides) to the correspondent azimuth interval
        if (SNR12 >= SNRmin || SNR21 >= SNRmin) 
            kk = kk + 1;
            ndeg12 = floor((AzimSta1Sta2 + deggrid/2)/deggrid)+1;
            if ndeg12 == degnum + 1
                ndeg12 = 1;
            end
            degCFnum(ndeg12) = degCFnum(ndeg12) + 1;
            degCFamp(ndeg12,degCFnum(ndeg12)) = SigAmp12;
            degStaDist(ndeg12,degCFnum(ndeg12)) = gfcn.StaDist;

            ndeg21 = floor((AzimSta2Sta1 + deggrid/2)/deggrid)+1;
            if ndeg21 == degnum + 1
                ndeg21 = 1;
            end
            degCFnum(ndeg21) = degCFnum(ndeg21) + 1;
            degCFamp(ndeg21,degCFnum(ndeg21)) = SigAmp21;
            degStaDist(ndeg21,degCFnum(ndeg21)) = gfcn.StaDist;
        end

        % plot CFs within certain azimuth angle range
        if (SNR12 >= SNRmin || SNR21 >= SNRmin)

            azimdev12 = mod(AzimSta1Sta2 - AzimCenter,360);
            azimdev21 = mod(AzimSta2Sta1 - AzimCenter,360);
            if azimdev12 > 180
                azimdev12 = azimdev12 - 360;
            end
            if azimdev21 > 180
                azimdev21 = azimdev21 - 360;
            end

            % amplitude normalization
            maxamp = max(max(gfcn.GreenFcn(1,:)),max(gfcn.GreenFcn(2,:)));

            % plot CFs for station pair with azimuth angle within the interval (1) AzimCenter + [-AzimHalfBand  AzimHalfBand] or 
            % (2) AzimCenter + 180 + [-AzimHalfBand  AzimHalfBand], the positive time part correspond to interval (1)
            % and the negative time part correspond to interval (2)
            if maxamp > 0
                if abs(azimdev12) <= AzimHalfBand || abs(azimdev21) <= AzimHalfBand
                    if abs(azimdev12) < AzimHalfBand 
                        plot(gfcn.Time, gfcn.StaDist+AmpScaleRatio*gfcn.GreenFcn(1,:)/maxamp,'k-', -gfcn.Time, ...
                            gfcn.StaDist+AmpScaleRatio*gfcn.GreenFcn(2,:)/maxamp,'k-', 'LineWidth', 0.5)
                        hold on
                    elseif abs(azimdev21) < AzimHalfBand
                        plot(gfcn.Time, gfcn.StaDist+AmpScaleRatio*gfcn.GreenFcn(2,:)/maxamp,'k-', -gfcn.Time, ...
                            gfcn.StaDist+AmpScaleRatio*gfcn.GreenFcn(1,:)/maxamp,'k-', 'LineWidth', 0.5)
                        hold on
                    end
                end
            end
        end
        
    end % end of if gfcn.StaDist > 0

end

kk

% plot lines of signal and noise windows
hold on
plot([0 ylimData(2)/SigWinV(1)], [0 ylimData(2)], 'r--', 'LineWidth', 1); % signal window line
hold on
plot([0 ylimData(2)/SigWinV(2)], [0 ylimData(2)], 'r--', 'LineWidth', 1); % signal window line
hold on;
plot(NoiseWinLength+[0 ylimData(2)/SigWinV(1)], [0 ylimData(2)], 'm--', 'LineWidth', 1); % noisal window line
hold on
plot(-[0 ylimData(2)/SigWinV(1)], [0 ylimData(2)], 'r--', 'LineWidth', 1);
hold on
plot(-[0 ylimData(2)/SigWinV(2)], [0 ylimData(2)], 'r--', 'LineWidth', 1);
hold on;
plot(-NoiseWinLength-[0 ylimData(2)/SigWinV(1)], [0 ylimData(2)], 'm--', 'LineWidth', 1);
hold on;
plot([xlimData(1) xlimData(2)], [0 0], 'b--', 'LineWidth', 1);

% set(gca, 'PlotBoxAspectRatio', [2 1 1]);
ylim(ylimData);
 xlim(xlimData);
titlestr = ['Period: ', num2str(Ts), '-', num2str(Te), ' s'];
title(titlestr, 'FontSize', 12);
xlabel('t (sec)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Interstation Distance (km)', 'FontSize', 16, 'FontWeight', 'bold');
% text(-140, 1.0, 'C_B_A(-t)', 'FontSize', 16, 'FontWeight', 'bold');
% text(100, 1.0, 'C_A_B(t)', 'FontSize', 16, 'FontWeight', 'bold');
box on
set(gca, 'FontSize', 16, 'FontWeight', 'bold');

% obtain the average amplitude of CFs in the signal window in each azimuth angle interval
degCFampOrig = zeros(degnum,1);
AzimCFsAmp = zeros(1, degnum);
for i = 1:degnum
    if degCFnum(i) > 0
        for j = 1:degCFnum(i)
            degCFampOrig(i,j) = degCFamp(i,j);
            degCFamp(i,j) = degCFamp(i,j)*sqrt(degStaDist(i,j));
        end
         AzimCFsAmp(i) = sum(degCFamp(i,1:degCFnum(i)))/degCFnum(i); % averaging
    end
end
AzimCFsAmp = AzimCFsAmp/max(AzimCFsAmp);

% plot the original amplitude of CFs in certain azimuth interval given by index k
% k = 35;
% figure
% plot(degStaDist(k,1:degCFnum(k)), degCFamp(k,1:degCFnum(k)),'^', 'MarkerEdgeColor','r', 'MarkerFaceColor','r', 'MarkerSize',10);
% hold on
% plot(degStaDist(k,1:degCFnum(k)), degCFampOrig(k,1:degCFnum(k)),'o', 'MarkerEdgeColor','b', 'MarkerFaceColor','b', 'MarkerSize',10);

% plot the azimuthal distribution of CF amplitude in a pie chart
figure
r = (0:1)';
r0 = ones(length(r),1);
theta = ([centerdeg 360] - deggrid/2);
theta = (-theta-90)/57.3;

% plot the azimuthal dependent image of ambient noise energy
r = (0.4:0.6:1)';

r0 = ones(length(r),1);
theta = ([centerdeg 360] - deggrid/2);
theta = (-theta-90)/57.3;
X = r*cos(theta);
Y = r*sin(theta);
C = r0*[AzimCFsAmp AzimCFsAmp(1)];


h = pcolor(X,Y,C); caxis([0 1]);
titlestr = ['Period: ', num2str(Ts), '-', num2str(Te), ' s'];
title(titlestr, 'FontSize', 12);
axis equal tight
axis off
set(h,'LineStyle','none');        
colorbar
hold on; plot(0,0,'^k', 'MarkerSize', 10, 'MarkerFaceColor', 'k');

% figure;
% for i = 1:degnum
%     thetapatch = theta(i:(i+1));
%     if AzimCFsAmp(i) >= 0
%         X = r*cos(thetapatch);
%         Y = r*sin(thetapatch);
%         C = r0*[AzimCFsAmp(i) AzimCFsAmp(i)];
%         caxis([0 1]); colormap(jet);
%         h = pcolor(X,Y,C);
%         axis equal tight
%         axis off
%         set(h,'LineStyle','none')     
%         hold on
%     end
% end
% colorbar

clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RdGreenFcn(greenfcnfile)
global gfcn
fgfcn = fopen(greenfcnfile, 'r');
gfcn.Lon1 = fscanf(fgfcn, '%f', 1);
gfcn.Lat1 = fscanf(fgfcn, '%f', 1);
temp = fgetl(fgfcn);
gfcn.Lon2 = fscanf(fgfcn, '%f', 1);
gfcn.Lat2 = fscanf(fgfcn, '%f', 1);
temp = fgetl(fgfcn);
GreenFcn = fscanf(fgfcn,'%f',[3,inf]);
gfcn.PtNum = size(GreenFcn, 2);
gfcn.Time = GreenFcn(1,:);

% % amplitude normalization
 maxamp = max(max(GreenFcn(2,:)),max(GreenFcn(3,:)));
 if maxamp > 0
     GreenFcn(2,:) = GreenFcn(2,:)/maxamp;
     GreenFcn(3,:) = GreenFcn(3,:)/maxamp;
 end

gfcn.GreenFcn = GreenFcn(2:3,:);
fclose(fgfcn);
