%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% single broadband cross-correlation function

% - by Huajian Yao, 2014 Dec 23, USTC
% hjyao@ustc.edu.cn   huajianyao@gmail.com

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% - modified by Qiushi Zhai, 2016 Jul 08, USTC
% zqs2010@mail.ustc.edu.cn   zjzzqs@gmail.com

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CFcnSB = CrossCorrelation(sta1, sta2, PeriodBand, MaxLagTime, indexCorrMethod)
% sta1 & sta2: input sac data structure
% PeriodBand: period bands [StartT1 EndT1; StartT2 EndT2; ...] for cross-correlation
%             (first filter the data in each band, then normalization, then cross correlation)
% MaxLagTime: cross correlation time lag: between [-MaxLagTime MaxLagTime] s
% indexCorrMethod: cross-correlation method: =1 for one-bit ; =2 for temporal normalization
% CFcnSB: output CF function in a broad band

MonthDays = [31 29 31 30 31 30 31 31 30 31 30 31; 31 28 31 30 31 30 31 31 30 31 30 31];

% one bit cross correlation
SampleF = round(1/sta1.DELTA);
SampleT = sta1.DELTA;

if SampleF == 0
    error('Error! Sampling Frequency of data should not be zero!');
end

% if the first point of the two wain trains start at differnet year (e.g., 2003, 2004)
% Here we only deal with the case of one year difference
if (sta2.NZYEAR - sta1.NZYEAR) == 1
    if mod(sta1.NZYEAR, 4) == 0
        YearDays = 366;
    else
        YearDays = 365;
    end
    sta2.NZJDAY = sta2.NZJDAY + YearDays;
elseif (sta1.NZYEAR - sta2.NZYEAR) == 1
    if mod(sta2.NZYEAR, 4) == 0
        YearDays = 366;
    else
        YearDays = 365;
    end
    sta1.NZJDAY = sta1.NZJDAY + YearDays;
elseif abs(sta1.NZYEAR - sta2.NZYEAR) > 1
    display('Year Error! ');
    exit
end

DeltaTInitial = (sta2.NZJDAY - sta1.NZJDAY)*24*3600 + (sta2.NZHOUR - sta1.NZHOUR)*3600 + (sta2.NZMIN - sta1.NZMIN)*60 + ...
    (sta2.NZSEC - sta1.NZSEC) + (sta2.NZMSEC - sta1.NZMSEC)/1000;
DeltaTInitial = round(DeltaTInitial*SampleF);

% let the first point of the two wave trains start at same time (the later one)
if DeltaTInitial > 0
    if (sta1.NPTS-DeltaTInitial) > 0
        sta1.DATA1(1:(sta1.NPTS-DeltaTInitial)) = sta1.DATA1((DeltaTInitial+1):sta1.NPTS);
        sta1.DATA1((sta1.NPTS-DeltaTInitial+1):sta1.NPTS) = 0;
    else
        sta1.DATA1(1:sta1.NPTS) = 0;
    end
    DeltaTInitial = 0;
elseif DeltaTInitial < 0
    DeltaTInitial = abs(DeltaTInitial);
    if (sta2.NPTS-DeltaTInitial) > 0
        sta2.DATA1(1:(sta2.NPTS-DeltaTInitial)) = sta2.DATA1((DeltaTInitial+1):sta2.NPTS);
        sta2.DATA1((sta2.NPTS-DeltaTInitial+1):sta2.NPTS) = 0;
    else
        sta2.DATA1(1:sta2.NPTS) = 0;
    end
    DeltaTInitial = 0;
end

PointNum = min(sta2.NPTS, sta1.NPTS);
MaxShiftNum = round(MaxLagTime/SampleT);
MinShiftNum = -MaxShiftNum;
ShiftNum = MaxShiftNum - MinShiftNum + 1;  % also ShiftNum = 2*MaxTravT*SampleF + 1

NumPB = size(PeriodBand,1); % number of period band

% GFcnMB = zeros(NumPB, MaxShiftNum+1, 3); % multiband Green's function: output
% CFcnMB = zeros(NumPB, MaxShiftNum+1, 3); % multiband cross-correlation function: output

CFcnSB = zeros(ShiftNum, NumPB);

%% noise cross-correlation for multi-frequency bands
seisdata1_stack=zeros(1,PointNum);
seisdata2_stack=zeros(1,PointNum);
for np = 1:NumPB
    
    % obtain the Start and End period for the period band
    StartT = PeriodBand(np, 1);
    EndT = PeriodBand(np, 2);
    % bandpass filter and one-bit normalize the waveform
    LowF = (2/SampleF)/EndT;
    HighF = (2/SampleF)/StartT;
    [B, A] = butter(2,[LowF, HighF]);
    seisdata1 = filtfilt(B, A, sta1.DATA1);
    seisdata2 = filtfilt(B, A, sta2.DATA1);
    if indexCorrMethod == 1  % one-bit normalization
        seisdata1 = sign(seisdata1);
        seisdata2 = sign(seisdata2);
    elseif indexCorrMethod == 2 % temporal normalization by dividing the smooth amplitude (from running average)
        winsize = round(EndT*2*SampleF); % set running window length for temporal normalization for each band
        if mod(winsize,2) == 0
            winsize = winsize + 1;
        end
        shiftpt = round((winsize+1)/2);
        
        tempamp = [ones(1,winsize)*abs(seisdata1(1)) abs(seisdata1) ones(1,winsize)*abs(seisdata1(end))];
        tempamp = filter(ones(1,winsize)/winsize, 1, tempamp); % running average of amplitude
        tempamp2 = tempamp((shiftpt+winsize):(shiftpt+winsize+sta1.NPTS-1));
        KK = find(tempamp2 >0);
        JJ = isnan(tempamp2);
        seisdata1(JJ) = 0;
        %figure; hold off; subplot(2,1,1); plot(seisdata1(1:10000)); hold on; plot(tempamp((shiftpt+winsize):(shiftpt+winsize+10000-1)),'r');
        seisdata1(KK) = seisdata1(KK)./tempamp2(KK); %where the NaN will occur?
        %subplot(2,1,2); plot(seisdata1(1:10000));
        clear tempamp tempamp2 KK JJ
        
        tempamp = [ones(1,winsize)*abs(seisdata2(1)) abs(seisdata2) ones(1,winsize)*abs(seisdata2(end))];
        tempamp = filter(ones(1,winsize)/winsize, 1, tempamp); % running average of amplitude
        tempamp2 = tempamp((shiftpt+winsize):(shiftpt+winsize+sta2.NPTS-1));
        KK = find(tempamp2 >0);
        JJ = isnan(tempamp2);
        seisdata2(JJ) = 0;
        %figure; hold off; subplot(2,1,1); plot(seisdata2(1:10000)); hold on; plot(tempamp((shiftpt+winsize):(shiftpt+winsize+10000-1)),'r');
        seisdata2(KK) = seisdata2(KK)./tempamp2(KK);
        % subplot(2,1,2); plot(seisdata2(1:10000));
        clear tempamp tempamp2 KK JJ
    end
    %stack all period bands
    seisdata1_stack(1:PointNum)=seisdata1_stack(1:PointNum)+seisdata1(1:PointNum);
    seisdata2_stack(1:PointNum)=seisdata2_stack(1:PointNum)+seisdata2(1:PointNum);
end

%display('Onebit normalization:');
% one-bit cross-correlation
OneBitCross = xcorr(seisdata2_stack(1:PointNum), seisdata1_stack(1:PointNum), MaxShiftNum);

%display('Onebit Cross correlation:');

% band-pass filter the correlation function in [LowF, HighF]

% obtain the Start and End period for the period band
StartT = PeriodBand(1, 1);
EndT = PeriodBand(np, 2);
% bandpass filter and one-bit normalize the waveform
LowF = (2/SampleF)/EndT;
HighF = (2/SampleF)/StartT;
[B, A] = butter(2,[LowF, HighF]);

OneBitCross = filtfilt(B, A, OneBitCross);

%     % Differentiate the correlation function to get the Green's fcn
%     GreenFcn = zeros(ShiftNum,1);
%
%     for i = 2:(ShiftNum-1)
%         GreenFcn(i) = (OneBitCross(i+1) - OneBitCross(i-1))/2;
%     end
%     GreenFcn(1) = OneBitCross(2) - OneBitCross(1);
%     GreenFcn(ShiftNum) = OneBitCross(ShiftNum) - OneBitCross(ShiftNum - 1);
%     % if max(GreenFcn>0)
%     %     GreenFcn = GreenFcn/max(GreenFcn);
%     % end
%
%     GFcn = zeros(MaxShiftNum+1, 3); % Green's function
%     GFcn(:, 1) = (0:1:MaxShiftNum)'/SampleF;
%     GFcn(:, 2) = - GreenFcn((MaxShiftNum+1):ShiftNum);
%     GFcn(:, 3) = GreenFcn((MaxShiftNum+1):-1:1);
%
%     CFcn = zeros(MaxShiftNum+1, 3); % cross correlation function
%     CFcn(:, 1) = (0:1:MaxShiftNum)'/SampleF;
%     CFcn(:, 2) = OneBitCross((MaxShiftNum+1):ShiftNum);
%     CFcn(:, 3) = OneBitCross((MaxShiftNum+1):-1:1);

%     GFcnMB(np,:,:) = GFcn;
%     CFcnMB(np,:,:) = CFcn;

CFcnSB = OneBitCross';
CFtime = ((-MaxShiftNum):1:MaxShiftNum)'/SampleF;

end