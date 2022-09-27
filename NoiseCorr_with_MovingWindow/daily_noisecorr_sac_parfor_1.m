%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parfor codes for  NoiseCorr_SAC_v4.m

% - by Huajian Yao, 2014 Dec 23, USTC
% hjyao@ustc.edu.cn   huajianyao@gmail.com

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% - modified by Qiushi Zhai, 2016 Jul 08, USTC
% zqs2010@mail.ustc.edu.cn   zjzzqs@gmail.com

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function noisecorr_sac_parfor(k,i_j,parameter_file_tmp,dayrange,daily)

%%
load(parameter_file_tmp);
i=i_j(1);
j=i_j(2);
stainfo(1) = allsta(i);
stainfo(2) = allsta(j);
display(['NO. ' num2str(k) ' of '  num2str(NumStaPair) ' station pairs : ' stainfo(1).name ' ' stainfo(2).name]);
clear i_j k;

ncf = 0;
tic

%%

MaxShiftNum = round(MaxLagTime*fsNew);
CFtime = ((-MaxShiftNum):1:MaxShiftNum)'/fsNew;
nptCF = length(CFtime);
CFdata = struct('year',0, 'day',0, 'NCF', zeros(nptCF, 1));

%% loop for year
for year= yearrange(1):yearrange(2);
    yearstr = num2str(year);
    %% loop for day
    for day = dayrange(1):dayrange(2) % loop for each julian day in a year
        
        % obtain the string for the day
        if day < 10
            daystr = ['00' num2str(day)];
        elseif day >= 10 && day < 100
            daystr = ['0' num2str(day)];
        elseif day >= 100
            daystr = num2str(day);
        end
        
        %% loop for hourxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxaAAAAAAAAAAAAAaaaaaaaaaaaaaaaaaa
        %for hour = 1:12
            
               %hourstr = num2str(hour);
            
        % specify the data file name (NOTE: USER HAS TO MODIFY THIS PART TO HAVE CORRECT INPUT DATA FILE !!!!!!!)
        % -------------------------------------------------------------------------------------------------------
        seisfile1 = [datadir allsta(i).name '/' yearstr '/' allsta(i).name '_' cmp1 '_' daystr '.sac'];
        seisfile2 = [datadir allsta(j).name '/' yearstr '/' allsta(j).name '_' cmp2 '_' daystr '.sac'];
        % -------------------------------------------------------------------------------------------------------
        %display[seisfile1];
        %display[seisfile2];
        
        % display(['read data ...' yearstr '  '  daystr]);
        % tic
        % read the specified sac data file, if this sac file does not exit,
        % the sac structure is empty (i.e., isempty(sta1) = 1)
        sta1 = readsac(seisfile1); % use this if readsacFS is not working
        sta2 = readsac(seisfile2); % use this if readsacFS is not working
        % toc
        
        % sta1 = readsacFS(seisfile1,0);
        % sta2 = readsacFS(seisfile2,0);
        
        % check whether the daily data is existing and longer than one hour
        sta1_dataok = (~isempty(sta1) && sta1.NPTS>100 *1/sta1.DELTA);
        sta2_dataok = (~isempty(sta2) && sta2.NPTS>100 *1/sta2.DELTA);
        
        
        if  sta1_dataok && sta2_dataok % sac data file for station 1 & 2 exists
            
            display(['Now processing: ' yearstr '  '  daystr '    ...... ']);
            display([sta1.DELTA,sta2.DELTA])
            %if sta1.DELTA ~= sta2.DELTA
            %    display('Error: data sampling rate is not the same!');
            %    break;
            %end
            
            % resampling sta1 waveform to a new sampling frequency
            if (1/sta1.DELTA) > fsNew
                DecimateR_org = (1/sta1.DELTA)/fsNew;
                DecimateR = round(DecimateR_org);
                if abs(DecimateR_org - DecimateR) > 0.001
                    display('Error: resampling frequency!');
                    break;
                end
                nn = floor(sta1.NPTS/DecimateR);
                if (nn*DecimateR+1) <= sta1.NPTS
                    ReSampleWave = decimate(sta1.DATA1(1:nn*DecimateR+1), DecimateR);
                else
                    ReSampleWave = decimate([sta1.DATA1(1:nn*DecimateR); sta1.DATA1(nn*DecimateR)], DecimateR);
                end
                sta1.DELTA = 1.0/fsNew;
                sta1.NPTS = length(ReSampleWave);
                sta1.DATA1 = ReSampleWave;
                clear ReSampleWave
            end
            
            % resampling sta2 waveform to a new sampling frequency
            if (1/sta2.DELTA) > fsNew
                DecimateR_org = (1/sta2.DELTA)/fsNew;
                DecimateR = round(DecimateR_org);
                if abs(DecimateR_org - DecimateR) > 0.001
                    display('Error: resampling frequency!');
                    break;
                end
                nn = floor(sta2.NPTS/DecimateR);
                if (nn*DecimateR+1) <= sta2.NPTS
                    ReSampleWave = decimate(sta2.DATA1(1:nn*DecimateR+1), DecimateR);
                else
                    ReSampleWave = decimate([sta2.DATA1(1:nn*DecimateR); sta2.DATA1(nn*DecimateR)], DecimateR);
                end
                sta2.DELTA = 1.0/fsNew;
                sta2.NPTS = length(ReSampleWave);
                sta2.DATA1 = ReSampleWave;
                clear ReSampleWave
            end
            
            % display('preprocess data:');
            % tic
            % create band pass filter
            SampleF = 1/sta1.DELTA;
            LowF = (2/SampleF)*freqrange(1);
            HighF = (2/SampleF)*freqrange(end);
            [B, A] = butter(2,[LowF, HighF]);
            
            sta1.DATA1 = detrend(sta1.DATA1 - mean(sta1.DATA1)); % demean and detrend the data
            % figure(1); subplot(2,1,1); hold off; plot(sta1.DATA1(1:1000));
            sta1.DATA1 = filtfilt(B, A, sta1.DATA1); % initially bandpass filter the data
            % figure(1); subplot(2,1,2); hold off; plot((1:1000)/SampleF,sta1.DATA1(1:1000)); waitforbuttonpress
            
            % bandpass filter the waveform and remove instrument response
            % (if IndexRmResp ~= 1, do not remove instrument response)
            bfwave1 = rmResp_bpfilter(sta1.DATA1, 1/sta1.DELTA, freqrange, allsta(i).respfile, IndexRmResp, IndexWhiteSpec,  SegHour);
            sta1.DATA1 = bfwave1;
            clear bfwave1;
            %   display(['   station 1 *** ' allsta(i).name ' ***']);
            
            sta2.DATA1 = detrend(sta2.DATA1 - mean(sta2.DATA1)); % detrend the data
            % figure(1); subplot(2,1,1); hold off; plot(sta2.DATA1(1:1000));
            sta2.DATA1 = filtfilt(B, A, sta2.DATA1); % initially bandpass filter the data
            % figure(1); subplot(2,1,2); hold off; plot((1:1000)/SampleF,sta2.DATA1(1:1000)); waitforbuttonpress
            
            % function call : read and remove instrument response; also band-pass filtering the wave trains
            bfwave2 = rmResp_bpfilter(sta2.DATA1, 1/sta2.DELTA, freqrange, allsta(j).respfile, IndexRmResp, IndexWhiteSpec, SegHour);
            sta2.DATA1 = bfwave2;
            clear bfwave2;
            %   display(['         station 2 *** ' allsta(j).name ' ***']);
            
            stadist = deg2km(distance(allsta(i).lat, allsta(i).lon, allsta(j).lat, allsta(j).lon)); % station distance (km)
            % toc
            
            % display('noise cross-correlation:');
            % tic
            % function call : single boardband noise cross-correlation with multiband time-domain normalization
            CFcnSB = CrossCorrelation(sta1, sta2, PeriodBand, MaxLagTime, indexCorrMethod);
            if sum(isnan(CFcnSB)) == 0 % no NaN data in CFcn
                ncf = ncf + 1;
                CFdata(ncf).NCF = CFcnSB;
                CFdata(ncf).year = year;
                CFdata(ncf).day = day;
                %CFdata(ncf).hour = hour;
            end
            % tocQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ
            
        end % end if: data existing
        
        %disp(ncf)
        
        %end %end for loop:hour
    end % end for loop: day
    
    
end % % end for loop: year

%% stack daily CF (normalized) of each freq band, then stack each band to form a broadband CF
%  save CFdata & CF info to mat files; save stacked broadband CF
%  data to ascii files which can be used for later dispersion
%  analysis
save('data,mat','CFdata')


if ncf > 0
    cross_info = [allsta(i).lon,allsta(i).lat, allsta(i).elev; allsta(j).lon, allsta(j).lat, allsta(j).elev];
    stackCF = zeros(length(CFtime),1);
    nnday = 0;
    for nd = 1:ncf
        dailyCFmax = max(CFdata(nd).NCF);
        if max(dailyCFmax) > 0
            stackCF1 = stackCF + CFdata(nd).NCF./dailyCFmax; % stack normalized daily CF for each period band
            stackCF = select_by_snr(stackCF1,stackCF);
            nnday = nnday + 1;
        end
    end
    disp(nnday)
    stackCF = stackCF./max(stackCF); % normalized CF in each frequency band
    
    %CFcnmatfile = [matCFdir '/' cmp1 cmp2 '_' allsta(i).name '-' allsta(j).name '.mat']; % output matCF data file name
    CFcnascfile = [asciiCFdir '/' cmp1 cmp2 '_' allsta(i).name '-' allsta(j).name '_day' num2str(daily) '_' num2str(nnday) 'd_' num2str(year) '.dat']; % ascii file name for saving broadband CF
    
    % if you don't want to save daily CF (usually very large storage
    % neededaAAAAAA
    % for saving all daily CFs for all station pairs), please comment the
    % following 3 lines with save(....).
%    save(CFcnmatfile, 'stainfo', 'yearrange', 'dayrange', 'PeriodBand', 'SampleF', 'MaxLagTime', 'cmp1', 'cmp2');    % save sta and noise CF info
%    save(CFcnmatfile, 'CFtime', 'stackCF', 'nnday', '-append');
%    save(CFcnmatfile, 'CFdata',  '-append');
    
    % save the final CF after stacking normalized daily CF
    CFcn = zeros(MaxShiftNum+1, 3);
    CFcn(:,1) = CFtime((MaxShiftNum+1):(2*MaxShiftNum+1));
    CFcn(:,2) = stackCF((MaxShiftNum+1):(2*MaxShiftNum+1));
    CFcn(:,3) = stackCF((MaxShiftNum+1):-1:1);
    save(CFcnascfile,'cross_info','CFcn','-ASCII');  % save station info and CF (broadband stack) to an ascii file
    
end
toc
clear CFdata CFtime stackCF stainfo cross_info CFcn

end
