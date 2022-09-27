%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extracting cross-correlation function (CF) from cross correlation of daily
% ambient seismic noise data (version 4.2)
% The current version deals with daily long sac format data and output daily
% cross-correlation function (in non-overlapping period bands) in mat format
% and the final stacked broadband CFs in ascii format. The input data are better
% to have the same sampling frequency.

% In version(4.1), we first remove the instrument response (using RESP files), 
% perform spectral whitening (using a running window approach), and band pass 
% filter the data in the frequency domain.
% Then we do one-bit or temporal normalization cross-correlation for 
% different bands separately. And then normalize the daily CFs in different
% bands and stack them together to form the broadband CFs. This may help to
% improve SNR of CFs in different bands than the normal one broadband prcoessing.
% Before using this code, please read the main function NoiseCorrMBand_SAC_v4
% very carefully. You have to change seisfile1 and seisfile2 in the main
% function to read data successfully.

% - by Huajian Yao, 2014 Dec 23, USTC
% hjyao@ustc.edu.cn   huajianyao@gmail.com

% Reference:
% Yao, H., van der Hilst R.D., and de Hoop, M.V., 2006. Surface-wave array
% tomography in SE Tibet from ambient seismic noise and two-station analysis
% : I - Phase velocity maps. Geophys.J. Int., Vol. 166(2), 732-744,
% doi: 10.1111/j.1365-246X.2006.03028.x.
% Yao, H., Gouedard, P., McGuire, J., Collins, J. and van der Hilst, R.D.,
% 2011. Structure of young East Pacific Rise lithosphere from ambient noise
% correlation analysis of fundamental- and higher-mode Scholte-Rayleigh waves,
% Comptes Rendues Geoscience de l'Acad��mie des Sciences, doi:10.1016/j.crte.2011.04.004

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In this version(4.2), we use multiple connecting bands for temporal normalization, 
% and stack the normalized traces to form one single broadband data, then do 
% cross-correlation. This aims to achieve more even energy distribution in
% both the time and frequency domain.
% In addition, parallel computing is available in this version. 
% - by Qiushi Zhai, 2016 Jul 08, USTC
% zqs2010@mail.ustc.edu.cn   zjzzqs@gmail.com

%------------------------------Important NOTES-----------------------------
% Before using this code, please read the main function NoiseCorrMBand_SAC_v4
% very carefully to set the parameters correctly. And you need to set 'seisfile1' 
% and 'seisfile2' in noisecorr_sac_parfor.m in order to read data successfully!
% If you don't want to save daily CF (usually very large storage needed
% for saving all daily CFs for all station pairs), please comment the
% following 3 lines (199-201) with save(....). But if you need to analyze
% daily CFs, please keep it as it is.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NoiseCorr_SAC_v4

%%---------Define the value of essential parameters in the programs---------
%__________________________________________________________________________
StaFile = './sta/sta_AXAS1_AXEC2.txt';   % file containing the station information
RespFile = 'OBSstaResp.txt'; % file containing corresponding station response file (RESP. format)
IndexRmResp = 0; % index of performing instrument response removal or not: =1 remove, will use RespFile; otherwise, RespFile is not used
datadir = '/work/li_chao/work/Axial_Seamount/Data/'; % folder containing all sac data for all stations (data must be named with Julian day ...)
PeriodBand = [0.1 15]; % multiple connecting period bands for temoporal normalization of ambient noise data
                                        % the final output CFs are broadband CFs from the smallest to largest period, here [2 40]s 
fsNew = 40; % data resampling frequency; original sample freq / fsNew must be an integer!!!

IndexWhiteSpec = 1; % = 1, spectrum whitening; otherwise use original spectrum
indexCorrMethod = 2; % = 1: one-bit cross correlation; = 2: temporal normalization cross-correlation
yearrange = [2015 2019]; % year of the data
dailyrange = [1 365]; % range of days in a year for CF calculation
dayeruption = 0 ;
MaxLagTime = 200;  % maximum time length (in seconds) for each side of CFs
SegHour = 2; % data segment length (in hour) for data preprocessing (removing response, whitenning)
cmp1 = 'Z'; % data component (Z, E, N) for the first station for the processing
cmp2 = 'Z'; % data component (Z, E, N) for the second station for the processing
outCFdir = './AXAS1_AXEC2/'; % output folder for CFs

IndexParfor=1; % =1, use parallel computing; otherwise don't use parallel computing.
n_parpool=4; % for parfor, use n_parpool cores in one CPU.(only be used in case of IndexParfor=1)
%__________________________________________________________________________
mkdir(outCFdir); % create output CF folder

%% read data file
allsta = struct('name', {}, 'net', {}, 'lat', {}, 'lon', {}, 'elev',{}, 'respfile', {});
i=0;
fstat = fopen(StaFile,'r');
while ~feof(fstat)
    name = fscanf(fstat,'%s',1);
    if ~strcmp(name,'')
        i=i+1;
        allsta(i).name = name;  %station name
        allsta(i).net = fscanf(fstat,'%s',1); %station network
        allsta(i).lat = fscanf(fstat,'%f',1); %station latitude
        allsta(i).lon = fscanf(fstat,'%f',1); %station longitude
        allsta(i).elev = fscanf(fstat,'%f',1); %station elevation
    else
        break
    end
    temp = fgetl(fstat);
end
NumSta = i;

%% read instrument file if we want to remove instrument response
if IndexRmResp == 1
    fresp = fopen(RespFile,'r');
    for i = 1:NumSta
        name = fscanf(fresp,'%s',1);
        allsta(i).respfile = fscanf(fresp,'%s',1);
    end
end


%% frequency band info
NumPB = size(PeriodBand,1); % number of period band
Tmax = max(PeriodBand(:,2)); % max period for filtering
Tmin = min(PeriodBand(:,1)); % min period for filtering
TRange = [0.75*Tmin Tmin Tmax 1.5*Tmax]; % period range for initial band pass filtering
freqrange = 1./TRange(end:-1:1); % freq range for initial band pass filtering:
% [f_low_cut f_low_pass f_high_pass f_high_cut]

%% create folders to store mat and ascii format CFs

matCFdir = [outCFdir cmp1 '-' cmp2 '/mat']; % output data folder
asciiCFdir = [outCFdir cmp1 '-' cmp2 '/ascii']; % output data folder
[s,mess,messid] = mkdir(matCFdir); % create station pair folder for saving mat CF data files
[s,mess,messid] = mkdir(asciiCFdir); % create station pair folder for saving ascii CF data files

%% cross-correlation for every combination of two stations (also include autocorrelation functions)

%prepare i j for parfor

k=0;
for i_sta=1:NumSta-1
    for j_sta=i_sta+1:NumSta
        k=k+1;
        i_j_sta(k,1)=i_sta;
        i_j_sta(k,2)=j_sta;
    end
end
clear i_sta j_sta k;
NumStaPair=length(i_j_sta(:,1));

%%% loop for year
for iyear = yearrange(1):yearrange(2);

    parameter_file_tmp=['parameter_file_tmp.mat'];
    save(parameter_file_tmp);
    
    parpool('local',n_parpool); % version > matlab2012a
    
    parfor day = dailyrange(1):dailyrange(2)
        day1 = day - 15;
        day2 = day + 15;
        daily = day
        if daily < dayeruption && day2 > dayeruption
            day2 = dayeruption;
        end
           
        if daily > dayeruption && day1 < dayeruption
            day1 = dayeruption;
        end
        
        dayrange = [day1,day2]
        
        for k=1:NumStaPair 
            daily_noisecorr_sac_parfor(k,i_j_sta(k,:),parameter_file_tmp,dayrange,daily,iyear);
        end
        
    end
    
    delete(gcp); % version > matlab2012a
    
end %%% end for loop : year
    
