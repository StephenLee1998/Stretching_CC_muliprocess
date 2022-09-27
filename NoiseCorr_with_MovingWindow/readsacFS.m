%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function HdrData =readsacFS(filename,plotornot)
% sacdata = readsacFS(filename,plotornot)
% Reads in data saved from within SAC and plots or not (plotornot = 1: plot)
% HdrData contains information about the earthquake and receiver, and data
% modified by FJS April 23th 1998
% last modified by Huajian Yao 2013, 2012

ppath=matlabpath;
    fid=fopen(filename,'rb');
    
if fid==-1
  % display([ 'File ',filename,' does not exist in current path ',pwd]);
  HdrData = [];
end

if ~isempty(HdrData)
    HdrFloats=fread(fid,70,'float32');
    HdrNhdr=fread(fid,15,'int32');
    HdrIhdr=fread(fid,20,'int32');
    HdrLhdr=fread(fid,5,'int32');
    HeaderStrings=str2mat(fread(fid,[8 24],'char'))';
    SeisData=fread(fid,HdrNhdr(10),'float32');
    fclose(fid);

    HdrData=struct(...
      'NPTS',HdrNhdr(10),...                    % number of points
      'DELTA',HdrFloats(1),...                  % sampling time
      'SCALE',HdrFloats(4),...                  % amplitude scaling factor
      'B',HdrFloats(6),...                      % begin time of record
      'E',HdrFloats(7),...                      % end time of record
      'O',HdrFloats(8),...                      % event origin time (seconds relative to reference recording time)
      'NZYEAR',HdrNhdr(1),...                   % year
      'NZJDAY',HdrNhdr(2),...                   % julian day
      'NZHOUR',HdrNhdr(3),...                   % hour
      'NZMIN',HdrNhdr(4),...                    % minute
      'NZSEC',HdrNhdr(5),...                    % second
      'NZMSEC',HdrNhdr(6),...                   % milisecond
      'KSTNM',deblank(HeaderStrings(1,:)),...   % station name
      'KCMPNM',deblank(HeaderStrings(21,:)),... % recording component
      'KNETWK',deblank(HeaderStrings(22,:)),... % station network      
      'KINST',deblank(HeaderStrings(24,:)),...  % generic name of recording instrument  
      'STLA',HdrFloats(32),...                  % station latitude
      'STLO',HdrFloats(33),...                  % station longitude
      'STEL',HdrFloats(34),...                  % station elevation
      'EVLA',HdrFloats(36),...                  % event latitude
      'EVLO',HdrFloats(37),...                  % event longitude
      'EVDP',HdrFloats(39),...                  % event depth
      'DIST',HdrFloats(51),...                  % epicentral distance in km
      'AZ',HdrFloats(52),...                    % azimuth
      'BAZ',HdrFloats(53),...                   % back azimuth
      'GCARC',HdrFloats(54),...                 % epicentral distance in degrees between source and receiver
      'DATA1',SeisData);                        % waveform data
    if plotornot==1
      plot(linspace(HdrData.B,HdrData.E,HdrData.NPTS),HdrData.DATA1);  
      title([filename]);
      xlabel([ 'Time (s)']);
    end
end