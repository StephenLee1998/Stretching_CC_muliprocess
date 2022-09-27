%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove instrument responses, whiten, and band-pass filter the wave

% - by Huajian Yao, 2014 Dec 23, USTC
% hjyao@ustc.edu.cn   huajianyao@gmail.com

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data Processing: remove instrument response, band-pass filtering
function bfwave = rmResp_bpfilter(seisdata, fs, freqrange, respfile, IndexRmResp, IndexWhiteSpec, SegHour)
% seisdata: data vector
% fs: sampling frequency
% freqrange = [f_low_cut f_low_pass f_high_pass f_high_cut] % filter freq band
%              f_low_cut < f_low_pass < f_high_pass < f_high_cut, e.g.,
%              [0.005 0.01 5 10] Hz
% respfile: response file name
% IndexRmResp = 1 : read RESP.* response file and remove response
%             else: do not remove instrument response
% IndexWhiteSpec: = 1, spectrum whitening; otherwise keep the original spectrum
% SegHour: data segment length (in hour) for processing the data
%          (remove instrument response, whitening, filtering,

% define instrument response information
resp = struct('Amp',0,...
    'NumZeros',0,...
    'NumPoles',0,...
    'Zeros',0,...
    'Poles',0,...
    'poly_num',0,...
    'poly_den',0);

%%
HighF = freqrange(3);  % high pass freq
LowF = freqrange(2);   % low pass freq
%% HighF must be less than fs/2, otherwise set to be 0.99*(fs/2)
if HighF > (fs/2)
    HighF = 0.99*(fs/2);
    display(['HighF error! HighF is reset to ', num2str(HighF)]);
end

%% read instrument response files
if IndexRmResp == 1 % read RESP.* response file
    [resp.Amp, resp.Numzeros, resp.Numpoles, resp.Zeros, resp.Poles] = Rd_InstruRespFile(respfile);
    [resp.poly_num, resp.poly_den] = zp2tf(resp.Zeros', resp.Poles', resp.Amp);
end

LowFMin = max(freqrange(1), 0); % lowest freq for response removal
HighFMax = min(freqrange(4), fs/2); % highest freq for reponse removal

SegLength = round(SegHour*3600*fs); % length of data parts for fft and processing (e.g., whitening)
if mod(SegLength,2) == 1 % ensure the even number of SegLength
    SegLength = SegLength + 1;
end
staNumPt = length(seisdata); % number of points in data
fftnum = ceil(staNumPt/SegLength); % number of segments for the data

%%
bfwave=zeros(1,staNumPt);

for k=0:(fftnum-1)
    %%
    if k ~=(fftnum-1)
        datapart = detrend(seisdata((1 + k*SegLength):(k + 1)*SegLength));
        fftlength=SegLength;
        fftdata=fft(datapart,fftlength);
        clear datapart
    else
        datapart = detrend(seisdata((1 + k*SegLength):staNumPt));
        fftlength = 2^(nextpow2(staNumPt - k*SegLength));
        fftdata=fft(datapart, fftlength);
        clear datapart
    end
    
    if max(abs(fftdata)) > 0 && fftlength > fs*3000 % to prevent the processing of zero or NaN data in the input data segment, more than half hour  data
        fftdata = reshape(fftdata, 1, fftlength);
        %%
        f(1:(fftlength/2+1)) = fs*(0:(fftlength/2))/fftlength;
        delta_f = fs/fftlength;
        
        
        MinFPoint = max(2, ceil(LowFMin/delta_f));
        MaxFPoint = min(fftlength/2, floor(HighFMax/delta_f));
        
        % remove instrument response
        if IndexRmResp == 1
            % remove instrument response: the first half frequency spectrum
            w(1:(fftlength/2+1)) = 2*pi*f(1:(fftlength/2+1));
            h = freqs(resp.poly_num, resp.poly_den, w); % obtain instrument response
            nn =  MinFPoint:MaxFPoint;
            % Y = XH -> X = Y/H -> X = Y*conj(H)/abs(H)^2
            h = h/max(abs(h)); % normalize the amplitude of the intrument response
            fftdata(nn) = fftdata(nn).*conj(h(nn))./(abs(h(nn)).^2 + 0.01);  % water-level deconvolution
            fftdata(1:MinFPoint) = 0;
            fftdata(MaxFPoint:(fftlength/2+1)) = 0;
            fftdata((fftlength/2+2):fftlength)=conj(fftdata((fftlength/2):-1:2)); % treat another half spectrum
        end
        
        %% spectrum whitenning (divide the smooth spectrum amplitude from running average)
        if IndexWhiteSpec == 1
            nn =  MinFPoint:MaxFPoint;
            datafamp = abs(fftdata(nn));
            winsize = max(round(0.02/delta_f), 11);
            if mod(winsize,2) == 0
                winsize = winsize + 1;
            end
            shiftpt = round((winsize+1)/2);
            datafampnew = [ones(1,winsize)*datafamp(1) datafamp ones(1,winsize)*datafamp(end)];
            datafampsmooth = filter(ones(1,winsize)/winsize,1,datafampnew);
            datafampsmooth2 = datafampsmooth((winsize+shiftpt):(winsize+shiftpt+length(nn)-1));
            KK = find(datafampsmooth2 > 0);
            JJ = isnan(datafampsmooth2);
            fftdata(nn(JJ)) = 0;
            fftdata(nn(KK)) = fftdata(nn(KK))./datafampsmooth2(KK) ;
            
            %         figure(1); hold off; subplot(2,1,1); plot(datafamp); hold on; plot(datafampsmooth2, 'r');
            %         subplot(2,1,2); plot(abs(fftdata(nn)));
        end
        
        %% band pass filtering
        LowPtN = round(LowF/delta_f);
        HighPtN = round(HighF/delta_f);
        nptdfs = round((LowF - LowFMin)/delta_f);
        if nptdfs >= 4
            nn = (LowPtN - nptdfs):(LowPtN-1);
            taperwin = hann(2*nptdfs-1)';
            fftdata(1:(LowPtN - nptdfs -1))=0;
            % figure(99); hold off; subplot(2,1,1); hold off; plot(abs(fftdata(nn)),'r');
            fftdata(nn) = taperwin(1:nptdfs).*fftdata(nn);
            % hold on; plot(abs(fftdata(nn)),'b--');
        end
        
        nptdfs = round((HighFMax - HighF)/delta_f);
        nn = (HighPtN + 1):(HighPtN + nptdfs);
        if nptdfs >= 4
            taperwin = hann(2*nptdfs-1)';
            % subplot(2,1,2); hold off; plot(abs(fftdata(nn)),'r');
            fftdata(nn)= taperwin(nptdfs:end).*fftdata(nn);
            fftdata((HighPtN + nptdfs + 1):(fftlength/2+1)) = 0;
            % hold on; plot(abs(fftdata(nn)),'b--');
        end
        
        fftdata((fftlength/2+2):fftlength)=conj(fftdata((fftlength/2):-1:2));
        
        %% time domain filtered data
        band_data_T=real(ifft(fftdata));
        
        %     figure(2); hold off; subplot(2,1,1), semilogx(f(1:(fftlength/2+1)), abs(fftdata(1:(fftlength/2+1)))); xlim([LowFMin HighFMax]);
        %     subplot(2,1,2), plot(band_data_T);
        %     waitforbuttonpress
        
        if k~=(fftnum-1)
            bfwave((1 + k*SegLength):(k + 1)*SegLength)=band_data_T;
        else
            bfwave((1 + k*SegLength):staNumPt)=band_data_T(1:(staNumPt - k*SegLength));
        end
    end
    
end