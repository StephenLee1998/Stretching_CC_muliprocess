% This program is to convert ascii data to sac format data.
% initial code by Song Luo, 2018/5/16

clear;

% indir = '/home/song/paper/Mianning/CFs_surface_RR';
%indir = './CFs_all/V-V/ascii';
%indir = '/work/li_chao/work/Axial_Seamount/NCFs/Distant_Station/AXEC2_AXBA1/Z-Z/ascii/'
%indir = '/work/li_chao/work/Axial_Seamount/NCFs/SmallArray/AXAS2_AXEC2/Z-Z/ascii'
indir = '/work/li_chao/work/Axial_Seamount/NCFs/Test/40days_stack/AXID1_AXBA1/Z-Z/ascii'
% indir = 'corofgn';
% indir = 'CFs';
% outdir = 'cor';
% if(~exist(outdir,'dir'))
%     mkdir(outdir)
% end

filelist = dir([indir,'/ZZ*.dat']);
% filelist = dir([indir,'/*/ZZ*.dat']);
% filelist = dir([indir,'/RR*.dat']);

for i=1:length(filelist)
    i
    AAA = load([filelist(i).folder,'/',filelist(i).name]);
    lons = AAA(1,1); lats = AAA(1,2);
    lonr = AAA(2,1); latr = AAA(2,2);
    t = AAA(3:end,1); ncfl = AAA(3:end,2); ncfr = AAA(3:end,3);
    t = [-flipud(t);t(2:end)]; ncf = [flipud(ncfl);ncfr(2:end)];
    sacfile = sachd();
    sacfile.FILENAME = [filelist(i).folder,'/',filelist(i).name,'.SAC'];
    sacfile.NPTS = length(ncf);
    sacfile.DELTA = t(2)-t(1);
    sacfile.B = t(1);
    sacfile.STLA = latr; sacfile.STLO = lonr;
    sacfile.EVLA = lats; sacfile.EVLO = lons;
    sacfile.OMARKER = 0;
    sacfile.NZYEAR = 2020; sacfile.NZJDAY = 1;
    sacfile.NZHOUR = 0; sacfile.NZMIN = 0; sacfile.NZSEC = 0; sacfile.NZMSEC = 0;
    sacfile.KSTNM = 'COR';
    sacfile.KCMPNM = 'ZZ';
    sacfile.DATA1 = ncf;
    sacfile.LCALDA = true;
    writesac(sacfile);
end

