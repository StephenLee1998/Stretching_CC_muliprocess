README for ambient noise cross-correlation codes (Version 4.2)

Huajian Yao (huajianyao@gmail.com)

NoiseCorr_SAC_v4.m  is the main code for ambient noise cross-correlation.
CFArrayAmpAnalysis.m is used to plot the ascii format cross-correlation functions.
All other .m files are the subroutines of NoiseCorr_SAC_v4.m

The current version deals with daily long sac format data and output daily
cross-correlation function in mat format and the final stacked broadband CFs 
in ascii format. The input data need to have the same sampling frequency.

We first remove the instrument response (using RESP files),  perform spectral 
whitening (using a running window approach), and band pass  filter the data 
in the frequency domain. Then we use multiple connecting bands for temporal 
normalization, and stack the normalized traces of each to form one single 
broadband data, then do cross-correlation. This aims to achieve more even energy distribution in both time and frequency domain.In addition, parallel computing is available in this version.

%------------------------------Important NOTES-----------------------------
Before using this code, please read the main function NoiseCorrMBand_SAC_v4
very carefully to set the parameters correctly. And you need to set 'seisfile1' 
and 'seisfile2' in noisecorr_sac_parfor.m in order to read data successfully!
If you don't want to save daily CF (usually very large storage needed
for saving all daily CFs for all station pairs), please comment the
following 3 lines (199-201) with save(....). But if you need to analyze
daily CFs, please keep it as it is.

Reference:
Yao, H., van der Hilst R.D., and de Hoop, M.V., 2006. Surface-wave array
tomography in SE Tibet from ambient seismic noise and two-station analysis
: I - Phase velocity maps. Geophys.J. Int., Vol. 166(2), 732-744,
doi: 10.1111/j.1365-246X.2006.03028.x.
Yao, H., Gouedard, P., McGuire, J., Collins, J. and van der Hilst, R.D.,
2011. Structure of young East Pacific Rise lithosphere from ambient noise
correlation analysis of fundamental- and higher-mode Scholte-Rayleigh waves,
Comptes Rendues Geoscience de l'Acad®¶mie des Sciences, doi:10.1016/j.crte.2011.04.004
