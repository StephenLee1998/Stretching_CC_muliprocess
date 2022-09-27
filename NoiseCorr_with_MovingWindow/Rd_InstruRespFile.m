%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read Instrument Response File

% - by Huajian Yao, 2014 Dec 23, USTC
% hjyao@ustc.edu.cn   huajianyao@gmail.com

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Normfactor, Numzeros, Numpoles, Respzero, Resppole] = Rd_InstruRespFile(RespFile)
%%
fname = fopen(RespFile,'r');
% skip 1 - 18 line
for i = 1:18
    temp1 = fgetl(fname);
end

%read line 19
temp1 = fscanf(fname,'%s',4);
Normfactor = fscanf(fname,'%f',1);
temp1 = fgetl(fname);

%read line 20
temp1 = fgetl(fname);
%read line 21
temp1 = fscanf(fname,'%s',4);
Numzeros = fscanf(fname,'%f',1);
temp1 = fgetl(fname);
%read line 22
temp1 = fscanf(fname,'%s',4);
Numpoles = fscanf(fname,'%f',1);
temp1 = fgetl(fname);
%read line 23, 24: zeroes header
temp1 = fgetl(fname);
temp1 = fgetl(fname);
%read zeros
Respzero = zeros(1, Numzeros);
for i = 1:Numzeros
    temp1 = fscanf(fname,'%s',1);
    temp = fscanf(fname,'%d',1);
    realpart = fscanf(fname,'%e',1);
    imagpart = fscanf(fname,'%e',1);
    Respzero(i) = complex(realpart, imagpart);
    temp1 = fgetl(fname);
end
%read 2 lines: poles header
temp1 = fgetl(fname);
temp1 = fgetl(fname);
%read poles
Resppole = zeros(1, Numpoles);
for i = 1:Numpoles
    temp1 = fscanf(fname,'%s',1);
    temp = fscanf(fname,'%d',1);
    realpart = fscanf(fname,'%e',1);
    imagpart = fscanf(fname,'%e',1);
    Resppole(i) = complex(realpart, imagpart);
    temp1 = fgetl(fname);
end
fclose(fname);