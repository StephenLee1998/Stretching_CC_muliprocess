%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the signal-noise ratio to selsct the daily cross-correlation
% function for stack  -by li_chao 2021.11.20 Friday NJU 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [CFs_stack] = Select_by_snr(StackCF,CFdata) 
StackCF1=StackCF ;
StackCF2=StackCF + CFdata; 
s1=max(abs(StackCF1));s2=max(abs(StackCF2));
n1=rms(s1);n2=rms(s2);
snr1=s1/n1;snr2=s2/n2;
if snr1<snr2 
    CFs_stack = StackCF2;
else
    CFs_stack = StackCF1;
end
end