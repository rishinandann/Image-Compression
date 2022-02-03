clc;
close all;
clear all;

wname = 'db4';
[LoD,HiD,LoR,HiR] = wfilters(wname);
r = [1 1 1 1 1 0 0 0 0 1];  %Input sequence
N = length(r);
[L,H] = dwt_proc(r,LoD,HiD,N);  %Calling DWT function
op = idwt_proc(L,H,LoR,HiR,N);  %Calling IDWT function 

function [al,dh] = dwt_proc(r,LoD,HiD, N)

shift = -7;

al = cconv(LoD,r,N);
al = circshift(al,shift);
al = downsample(al,2);


dh= cconv(HiD,r,N);
dh= circshift(dh,shift);
dh = downsample(dh,2);

end

function seq = idwt_proc(s1,s2,LoR,HiR, N)

shift = 0;
al = upsample(s1,2);
al = cconv(LoR,al,N);
al = circshift(al,shift);

dh = upsample(s2,2);
dh = cconv(HiR,dh,N);
dh= circshift(dh,shift);

seq = al + dh;
end
