clc;
close all;
clear all;

M = 8;
DCT_trans = zeros(M); %DCT coefficients

i = 0;
for j=0:M-1
    DCT_trans(i+1,j+1)=sqrt(1/M) * cos((2*j +1) * i * pi / (2*M));
end

for i = 1: M -1
    for j = 0 : M-1
        DCT_trans(i+1,j+1) = sqrt(2/M) * cos((2*j +1) * i * pi / (2*M));
    end
end
