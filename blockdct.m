clc;
clear all;
close all;
I = imread('airfield512x512.tif');
A = double(I);
%I_DCT = dct2(A);

%Apply DCT 8*8 blocks
[h ,w]=size(A);
for a=1:8:h
    for b=1:8:w
        I_8(a:a+7,b:b+7) = dct_2d(I(a:a+7,b:b+7));
    end
end

%Display DCT coefficients of Image array
% figure;
% subplot(131); imshow(I); title('Original Image');
% subplot(132); imshow(I_8); title('After 8*8');
% subplot(133); imshow(log(abs(I_8)),[]); title('log(abs)');

%Apply inverse DCT 8*8
[h,w] = size(I_8);
for a=1:8:h
    for b=1:8:h
        I_8_inv2(a:a+7,b:b+7) = idct_2d(I_8(a:a+7,b:b+7));
    end
end

%displaying

figure;
subplot(211);imshow(I),title('Original Image');
subplot(212);imshow(I_8_inv2,[]),title('Reconstructed Image');

function DCT = dct_2d(image)
image = double(image);
N = length(image);
DCT = zeros(N);
cu = sqrt(2/N);
cv = sqrt(2/N);
A = ones (1,N); A(1) = 1/sqrt(2);
for u = 0:N-1
    for v=0:N-1
        mysum = 0;
        for x=0:N-1
            for y=0:N-1
                mysum=mysum+A(x+1)*A(y+1) * cos(pi*u*(2*x+1)/(2*N))*cos(pi*v*(2*y+1)/(2*N))*image(x+1,y+1);
            end
        end
        DCT(u+1,v+1) = cu*cv*mysum;
    end
    disp([num2str(u) 'of' num2str(N-1)])
end
end

function idct = idct_2d(image)
N=length(image);
idct=zeros(N);
cu = sqrt(2/N);
cv=sqrt(2/N);
A=ones(1,N); A(1) = 1/sqrt(2);
for u =0:N-1
    for v=0:N-1
        mysum=0;
        for x=0:N-1
            for y=0:N-1
                mysum=mysum+cu*cv*A(x+1)*A(y+1) * cos(pi*u*(2*x+1)/(2*N))*cos(pi*v*(2*y+1)/2*N)*image(x+1,y+1);
            end
        end
        idct(u+1,v+1) = mysum;
    end
    disp([num2str(u) 'of' num2str((N-1))]);
end
end

        




