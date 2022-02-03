clear;
close;
clear all;

f1 = @(block_struct) dct2(block_struct.data);
f2 = @(block_struct) idct2(block_struct.data);


Im = imread('boats512x512.tif');

J = blockproc(Im,[8 8],f1);
Qi = round(J);

K = blockproc(J, [8 8], f2)/255;
imwrite(K,'Compressed image.tif');
MSE_restored = sum((double(Im(:))- K(:).^2) / numel(Im(:)));  % MSE of restored image

 
figure, imshow(Im), title('Original Image');
figure , imshow(K), title('Reconstructed Image');