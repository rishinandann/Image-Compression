clc;
clear;
close all;

A=imread('harbour512x512.tif');
inp_img_set = {'boats512x512.tif', 'peppers512x512.tif', 'harbour512x512.tif'}; % Reading the 3 images

N=512; %image size

mse = zeros(3,10);  %Creating 3x10 array for MSE
entropy = zeros(3,10);  %Creating 3x10 array for entropy
PSNR_immse=zeros(3,10);  %Creating 3x10 array for PSNR
mse_dwt = zeros(3,10);   %Creating 3x10 array to store mse of quantized and unquantized coefficients
quantization = zeros(1,10);  
for im_idx = 1:numel(inp_img_set)
    
    Im = imread(inp_img_set{im_idx});
    %Analysis
    [LL, HL, LH, HH] = analysis_process(Im, N);     %Scale 1
    [LLLL,HLLL,LHLL,HHLL] = analysis_process(LL,N/2); %Scale 2
    [LLLLLL,HLLLLL,LHLLLL,HHLLLL] = analysis_process(LLLL,N/4); %Scale 3
    [LLLLLLLL,HLLLLLLL,LHLLLLLL,HHLLLLLL] = analysis_process(LLLLLL,N/8); %Scale 4

    %wavelet coeffs for scale 4
    a4 = [LLLLLLLL,HLLLLLLL;LHLLLLLL,HHLLLLLL];
    a3 = [a4,HLLLLL;LHLLLL,HHLLLL];
    a2 = [a3,HLLL;LHLL,HHLL];
    a1 = [a2,HL; LH, HH];  %Combining wavelet coefficients into one image
    figure;
    imshow(uint8(a1));
    title('Scale 4 Wavelet coefficients for' + string(inp_img_set(im_idx)));

    for q_lvl=0:9

        %Quantize
        LLLLLLLL_q = int32(round((LLLLLLLL)/(2^q_lvl))*(2^(q_lvl)));
        HLLLLLLL_q = int32(round((HLLLLLLL)/(2^q_lvl))*(2^(q_lvl)));
        LHLLLLLL_q = int32(round((LHLLLLLL)/(2^q_lvl))*(2^(q_lvl)));
        HHLLLLLL_q = int32(round((HHLLLLLL)/(2^q_lvl))*(2^(q_lvl)));

        HLLLLL_q = int32(round((HLLLLL)/(2^q_lvl))*(2^(q_lvl)));
        LHLLLL_q = int32(round((LHLLLL)/(2^q_lvl))*(2^(q_lvl)));
        HHLLLL_q = int32(round((HHLLLL)/(2^q_lvl))*(2^(q_lvl)));

        HLLL_q = int32(round((HLLL)/(2^q_lvl))*(2^(q_lvl)));
        LHLL_q = int32(round((LHLL)/(2^q_lvl))*(2^(q_lvl)));
        HHLL_q = int32(round((HHLL)/(2^q_lvl))*(2^(q_lvl)));

        HL_q = int32(round((HL)/(2^q_lvl))*(2^(q_lvl)));
        LH_q = int32(round((LH)/(2^q_lvl))*(2^(q_lvl)));
        HH_q = int32(round((HH)/(2^q_lvl))*(2^(q_lvl)));

        %Quantized wavelet coeffs for scale 4
        a4_q = [LLLLLLLL_q,HLLLLLLL_q;LHLLLLLL_q,HHLLLLLL_q];
        a3_q = [a4_q,HLLLLL_q;LHLLLL_q,HHLLLL_q];
        a2_q = [a3_q,HLLL_q;LHLL_q,HHLL_q];
        a1_q = [a2_q,HL_q; LH_q, HH_q];
        
       
        %Entropy
        entropy(im_idx, q_lvl+1) = entropy(im_idx, q_lvl+1) + calc_scaled_entropy(LLLLLLLL_q(:));
        entropy(im_idx, q_lvl+1) = entropy(im_idx, q_lvl+1) + calc_scaled_entropy(HLLLLLLL_q(:));
        entropy(im_idx, q_lvl+1) = entropy(im_idx, q_lvl+1) + calc_scaled_entropy(LHLLLLLL_q(:));
        entropy(im_idx, q_lvl+1) = entropy(im_idx, q_lvl+1) + calc_scaled_entropy(HHLLLLLL_q(:));

        entropy(im_idx, q_lvl+1) = entropy(im_idx, q_lvl+1) + calc_scaled_entropy(HLLLLL_q(:));
        entropy(im_idx, q_lvl+1) = entropy(im_idx, q_lvl+1) + calc_scaled_entropy(LHLLLL_q(:));
        entropy(im_idx, q_lvl+1) = entropy(im_idx, q_lvl+1) + calc_scaled_entropy(HHLLLL_q(:));

        entropy(im_idx, q_lvl+1) = entropy(im_idx, q_lvl+1) + calc_scaled_entropy(HLLL_q(:));
        entropy(im_idx, q_lvl+1) = entropy(im_idx, q_lvl+1) + calc_scaled_entropy(LHLL_q(:));
        entropy(im_idx, q_lvl+1) = entropy(im_idx, q_lvl+1) + calc_scaled_entropy(HHLL_q(:));

        entropy(im_idx, q_lvl+1) = entropy(im_idx, q_lvl+1) + calc_scaled_entropy(HL_q(:));
        entropy(im_idx, q_lvl+1) = entropy(im_idx, q_lvl+1) + calc_scaled_entropy(LH_q(:));
        entropy(im_idx, q_lvl+1) = entropy(im_idx, q_lvl+1) + calc_scaled_entropy(HH_q(:));

        entropy(im_idx, q_lvl+1) = entropy(im_idx, q_lvl+1)/(512*512);

        %Synthesis
        s4_q = synthesis_process(LLLLLLLL_q,HLLLLLLL_q,LHLLLLLL_q,HHLLLLLL_q, N/8);
        s3_q = synthesis_process(s4_q,HLLLLL_q,LHLLLL_q,HHLLLL_q,N/4);
        s2_q = synthesis_process(s3_q,HLLL_q,LHLL_q,HHLL_q,N/2);
        s1_q = synthesis_process(s2_q,HL_q,LH_q,HH_q,N);

        %displaying
        if q_lvl == 0
            figure;
        end
        subplot(2,5,q_lvl+1);
        imshow(uint8(s1_q));
        title('q = ' + string(q_lvl));

        s1_q = uint8(s1_q);
        quantization(1,q_lvl+1) = 2.^q_lvl;
        mse(im_idx, q_lvl+1) = immse(Im,s1_q);
        mse_dwt(im_idx,q_lvl+1) = immse(int32(a1),a1_q);
        PSNR_immse(im_idx,q_lvl+1) = 10*log(255^2/mse(im_idx,q_lvl+1));
    end
end

 figure;
    plot(quantization, mse(3,:),'r');
    hold on;
    plot(quantization, mse_dwt(3,:),'g');
    hold off;
    legend('MSE between original and reconstructed image','MSE between original and Quantized wavelet coefficients');
%     hold on
%     plot(entropy_avg, PSNR_mat);
%     legend('PSNR_mat')
%     hold on
%     plot(entropy_avg, PSNR_immse);
%     legend('PSNR','PSNR-mat','PSNR-immse')
    xlabel('Quantization Step Size')
    ylabel('Distortion')
    title('MSE Error comparison - FWT')

     figure;
    plot(entropy(1,:), PSNR_immse(1,:),'r');
    hold on;
    plot(entropy(2,:), PSNR_immse(2,:),'b');
    hold on;
    plot(entropy(3,:), PSNR_immse(3,:),'g');
    hold off;
    legend('Boats','Peppers','Harbour');
%     hold on
%     plot(entropy_avg, PSNR_mat);
%     legend('PSNR_mat')
%     hold on
%     plot(entropy_avg, PSNR_immse);
%     legend('PSNR','PSNR-mat','PSNR-immse')
    xlabel('bitrate')
    ylabel('PSNR')
    title('Distortion v/s bitrate')
    
    %ANALYSIS PROCESS 
function [LL, HL, LH, HH] = analysis_process(A, N)
    wname = 'db4';
    [LoD,HiD,LoR,HiR] = wfilters(wname);
    
    Ain = A;

    dwt_mat_l = zeros(N,N/2);
    dwt_mat_h = zeros(N,N/2);

    dwt_mat_ll = zeros(N/2);
    dwt_mat_lh = zeros(N/2);
    dwt_mat_hl = zeros(N/2);
    dwt_mat_hh = zeros(N/2);

    % Analysis
    for i=1:N
        row = Ain(i,:);
        [dwt_mat_l(i,:), dwt_mat_h(i,:)] = dwt_proc(row,LoD,HiD, N);
    end

    for i=1:N/2
        col = dwt_mat_l(:,i)';
        [dwt_mat_ll(:,i), dwt_mat_hl(:,i)] = dwt_proc(col,LoD,HiD, N);

        col = dwt_mat_h(:,i)';
        [dwt_mat_lh(:,i), dwt_mat_hh(:,i)] = dwt_proc(col,LoD,HiD, N);
    end
    
%     figure;
%     dwt_mat = [dwt_mat_ll, dwt_mat_hl; dwt_mat_lh, dwt_mat_hh];
%     imshow(uint8(dwt_mat))
    
    LL = dwt_mat_ll;
    HL = dwt_mat_hl;
    LH = dwt_mat_lh;
    HH = dwt_mat_hh;
end

%SYNTHESIS PROCESS
function synthesis_img = synthesis_process(LL, HL, LH, HH, N)
    % Synthesis
    wname = 'db4';
    [LoD,HiD,LoR,HiR] = wfilters(wname);
    
    dwt_mat_ll = LL;
    dwt_mat_hl = HL;
    dwt_mat_lh = LH;
    dwt_mat_hh = HH;
    
    idwt_mat_l = zeros(N,N/2);
    idwt_mat_h = zeros(N,N/2);

    idwt_mat = zeros(N);

    for i = 1:N/2
        s1 = dwt_mat_ll(:,i)';
        s2 = dwt_mat_hl(:,i)';
        idwt_mat_l(:,i) = idwt_proc(s1,s2,LoR,HiR,N)';

        s1 = dwt_mat_lh(:,i)';
        s2 = dwt_mat_hh(:,i)';
        idwt_mat_h(:,i) = idwt_proc(s1,s2,LoR,HiR,N)';
    end

    %figure;
    %imshow(uint8([idwt_mat_l,idwt_mat_h]));

    for i = 1:N
        s1 = idwt_mat_l(i,:);
        s2 = idwt_mat_h(i,:);
        idwt_mat(i,:) = idwt_proc(s1,s2,LoR,HiR,N);
    end
    
    synthesis_img = idwt_mat;
end

%DWT
function [al,dh] = dwt_proc(r,LoD,HiD, N)

shift = -7;
al = cconv(LoD,r,N);   %Circular convolution
al = circshift(al,shift); 
al = downsample(al,2);  

dh= cconv(HiD,r,N);
dh= circshift(dh,shift); % Circular Convolution
dh = downsample(dh,2);

end

%IDWT
function seq = idwt_proc(s1,s2,LoR,HiR, N)

shift = 0;
al = upsample(s1,2);
%al = filter(LoR,1,al);
al = cconv(LoR,al,N);
al = circshift(al,shift);

dh = upsample(s2,2);
%dh = filter(HiR,1,dh);
dh = cconv(HiR,dh,N);
dh= circshift(dh,shift);

seq = al + dh;

end

%FUNCTION TO CALCULATE ENTROPY
function entropy = calc_scaled_entropy(stream)
    stream = double(stream);
    [pr,symbols] = hist(stream,unique(stream));
    pr = pr/sum(pr);
    
    entropy = 0;
    for j = 1:numel(pr)
        entropy = entropy - pr(j)*log2(pr(j));
    end
    entropy = entropy*numel(stream);
end




