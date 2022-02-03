clc;
clear;
close all;

%inp_img_set = {'airfield512x512.tif'};
entropy_avg = zeros(3,10); %Initializing entropy average matrix 
PSNR_immse = zeros(3,10);  %Initializing PSNR immse matrix
distortion_immse = zeros(3,10);   %Initializing distortion immse matrix
inp_img_set = {'boats512x512.tif','harbour512x512.tif', 'peppers512x512.tif'};
for im_idx = 1:numel(inp_img_set)
    Im = imread(inp_img_set{im_idx});

    
    distortion = zeros(1,10);
    
    distortion_dct = zeros(1,10);
    PSNR = zeros(1,10);
  
    PSNR_mat = zeros(1,10);

    %Apply DCT 8*8 blocks
    [h ,w]=size(Im);
    for a=1:8:h
        for b=1:8:w
            DCT_8x8(a:a+7,b:b+7) = dct2(Im(a:a+7,b:b+7));  %Calling dct2 block function
        end
    end
    
    for q_lvl = 0:9
        % Quantize
        DCT_8x8_quant = int32(round((DCT_8x8)/(2^q_lvl))*(2^(q_lvl)));

        % VLC avg bitrate
        entropy = zeros(1,64);
        for pos = 0:63
            i = 1;
            stream = zeros(1,4096);
            for a=1:8:h
                for b=1:8:w
                    stream(i) = DCT_8x8_quant(a + int8(pos/8), b + mod(pos,8));
                    i = i + 1;
                end
            end
            [pr,symbols] = hist(stream,unique(stream));
            pr = pr/sum(pr);
            entropy(pos+1) = 0;
            for j = 1:numel(pr)
                entropy(pos+1) = entropy(pos+1) - pr(j)*log2(pr(j));
            end
        end
        
        entropy_avg(im_idx,q_lvl + 1) = mean(entropy);
      
       
        %Apply inverse DCT 8*8
        [h,w] = size(DCT_8x8_quant);
        for a=1:8:h
            for b=1:8:h
                IDCT_8x8(a:a+7,b:b+7) = idct2(DCT_8x8_quant(a:a+7,b:b+7));  %Calling idct2 block function
            end
        end

        %displaying
        if q_lvl == 0
            figure;
        end
        %subplot(211);imshow(Im);
        subplot(2,5,q_lvl+1);imshow(IDCT_8x8,[]);
        title('q = ' + string(q_lvl))

        %distortion(q_lvl+1) = sum(sum((Im - uint8(IDCT_8x8)).^2))/numel(Im);
        distortion_immse(im_idx,q_lvl+1) = immse(Im, uint8(IDCT_8x8));
        distortion_dct(q_lvl+1) = sum(sum((int32(DCT_8x8) - DCT_8x8_quant).^2))/numel(DCT_8x8);
        %PSNR(q_lvl+1) = 10*log(255^2/distortion(q_lvl+1));
        PSNR_immse(im_idx,q_lvl+1) = 10*log(255^2/distortion_immse(im_idx,q_lvl+1));
        %PSNR_mat(q_lvl+1) = psnr(Im, uint8(IDCT_8x8));

    end
end

   figure;
    plot(entropy_avg(1,:), PSNR_immse(1,:),'r');
    hold on;
    plot(entropy_avg(2,:), PSNR_immse(2,:),'b');
    hold on;
    plot(entropy_avg(3,:), PSNR_immse(3,:),'g');
    hold off;
   
%     hold on
%     plot(entropy_avg, PSNR_mat);
%     legend('PSNR_mat')
%     hold on
%     plot(entropy_avg, PSNR_immse);
%     legend('PSNR','PSNR-mat','PSNR-immse')
    xlabel('bitrate')
    ylabel('PSNR')
    title('Distortion v/s bitrate')