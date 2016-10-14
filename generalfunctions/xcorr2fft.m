function [shiftx,shifty,c] = xcorr2fft(image1,image2)
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   xcorr2fft computes offsets between images image1 and image2 based
    %   on Phase Correlation method. image1 & image2 are assumed to have
    %   the same dimensions.
    %   
    %   Written by: Sebastian J Streichan, EMBL, February 29, 2012
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    % %% This is an example script to use xcorr2fft. 
    % A = zeros(1800,1800);
    % B = A;
    % 
    % index = 100:300;
    % 
    % A(index,index) = magic(length(index));
    % shiftx = 5;
    % shifty = 5;
    % B(index+shiftx,index+shifty) = magic(length(index));
    % 
    % A = A  + normrnd(0,10000,size(A));
    % B = B  + normrnd(0,10000,size(B));
    % 
    % tic
    % [shiftx,shifty] = xcorr2fft(A,B)
    % toc
    % 
    % % Compare the speed to xcorr2. 
    % % tic 
    % % C = xcorr2(A,B);
    % % toc
    % %% End of example script
    
    
    F     = fftn(image1);
    Fc    = conj(fftn(image2));
    eps   = 0.00001;
    R     = F.*Fc; % DO not normalise!
    % y = ifft(x) can be inlined with y = conj(fft(conj(x)))/length(x).
    c     = ifftn(R); % Inverse fft of Phase correlation gives cross correlation map. 
    [~,i] = max(c(:));
    [I,J] = ind2sub(size(c),i);


    if abs(I-1)<abs(size(image1,1)-I+1)
       shiftx = -I+1;
    else
       shiftx =  size(image1,1)-I+1;
    end


    if abs(J-1)<abs(size(image1,2)-J+1)
        shifty = -J+1;
    else
        shifty = size(image1,2)-J+1;
    end
    
    c = c(I,J);

end