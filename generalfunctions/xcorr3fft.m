function [shift] = xcorr3fft(image1,image2)
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   xcorr2fft computes offsets between images image1 and image2 based
    %   on Phase Correlation method. image1 & image2 are assumed to have
    %   the same dimensions.
    %   
    %   Written by: Sebastian J Streichan, EMBL, February 29, 2012
    %   Extended to 3D and bug corrected by: Stefan Gunther, EMBL, March, 20, 2012
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
   % %% This is an example script to use xcorr2fft_3D. 
   %
   % ballpos=[40,30,35];
   % ballradius=5;
   % ball = fspecial3('gaussian',2*ballradius+1);
   % 
   % A = zeros(100,100,100);
   % B = A;
   % 
   % shift = [10,10,10];
   % 
   % A(ballpos(1)-ballradius:ballpos(1)+ballradius,ballpos(2)-ballradius:ballpos(2)+ballradius,ballpos(3)-ballradius:ballpos(3)+ballradius)=ball;
   % ballpos=ballpos+shift;
   % B(ballpos(1)-ballradius:ballpos(1)+ballradius,ballpos(2)-ballradius:ballpos(2)+ballradius,ballpos(3)-ballradius:ballpos(3)+ballradius)=ball;
   % 
   % %figure, imshow(max(A,[],3),[])
   % %figure, imshow(max(B,[],3),[])
   % 
   % A = A  + normrnd(0,0.0005,size(A));
   % B = B  + normrnd(0,0.0005,size(B));
   % 
   % %figure, imshow(A(:,:,35),[])
   % %figure, imshow(B(:,:,45),[])
   %
   % xcorr2fft_3D(A,B)
   %
   % %% End of example script
    
        
    F     = fftn(image1);
    Fc    = conj(fftn(image2));
%    eps   = 0.00001;


%    R     = F.*Fc./(abs(F.*Fc)+eps); % Phase correlation, Compute fft of image1 and conjugate fft of image2, elementwise multiply and normalise. 
    R     = F.*Fc; % Phase correlation, Compute fft of image1 and conjugate fft of image2, elementwise multiply and normalise. 
    c     = ifftn(R); % Inverse fft of Phase correlation gives cross correlation map. 
    [~,i] = max(c(:));
    
    [I,J,K] = ind2sub(size(c),i);


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

    if abs(K-1)<abs(size(image1,3)-K+1)
        shiftz = -K+1;
    else
        shiftz = size(image1,3)-K+1;
    end
    
    shift=[shiftx,shifty,shiftz];

end

