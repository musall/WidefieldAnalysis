function [output_final, Greg_final, angle, error_curve] = Widefield_dftregistration_rotationway(buf1ft_row,buf2ft_row,usfac)

error_curve = [];
count = 1;
for rotation_angle = -25 : 0.2 : 25
    
    diagonal_length = ceil(sqrt(size(buf1ft_row,1)^2 + size(buf1ft_row,2)^2)) + 3;
    black_background = zeros(diagonal_length,diagonal_length);
    buf1ft = uint16(ImageAdd(black_background,buf1ft_row));
    buf2ft = uint16(ImageAdd(black_background,imrotate(buf2ft_row,rotation_angle)));
    
    
    buf1ft = fft2(buf1ft); %transform reference to fourier version
    buf2ft = fft2(buf2ft); %same as above
    
    % function [output Greg] = dftregistration(buf1ft,buf2ft,usfac);
    % Efficient subpixel image registration by crosscorrelation. This code
    % gives the same precision as the FFT upsampled cross correlation in a
    % small fraction of the computation time and with reduced memory
    % requirements. It obtains an initial estimate of the crosscorrelation peak
    % by an FFT and then refines the shift estimation by upsampling the DFT
    % only in a small neighborhood of that estimate by means of a
    % matrix-multiply DFT. With this procedure all the image points are used to
    % compute the upsampled crosscorrelation.
    % Manuel Guizar - Dec 13, 2007
    
    % Portions of this code were taken from code written by Ann M. Kowalczyk
    % and James R. Fienup.
    % J.R. Fienup and A.M. Kowalczyk, "Phase retrieval for a complex-valued
    % object by using a low-resolution image," J. Opt. Soc. Am. A 7, 450-458
    % (1990).
    
    % Citation for this algorithm:
    % Manuel Guizar-Sicairos, Samuel T. Thurman, and James R. Fienup,
    % "Efficient subpixel image registration algorithms," Opt. Lett. 33,
    % 156-158 (2008).
    
    % Inputs
    % buf1ft    Fourier transform of reference image,
    %           DC in (1,1)   [DO NOT FFTSHIFT]
    % buf2ft    Fourier transform of image to register,
    %           DC in (1,1) [DO NOT FFTSHIFT]
    % usfac     Upsampling factor (integer). Images will be registered to
    %           within 1/usfac of a pixel. For example usfac = 20 means the
    %           images will be registered within 1/20 of a pixel. (default = 1)
    
    % Outputs
    % output =  [error,diffphase,net_row_shift,net_col_shift]
    % error     Translation invariant normalized RMS error between f and g
    % diffphase     Global phase difference between the two images (should be
    %               zero if images are non-negative).
    % net_row_shift net_col_shift   Pixel shifts between images
    % Greg      (Optional) Fourier transform of registered version of buf2ft,
    %           the global phase difference is compensated for.  
    [Greg, output] = translation(buf1ft, buf2ft, usfac);
    
    if count == 1
        Greg_final = Greg;
        output_final = output;
        angle = rotation_angle;
        
    elseif count > 1
        if output(1) < output_final(1)
            output_final = output;
            Greg_final = Greg;
            angle = rotation_angle;
        end
    end
    
    error_curve = [error_curve, output(1)];
    count = count + 1;
end


% Then we need to fix the deviation induced by rotation. It's wired! Because
% the code do the rotation before the translation, there should not be such
% kind of error.
[Greg_final, output2] = translation(buf1ft, Greg_final, usfac);



Greg_final = uint16(abs(ifft2(Greg_final)));

Greg_final = Greg_final((size(Greg_final, 1)/2 - size(buf1ft_row, 1)/2) : (size(Greg_final, 1)/2 + size(buf1ft_row, 1)/2)-1, ...
    (size(Greg_final, 2)/2 - size(buf1ft_row, 2)/2) : (size(Greg_final, 2)/2 + size(buf1ft_row, 2)/2)-1);

plot(error_curve);
xlabel('Angle');
ylabel('Error');

end


function out=dftups(in,nor,noc,usfac,roff,coff)
% function out=dftups(in,nor,noc,usfac,roff,coff);
% Upsampled DFT by matrix multiplies, can compute an upsampled DFT in just
% a small region.
% usfac         Upsampling factor (default usfac = 1)
% [nor,noc]     Number of pixels in the output upsampled DFT, in
%               units of upsampled pixels (default = size(in))
% roff, coff    Row and column offsets, allow to shift the output array to
%               a region of interest on the DFT (default = 0)
% Recieves DC in upper left corner, image center must be in (1,1)
% Manuel Guizar - Dec 13, 2007
% Modified from dftus, by J.R. Fienup 7/31/06

% This code is intended to provide the same result as if the following
% operations were performed
%   - Embed the array "in" in an array that is usfac times larger in each
%     dimension. ifftshift to bring the center of the image to (1,1).
%   - Take the FFT of the larger array
%   - Extract an [nor, noc] region of the result. Starting with the
%     [roff+1 coff+1] element.

% It achieves this result by computing the DFT in the output array without
% the need to zeropad. Much faster and memory efficient than the
% zero-padded FFT approach if [nor noc] are much smaller than [nr*usfac nc*usfac]

[nr,nc]=size(in);
% Set defaults
if exist('roff')~=1, roff=0; end
if exist('coff')~=1, coff=0; end
if exist('usfac')~=1, usfac=1; end
if exist('noc')~=1, noc=nc; end
if exist('nor')~=1, nor=nr; end
% Compute kernels and obtain DFT by matrix products
kernc=exp((-i*2*pi/(nc*usfac))*( ifftshift([0:nc-1]).' - floor(nc/2) )*( [0:noc-1] - coff ));
kernr=exp((-i*2*pi/(nr*usfac))*( [0:nor-1].' - roff )*( ifftshift([0:nr-1]) - floor(nr/2)  ));
out=kernr*in*kernc;
return
end



function [ ima ] = ImageAdd( ima,imb )
% Add imb to the center of ima, the ima's size must larger than the imb's size

[m,n] = size(ima);
[a,b] = size(imb);

ima(floor(m/2 - a/2) : floor(m/2 + a/2)-1 , floor(n/2 - b/2) : floor(n/2 + b/2)-1) = imb;


end



function [Greg, output] = translation(buf1ft, buf2ft, usfac)      % The buf1ft and buf2ft need FFT2 before using

% Default usfac to 1
if exist('usfac')~=1, usfac=1; end

% Compute error for no pixel shift
if usfac == 0
    CCmax = sum(sum(buf1ft.*conj(buf2ft)));
    rfzero = sum(abs(buf1ft(:)).^2);
    rgzero = sum(abs(buf2ft(:)).^2);
    error = 1.0 - CCmax.*conj(CCmax)/(rgzero*rfzero);
    error = sqrt(abs(error));
    diffphase=atan2(imag(CCmax),real(CCmax));
    output=[error,diffphase];
    
    % Whole-pixel shift - Compute crosscorrelation by an IFFT and locate the
    % peak
elseif usfac == 1
    [m,n]=size(buf1ft);
    CC = ifft2(buf1ft.*conj(buf2ft));
    [max1,loc1] = max(CC);
    [max2,loc2] = max(max1);
    rloc=loc1(loc2);
    cloc=loc2;
    CCmax=CC(rloc,cloc);
    rfzero = sum(abs(buf1ft(:)).^2)/(m*n);
    rgzero = sum(abs(buf2ft(:)).^2)/(m*n);
    error = 1.0 - CCmax.*conj(CCmax)/(rgzero(1,1)*rfzero(1,1));
    error = sqrt(abs(error));
    diffphase=atan2(imag(CCmax),real(CCmax));
    md2 = fix(m/2);
    nd2 = fix(n/2);
    if rloc > md2
        row_shift = rloc - m - 1;
    else
        row_shift = rloc - 1;
    end
    
    if cloc > nd2
        col_shift = cloc - n - 1;
    else
        col_shift = cloc - 1;
    end
    output=[error,diffphase,row_shift,col_shift];
    
    % Partial-pixel shift
else
    
    % First upsample by a factor of 2 to obtain initial estimate
    % Embed Fourier data in a 2x larger array
    [m,n]=size(buf1ft);
    mlarge=m*2;
    nlarge=n*2;
    if isa(buf1ft, 'gpuArray')
        CC=zeros(mlarge,nlarge, 'gpuArray');
    else
        CC=zeros(mlarge,nlarge);
    end
    
    % Transforming the FFT into Polar Coordinate
    %function [r,g,b] = transformImage(A, Ar, Ac, Nrho, Ntheta, Method, Center, Shape)
    
    % Inputs:   A       the input image
    %           Nrho    the desired number of rows of transformed image
    %           Ntheta  the desired number of columns of transformed image
    %           Method  interpolation method (nearest,bilinear,bicubic)
    %           Center  origin of input image
    %           Shape   output size (full,valid)
    %           Class   storage class of A
    %     SizeX = size(buf1ft, 1);
    %     SizeY = size(buf1ft, 2);
    %     L1 = transformImage(buf1ft, SizeX, SizeY, SizeX, SizeY, 'nearest', size(buf1ft) / 2, 'valid');
    %     L2 = transformImage(buf2ft, SizeX, SizeY, SizeX, SizeY, 'nearest', size(buf2ft) / 2, 'valid');
    
    
    
    
    
    %     CC_polar(m+1-fix(m/2):m+1+fix((m-1)/2),n+1-fix(n/2):n+1+fix((n-1)/2)) = ...
    %         fftshift(L1).*conj(fftshift(L2));          % The polar coordinate one
    
    
    
    CC(m+1-fix(m/2):m+1+fix((m-1)/2),n+1-fix(n/2):n+1+fix((n-1)/2)) = ...
        fftshift(buf1ft).*conj(fftshift(buf2ft));   % The rectangular coordinate one
    
    % Compute crosscorrelation and locate the peak
    CC = ifft2(ifftshift(CC)); % Calculate cross-correlation
    [max1,loc1] = max(CC);
    [max2,loc2] = max(max1);
    rloc=loc1(loc2);cloc=loc2;
    CCmax=CC(rloc,cloc);
    
    % Obtain shift in original pixel grid from the position of the
    % crosscorrelation peak
    [m,n] = size(CC); md2 = fix(m/2); nd2 = fix(n/2);
    if rloc > md2
        row_shift = rloc - m - 1;
    else
        row_shift = rloc - 1;
    end
    if cloc > nd2
        col_shift = cloc - n - 1;
    else
        col_shift = cloc - 1;
    end
    row_shift=row_shift/2;
    col_shift=col_shift/2;
    
    % If upsampling > 2, then refine estimate with matrix multiply DFT
    if usfac > 2
        %%% DFT computation %%%
        % Initial shift estimate in upsampled grid
        row_shift = round(row_shift*usfac)/usfac;
        col_shift = round(col_shift*usfac)/usfac;
        dftshift = fix(ceil(usfac*1.5)/2); %% Center of output array at dftshift+1
        % Matrix multiply DFT around the current shift estimate
        CC = conj(dftups(buf2ft.*conj(buf1ft),ceil(usfac*1.5),ceil(usfac*1.5),usfac,...
            dftshift-row_shift*usfac,dftshift-col_shift*usfac))/(md2*nd2*usfac^2);
        % Locate maximum and map back to original pixel grid
        [max1,loc1] = max(CC);
        [max2,loc2] = max(max1);
        rloc = loc1(loc2); cloc = loc2;
        CCmax = CC(rloc,cloc);
        rg00 = dftups(buf1ft.*conj(buf1ft),1,1,usfac)/(md2*nd2*usfac^2);
        rf00 = dftups(buf2ft.*conj(buf2ft),1,1,usfac)/(md2*nd2*usfac^2);
        rloc = rloc - dftshift - 1;
        cloc = cloc - dftshift - 1;
        row_shift = row_shift + rloc/usfac;
        col_shift = col_shift + cloc/usfac;
        
        % If upsampling = 2, no additional pixel shift refinement
    else
        rg00 = sum(sum( buf1ft.*conj(buf1ft) ))/m/n;
        rf00 = sum(sum( buf2ft.*conj(buf2ft) ))/m/n;
    end
    error = 1.0 - CCmax.*conj(CCmax)/(rg00*rf00);
    error = sqrt(abs(error));
    diffphase=atan2(imag(CCmax),real(CCmax));
    % If its only one row or column the shift along that dimension has no
    % effect. We set to zero.
    if md2 == 1,
        row_shift = 0;
    end
    if nd2 == 1,
        col_shift = 0;
    end
    output=[error,diffphase,row_shift,col_shift];
end

% Compute registered version of buf2ft
if (nargout > 1)&&(usfac > 0),
    [nr,nc]=size(buf2ft);
    Nr = ifftshift([-fix(nr/2):ceil(nr/2)-1]);
    Nc = ifftshift([-fix(nc/2):ceil(nc/2)-1]);
    [Nc,Nr] = meshgrid(Nc,Nr);
    Greg = buf2ft.*exp(i*2*pi*(-row_shift*Nr/nr-col_shift*Nc/nc));
    Greg = Greg*exp(i*diffphase);
elseif (nargout > 1)&&(usfac == 0)
    Greg = buf2ft*exp(i*diffphase);
end
end

