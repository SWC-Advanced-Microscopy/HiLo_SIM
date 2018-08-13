classdef hilospeckle < handle


    properties

        eta %double
        etaSuggested = 0.0 %double

        correctShotNoise=true %boolean
        etaPreview %boolean
        useSuggestedEta=true %boolean

        targetThickness=1.0 %The "thickenss" of the optical slice to attempt

        % Set some sort of default values for the imaging parameters
        illuminationWavelength=488% in nanometers
        detectionWavelength=520 % in nanometers

        illuminationNA=0.5 
        detectionNA=0.5

        cameraGain=0.15 % ADU/electron
        readoutNoise=1.64 % RMS ADU
        magnification=1 
        pixelSize=20 % In microns


        bandPassFilterVolume %float

        showIntermediateImages=false %boolean

        % Likely not needed
        %showProcessingTime=false %boolean
        %index1 %int
        %index2 %int
        %numberImagesToProcess %int
        %stackSize %int
        %maxN % The maximum image dimension. Doesn't need to be passed around
        
        % These are the input stacks and output stacks
        uniformIm % The uniform input image
        speckleIm % The speckle input image
        HiloFinal % The output of the processing

        doGPU=true
    end

    properties (Hidden)
        verbose = false
        padded %boolean - true if images were padded
        impLo %The "Lo" image
        impHi %The "Hi" image
        differenceImage

        originalWidth % Original image width before padding
        originalHeight % Original image height before padding

    end


    methods
        function obj=hilospeckle(uniformIm, speckleIm)
            if ~isequal(size(uniformIm), size(speckleIm))
                fprintf('uniformIm and speckleIm must be the same size. Not running.\n')
                return
            end

            obj.uniformIm = uniformIm;
            obj.speckleIm = speckleIm;
            obj.runHiLo;
            
            obj.showResults

        end %hilospeckle


        function showResults(obj)
            if size(obj.uniformIm,3)>1
                fprintf('Displaying the first frame only\n')
            end
            ii=1;
            clf
            subplot(2,2,1),imagesc(obj.uniformIm(:,:,ii))
            subplot(2,2,2),imagesc(obj.speckleIm(:,:,ii))
            subplot(2,2,3),imagesc(obj.differenceImage(:,:,ii))
            subplot(2,2,4),imagesc(obj.HiloFinal(:,:,ii))
            colormap gray
        end

        function OUT=doFFT(obj, thisImage)
            %Supposed to return an object of class FHT
            paddedImage = obj.padImage(thisImage);

            % Run an FFT
            OUT = fftshift(fft2(paddedImage));
        end %doFFT


        function OUT=doIFFT(obj, thisImage)
            %Supposed to return an object of class FHT
            %paddedImage = obj.padImage(thisImage);
            % Run the inverse FFT
            OUT = ifft2(ifftshift(thisImage));
            if obj.doGPU
                OUT = gather(OUT);
            end
        end %doIFFT

        
        function out=padImage(obj,ip)
            % Only called by doFFT
            obj.originalWidth = size(ip,2);
            obj.originalHeight = size(ip,1);
            maxN = max([obj.originalWidth, obj.originalHeight]);
            
            out=ip;
            padded=false; %DO NOT PAD FOR NOW
            return
            ii = 2;
            while ii < maxN
                ii = ii*2;
            end

            % Don't pad if not needed
            if ii == maxN && obj.originalWidth == obj.originalHeight
                obj.padded = false; %only used by doFFT
                out=ip;
                return
            end

            maxN = ii;

            % pad bottom and right edge with the mean value of the image (TODO: is that what MATLAB wants for FFT?)
            out = padarray(ip,[maxN,maxN]-size(ip), mean(ip(:)), 'post');
            obj.padded = true; %only used by doFFT
        end %close pad


        function gaussianFilter = createGaussianFilter(obj,sigmaReal, sizeX, sizeY)
            % Tested and produces something that looks plausible
            % Return the FFT of a Gaussian
            thisGauss = zeros(sizeY,sizeX);

            %Not idiomatic
            for x = 1:sizeX
                for y = 1:sizeY
                    thisGauss(x,y) = (x + 0.5 - (sizeX / 2)) * (x + 0.5 - (sizeX / 2)) + (y + 0.5 - (sizeY / 2)) * (y + 0.5 - (sizeY / 2));
                end
            end

            thisGauss = -thisGauss / (2*sigmaReal^2);
            thisGauss = exp(thisGauss);
            thisGauss = thisGauss/sum(thisGauss(:));

            if obj.verbose
                fprintf('Built Gaussian kernel of size %d by %d\n', size(thisGauss))
            end

            % Now do the FFT of the gaussian 
            gaussianFilter = obj.doFFT(thisGauss);
        end



        function bandPassFilterPixels = createBandPassFilterFourier(~,kBigSigma, ratio, sizeX, sizeY)
            % Tested and produces something that looks plausible
            % Create a 2-D bandpass filter
            bigFilter_k   = -1.0 * (2*pi / sizeX) * (2*pi / sizeY) * kBigSigma^2;
            smallFilter_k = bigFilter_k * ratio^2;
            radial2Pixels = zeros(sizeY, sizeX);

            for x = 1:sizeX
                for y = 1:sizeY
                    radial2Pixels(y,x) = ... 
                      (x + 0.5 - sizeX/2) * (x + 0.5 - sizeX/2) + (y + 0.5 - sizeY/2) * (y + 0.5 - sizeY/2);
                 end
            end

            bandPassFilterPixels = exp(bigFilter_k * radial2Pixels) - ...
                                   exp(smallFilter_k * radial2Pixels);

        end % createBandPassFilterFourier



        function pixel = createCameraOTF(~, imWidth, imHeight)
            % Tested and produces something that looks plausible
            pixel = ones(imHeight, imWidth);

            for x = 1:imWidth
                for y = 1:imHeight

                    x1 = (2 * x - imWidth) * pi / imWidth;
                    y1 = (2 * y - imHeight) * pi / imHeight;

                    if x1 ~= 0 && y1 ~= 0
                        pixel(y,x) = sin(x1) * sin(y1) / x1 / y1;

                    elseif x1 == 0 && y1 ~= 0
                        pixel(y,x) = sin(y1) / y1;

                    elseif x1 ~= 0 && y1 == 0
                        pixel(y,x) = sin(x1) / x1;
                    %otherwise it's a one
                    end

               end % for y
            end % for x

            pixel(pixel<0) = 0;

        end % createCameraOFT


        function pixel = createImagingOTF(obj, numericalAperture, wavelengthInNM, sizeX,  sizeY)
            % Tested but tends to yield zero arrays
            BandWidth = 2 * numericalAperture / (wavelengthInNM * 1E-9);
            scaleUnits = obj.magnification / (obj.pixelSize * 1E-6) / BandWidth;

            pixel = zeros(sizeY,sizeX);

            for x = 1:sizeX
                for y = 1:sizeY
                    pixel(y,x) = sqrt((x + 0.5 - sizeX / 2) * (x + 0.5 - sizeX / 2) / (sizeX / 2 - 1) / (sizeX / 2 - 1) + (y + 0.5 - sizeY / 2) * (y + 0.5 - sizeY / 2) / (sizeY / 2 - 1) / (sizeY / 2 - 1)) * scaleUnits;
                end
            end
            
            pixel(pixel>1) = 1;
            pixel = 0.6366197723675814 * acos(pixel(y,x)) - pixel(y,x) * sqrt(1 - pixel.^2);
            pixel(pixel<0) = 0;

            if all(pixel(:)==0)
                fprintf('output of createImagingOTF is all zeros\n')
            end

        end % createImagingOTF


        function eta = findEta(~, cameraOTFFilter, detectionOTF, illuminationOTF, bandPassFilterFourier)
            % Not tested. No idea what the eta is.
            denominator = (bandPassFilterFourier .* detectionOTF .* cameraOTFFilter).^2 * abs(illuminationOTF);
            denominator = sum(denominator(:));
            numerator = sum(illuminationOTF(:));
            fudgeFactor = 1.2;
            eta = sqrt(numerator / denominator) * fudgeFactor;
        end %findEta


        function volume = findFilterVolume(~, thisFilter1, thisFilter2)
            % Not tested
            volume = thisFilter1.^2 .* thisFilter2.^2;  % Should be correct 
            volume = sum(volume(:)) / numel(thisFilter1); % Should be correct 
        end

        function filteredImage = applyFilter(obj, ftTarget, ftFilter)
            % Not tested, but ought to be correct (assuming I remember something about FFT)

            % Looks like this multiplies images in the fourier domain then does an iFFT then crops
            % Inputs were of class FHT and output was class ImagePlus

            if obj.doGPU
                A=gpuArray(ftTarget);
                B=gpuArray(ftFilter);
            else
                A=ftTarget;
                B=ftFilter;
            end

            filteredImage = A .* B; %ftTarget.multiply(ftFilter);
            filteredImage = obj.doIFFT(filteredImage);

            % filtered.swapQuadrants(); MAY NEED TO DO THIS

            % Return the original image without the padded extra regions
            filteredImage = filteredImage(1:obj.originalHeight, 1:obj.originalWidth);
            filteredImage = fftshift(filteredImage);
        end


    end

end
