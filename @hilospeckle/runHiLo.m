function runHiLo(obj)

        if ~isequal(size(obj.uniformIm), size(obj.speckleIm))
            fprintf('uniformIm and speckleIm must be the same size. Not running.\n')
            return
        end

        if isequal(obj.uniformIm,obj.speckleIm)
            fprintf('Uniform and speckleIm are the same. Taking shortcuts\n')
            doShortcut=true;
        else
            doShortcut=false;
        end

        imageWidth = size(obj.uniformIm,2);
        imageHeight = size(obj.uniformIm,1);

        waveletGaussiansRatio = 2;
        sigmaW = obj.targetThickness;
        sigmaLo = 2.86 * sigmaW * (1 + waveletGaussiansRatio / 2);

        tic
        fprintf('Creating filters... ');

        maxN = max(size(obj.uniformIm));

        % TODO - nasty this is done twice
        paddedWidth = 2;

        while paddedWidth < maxN
            paddedWidth = paddedWidth * 2;
        end


        %maxN = paddedWidth;
        maxN = imageHeight; %DO NOT PAD FOR NOW
        
        % Create filters
        gaussianLP = obj.createGaussianFilter(sigmaLo, maxN, maxN);
        gaussianBPNarrow = obj.createGaussianFilter(2 * sigmaW, maxN, maxN);
        gaussianBPWide = obj.createGaussianFilter(2 * waveletGaussiansRatio * sigmaW, maxN, maxN);

        cameraOTF = obj.createCameraOTF(imageWidth, imageHeight);
        bandPassFilterFourier = obj.createBandPassFilterFourier(sigmaW, waveletGaussiansRatio, imageWidth, imageHeight);
        illuminationOTF = obj.createImagingOTF(obj.illuminationNA, obj.illuminationWavelength, imageWidth, imageHeight);
        detectionOTF = obj.createImagingOTF(obj.detectionNA, obj.detectionWavelength, imageWidth, imageHeight);
        suggestedEta = obj.findEta(cameraOTF, detectionOTF, illuminationOTF, bandPassFilterFourier);

        if obj.correctShotNoise
            obj.bandPassFilterVolume = obj.findFilterVolume(cameraOTF, bandPassFilterFourier);
        end


        obj.HiloFinal = zeros(size(obj.uniformIm));
        fprintf(' completed in %0.1f s\n', toc)
        tic
        for ii = 1:size(obj.uniformIm,3) %Loop through the stack

            fprintf('Processing %d/%d\n', ii, size(obj.uniformIm,3))

            % Pull out this frame and make them floats
            ipU = single(obj.uniformIm(:,:,ii));
            ipS = single(obj.speckleIm(:,:,ii));

            impUNoNorm = ipU; % Uniform temporary image not normalized
            impSNoNorm = ipS; % Structured temporary image not normalized
            
            impULocalMean = ipU; % Local mean uniform
            impULocalMean = obj.applyFilter(obj.doFFT(impULocalMean), gaussianLP);
            
            impULocalMean(impULocalMean<0.001) = 0.001;

            if doShortcut
                impSLocalMean = impULocalMean;
            else
                impSLocalMean = ipS; % Local mean structured
                impSLocalMean = obj.applyFilter(obj.doFFT(impSLocalMean), gaussianLP);
                impSLocalMean(impSLocalMean<0.001) = 0.001;
            end


            impU = ipU./impULocalMean;

            if doShortcut
                impS = impU;
            else
                impS = ipS./impSLocalMean;
            end


            impDummy1 = impU - impS;

            obj.differenceImage=real(impDummy1);

            if doShortcut
                obj.impLo = zeros(size(impDummy1));
            else
                obj.impLo = obj.applyFilter(obj.doFFT(impDummy1), gaussianBPNarrow);
                impDummy1 = obj.applyFilter(obj.doFFT(impDummy1), gaussianBPWide);
                obj.impLo = abs(obj.impLo-impDummy1); 

            end


            if obj.correctShotNoise
               readoutVariance = obj.readoutNoise^2;

               impDummy1 = (impUNoNorm*obj.cameraGain) + readoutVariance;
               impDummy1 = impDummy1 ./ impULocalMean;
               impDummy1 = impDummy1 ./ impULocalMean;


               impDummy2 = (impSNoNorm*obj.cameraGain) + readoutVariance;
               impDummy2 = impDummy2 ./ impSLocalMean;
               impDummy2 = impDummy2 ./ impSLocalMean;

               impShotNoise = impDummy1 + impDummy2;
               impShotNoise = impShotNoise * obj.bandPassFilterVolume;

               obj.impLo = obj.impLo.^2 - impShotNoise;

               impDummy1 = obj.impLo;
               impDummy1(impDummy1<0)=0;

               impDummy2 = obj.impLo;
               impDummy2(impDummy2>0)=0;
               impDummy2 = abs(impDummy2);

               obj.impLo = sqrt(impDummy1) - sqrt(impDummy2);
            end


            obj.impLo = obj.impLo .* impUNoNorm;
            obj.impLo = obj.applyFilter(obj.doFFT(obj.impLo), gaussianLP);
            impDummy1 = obj.applyFilter(obj.doFFT(impUNoNorm), gaussianLP);

            obj.impHi = impUNoNorm - impDummy1;



            if obj.useSuggestedEta
                obj.eta = suggestedEta;
            end

            obj.impLo = obj.impLo .* obj.eta;
            obj.HiloFinal(:,:,ii) = obj.impLo + obj.impHi;
            
        end %for ii = 1:size(obj.uniformIm,3)
        obj.HiloFinal(obj.HiloFinal<0)=0;
        obj.HiloFinal = uint16(obj.HiloFinal/max(obj.HiloFinal(:))*2^16);
        fprintf('\nCompleted in %0.1f s\n', toc)
end % hilo