function samplesOut = doFFTClean(samplesIn, fftThresh, useGPU)
    %DOFFTCLEAN Actually perform the FFT clean
    nBins = 20;
    nSkipMed = 4; % skip this many samples when estimating median
    nw = 3;       % frequency neighbors to set to zero

    if fftThresh == 0 || isempty(fftThresh)
        samplesOut = samplesIn;
        return;
    end

    samplesIn = jrclust.utils.tryGpuArray(samplesIn, useGPU);

    % get the next-largest power of 2 after nSamples (for padding)
    nSamples = size(samplesIn, 1);
    nSamplesPad = 2^nextpow2(nSamples);

    samplesOut = single(samplesIn);
    sampleMeans = mean(samplesOut, 1);
    samplesOut = bsxfun(@minus, samplesOut, sampleMeans); % center the samples

    if nSamples < nSamplesPad
        samplesOut = fft(samplesOut, nSamplesPad);
    else
        samplesOut = fft(samplesOut);
    end

    % find frequency outliers
    n1 = nSamplesPad/2; % previous power of 2
    viFreq = (1:n1)';
    vrFft = mean(bsxfun(@times, abs(samplesOut(1+viFreq,:)), viFreq), 2);
    n2 = round(n1/nBins);

    for iBin = 1:nBins
        vi1 = (n2*(iBin - 1):n2*iBin) + 1;
        if iBin == nBins
            vi1(vi1 > n1) = [];
        end

        % MAD transform
        vrFftMad = vrFft(vi1);
        vrFftMad = vrFftMad - median(vrFftMad(1:nSkipMed:end));
        vrFft(vi1) = vrFftMad / median(abs(vrFftMad(1:nSkipMed:end)));
    end

    % broaden spectrum
    vrFft = jrclust.utils.tryGather(vrFft);
    vrFft = abs(vrFft);
    
    isNoise = vrFft > fftThresh;
    vi_noise = find(isNoise);
    for i_nw = 1:nw
        viA = vi_noise - i_nw;
        viA(viA < 1)=[];

        viB = vi_noise + i_nw;
        viB(viB > n1)=[];

        isNoise(viA) = 1;
        isNoise(viB) = 1;
    end

    vi_noise = find(isNoise);
    samplesOut(1 + vi_noise, :) = 0;
    samplesOut(end - vi_noise + 1, :) = 0;

    % inverse transform back to the time domain
    samplesOut = real(ifft(samplesOut, nSamplesPad, 'symmetric'));

    if nSamples < nSamplesPad
        samplesOut = samplesOut(1:nSamples,:);
    end

    samplesOut = bsxfun(@plus, samplesOut, sampleMeans); % add mean back in
    % clean up GPU memory
    [samplesIn, samplesOut, sampleMeans, vrFftMad] = jrclust.utils.tryGather(samplesIn, samplesOut, sampleMeans, vrFftMad); %#ok<ASGLU>
    samplesOut = cast(samplesOut, jrclust.utils.trueClass(samplesIn)); % cast back to the original type
end