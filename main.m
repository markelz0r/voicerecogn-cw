path = './audio_out_test/';
aud_voice_path = strcat(path,'*.wav');
% disp(aud_voice_path);
files = dir(aud_voice_path);
% disp(files)
L = length (files);
% disp(L)
% fopen('list.txt','w');
for f=1:L
    
    %reading audio file and resampling it
    fs_target = 16000;
    [audion_in,fs_old] = audioread(strcat(path,'/',files(f).name));
    x = audion_in(:,1); %taking chanel 1 only
    x_res = resample(x,fs_target,fs_old); %resampling onto 16khz

    rndn=randn([1 5000])/5000;
    x_res(1:5000) = rndn;
    
    %pre-emphasis-filter
     pre_f=[1, -0.97];
     x_res = filter(pre_f, 1, x_res, [], 2);


    %applying hamming window
    samples_num = length(x_res);
    frame_length = 320;
    frame_num = (floor(samples_num/frame_length))*2;
    magSpecArr = zeros(frame_num, frame_length/2);

    i = 1;
    for frame = 1:frame_num
        if i==1 
            sample1 = (frame * frame_length) - (frame_length - 1);
            sample2 = (frame * frame_length);
%         else  
%             sample1 = (frame * frame_length) - (frame_length - 1) - (frame_length/2);
%             sample2 = (frame * frame_length) - (frame_length/2);
%         end     
        else
            sample1 = oldsample2 - (frame_length/2);
            sample2 = oldsample2 + (frame_length/2 -1);
        end  
        oldsample2 = sample2;
        ham = hamming(frame_length);
        tf = x_res(sample1:sample2);
        [magSpec] = magAndPhase(tf);
        magSpecArr(i,:) = magSpec;
        i=i+1;
    end    

    % Calculating the power spectrum 
    powerSpec = zeros(frame_num, 1);
    for magSpecArrIndex = 1:length(magSpecArr)
        powerSpec(magSpecArrIndex, 1) = sum(magSpecArr(magSpecArrIndex, :).^2);
    end
    
    
    %%Mel-Scale Filterbank%%
    melLowerBound = 2595 * log10((1 + 100/700)); % Lower bound = 100Hz
    melHigherBound = 2595 * log10((1 + 8000/700)); % Lower bound = 8000Hz
    numChan = 22;   % Number of channels
    % Generating 21 mel points
    melLinSpacedArr = floor(linspace(melLowerBound, melHigherBound, numChan));

    %Converting the mel points to frequency and finding the corresponding
    %sample number in each frame
    melScaleSamp = zeros(1, numChan);
    for melLinSpacedIndex = 1:length(melLinSpacedArr)
        freqMelScl = 700 * (10^(melLinSpacedArr(melLinSpacedIndex)/2595) - 1);
        melScaleSamp(melLinSpacedIndex) = floor((frame_length+1)*freqMelScl/fs_target);
    end


% Extracting channels (triangular filterbank)
    filterbank = zeros(numChan-2, length(frame_length/2));
    for channelNumber = 2:(numChan-1)

        prevMelPoint = melScaleSamp(channelNumber-1);
        midMelPoint = melScaleSamp(channelNumber);
        nextMelPoint = melScaleSamp(channelNumber+1);

        for lastMelToMidMel = prevMelPoint:midMelPoint

            filterbank(channelNumber-1, lastMelToMidMel) = ...
            (lastMelToMidMel - prevMelPoint) / (midMelPoint - prevMelPoint);

        end

        for midMelToNextMel = midMelPoint:nextMelPoint

            filterbank(channelNumber-1, midMelToNextMel) = ...
            (nextMelPoint - midMelToNextMel) / (nextMelPoint - midMelPoint);

        end

    end


    % Applying filterbank to signal    
    filterbank = filterbank';
    filteredSignal = zeros(numChan-2, frame_length/2);
    filteredFrames= zeros(length(magSpecArr), numChan-2);
    filterMean = zeros(1, numChan-2);
    logOfFilterBank = zeros(length(magSpecArr), numChan-2);
    dctResult = zeros(length(magSpecArr), (numChan-2));
    vocalTractFrames = zeros(length(magSpecArr), numChan-1 );
    
    for magSpecArrIndex = 1:length(magSpecArr)
        for filter1 = 1:size(filterbank, 2)
            signalFrame = magSpecArr(magSpecArrIndex, :)';
            filterSignal = filterbank(:, filter1);
            filteredSignal(filter1, :) = signalFrame.*filterSignal;
            filterMean(1, filter1) = mean(filteredSignal(filter1, :));
        end
        filteredFrames(magSpecArrIndex, :) = filterMean(1, :);
        logOfFilterBank(magSpecArrIndex, :) = log10(filteredFrames(magSpecArrIndex, :)); 
        dctResult(magSpecArrIndex, :) = dct(logOfFilterBank(magSpecArrIndex, :));        
        %Truncating the result of DCT to extract vocal tract informationx    
        vocalTractFrames(magSpecArrIndex, 1:((numChan-2)/2) ) = dctResult(magSpecArrIndex, 1:size(dctResult, 2)/2);        

    end
    
    
    defCoeff = calculateDefCoef(vocalTractFrames, numChan);
    
    vocalTractFrames(:, (numChan/2):(numChan-2)) = defCoeff;
    vocalTractFrames(:, numChan-1 ) = powerSpec;

    
    %Writing the mfcc files_______________________________________

    % Open file for writing:
    %fid = fopen('mfc_out/'+randString+'.mfc', 'w', 'ieee-be');
    fn = files(f).name(1:end-4);
    fid = fopen(strcat('mfc_out/test/',fn,'.mfc'),'w', 'ieee-be');

    vocalTractArr = vocalTractFrames;
    numVectors = length(vocalTractArr);
    vectorPeriod = 0.01*10000000;
    numDims = 21;
    parmKind = 9;
    
    % Write the header information% 
    fwrite(fid, numVectors, 'int32');    % number of vectors in file (4 byte int)
    fwrite(fid, vectorPeriod, 'int32');  % sample period in 100ns units (4 byte int)
    fwrite(fid, numDims * 4, 'int16');   % number of bytes per vector (2 byte int)
    fwrite(fid, parmKind, 'int16');      % code for the sample kind (2 byte int)

    % Write the data: one coefficient at a time:
    for ii = 1: numVectors     
        for j = 1:numDims      
            fwrite(fid, vocalTractArr(ii, j), 'float32');    
        end
    end
fclose('all');
end


function [defcoef] = calculateDefCoef(mefccCoef, numChan)
    defcoef = zeros(length(mefccCoef), ((numChan-2)/2 ) );    
    for mefccCoefFrameIndex = 1:length(mefccCoef)        
        for mefccCoefFeatIndex = 1:(numChan-2)/2
            if mefccCoefFeatIndex == 1
                defcoef(mefccCoefFrameIndex, mefccCoefFeatIndex) = ...
                    mefccCoef(mefccCoefFrameIndex, mefccCoefFeatIndex + 1) - ...
                    mefccCoef(mefccCoefFrameIndex, mefccCoefFeatIndex);
            elseif mefccCoefFeatIndex == size(mefccCoef, 2)
                defcoef(mefccCoefFrameIndex, mefccCoefFeatIndex) = ...
                    mefccCoef(mefccCoefFrameIndex, mefccCoefFeatIndex) - ...
                    mefccCoef(mefccCoefFrameIndex, mefccCoefFeatIndex - 1);
            else 
                defcoef(mefccCoefFrameIndex, mefccCoefFeatIndex) = ...
                    mefccCoef(mefccCoefFrameIndex, mefccCoefFeatIndex + 1) - ...
                    mefccCoef(mefccCoefFrameIndex, mefccCoefFeatIndex - 1);
            end                        
        end        
    end
    defcoef;
end

function [magSpec] = magAndPhase(shortTimeFrame)
frame_length = 320;
    ham = hamming(frame_length);
    ham_res = ham.*shortTimeFrame;
    dft = fft(ham_res);
    magSpecFull = abs(dft);
    magSpec = magSpecFull(1:(length(magSpecFull)/2));  
end
