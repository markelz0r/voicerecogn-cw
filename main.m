path = 'audio_out_training';
files = dir(strcat(path,'\*.wav'));
L = length (files);
for f=1:L
    
    %reading audio file and resampling it
    fs_target = 16000;
    [audion_in,fs_old] = audioread(strcat(path,'\',files(f).name));
    x = audion_in(:,1); %taking chanel 1 only
    x_res = resample(x,fs_target,fs_old); %resampling onto 16khz

    %pre-emphasis-filter
    B=[1, -0.97];
    x_res = filter(B, 1, x_res, [], 2);


    %applying hamming window
    samples_num = length(x_res);
    frame_length = 320;
    frame_num = floor(samples_num/frame_length);

    i = 1;
    for frame = 1:frame_num
        if i==1 
            sample1 = (frame * frame_length) - (frame_length - 1);
            sample2 = (frame * frame_length);
        else  
            sample1 = (frame * frame_length) - (frame_length - 1) - (frame_length/2);
            sample2 = (frame * frame_length) - (frame_length/2);
        end     
        ham = hamming(frame_length);
        tf = x_res(sample1:sample2);
        [magSpec] = magAndPhase(tf);
        magSpecArr(:,i) = magSpec;
        i=i+1;
    end    

    %%Filterbank%%
    numChan = 21;
    x = floor(linspace(1, frame_length/2, numChan));

    chanMean = zeros(1, numChan-1);
    filterBankArr = zeros(numChan-1,length(magSpecArr));
    vocalTractArr = zeros(((numChan-1)/2),length(magSpecArr));
    for magSpecArrIndex = 1:length(magSpecArr)
        for arrayIndex = 2:length(x)
            magSpecElem = magSpecArr(:, magSpecArrIndex);    
            firstSamp = x(arrayIndex-1);
            lastSamp = x(arrayIndex);
            channel = magSpecElem(firstSamp:lastSamp);
            chanMean(arrayIndex-1) = mean(channel);        
        end
        filterBankArr(:,magSpecArrIndex)= chanMean;
        logOfFilterBank = log10(chanMean); 
        dctResult = dct(logOfFilterBank);
        vocalTractArr(:,magSpecArrIndex) = dctResult(1:length(dctResult)/2);
    end

    plot(vocalTractArr);



  



%     %generating a filename
%     s = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789';
% 
%     %find number of random characters to choose from
%     numRands = length(s); 
% 
%     %specify length of random string to generate
%     sLength = 10;
% 
%     %generate random string
%     randString = s( ceil(rand(1,sLength)*numRands) );

     flist = fopen('list.txt','a');
     fwrite(flist, strcat('MFCCs/train/',files(f).name)); 


    % Open file for writing:
    %fid = fopen('mfc_out/'+randString+'.mfc', 'w', 'ieee-be');
    fn = files(f).name(1:end-4);
    fid = fopen(strcat('mfc_out/',fn,'.mfc'),'w', 'ieee-be');

    vocalTractArr = vocalTractArr';
    numVectors = length(vocalTractArr);
    vectorPeriod = 625;
    numDims = 10;
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

end



function [magSpec] = magAndPhase(shortTimeFrame)
frame_length = 320;
    ham = hamming(frame_length);
    ham_res = ham.*shortTimeFrame;
    dft = fft(ham_res);
    magSpecFull = abs(dft);
    magSpec = magSpecFull(1:(length(magSpecFull)/2));  
end







