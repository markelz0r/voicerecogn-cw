



%reading audio file and resampling it
fs_target = 16000;
[audion_in,fs_old] = audioread('audio_out_training\WAV - 2018-10-16 13-24-22.wav');
x = audion_in(:,1); %taking chanel 1 only
x_res = resample(x,fs_target,fs_old); %resampling onto 16khz

%pre-emphasis-filter
B=[1, -0.97];
x_res = filter(B, 1, x_res, [], 2);


%applying hamming window
samples_num = length(x_res);
frame_length = 320;
frame_num = floor(samples_num/frame_length);
magSpecArr;
i = 1;

for frame = 1:frame_num
    sample1 = (frame * frame_length) - (frame_length - 1);
    sample2 = (frame * frame_length);
    ham = hamming(frame_length);
    tf = x_res(sample1:sample2);
    [magSpec] = magAndPhase(tf);
    plot(magSpec);
    magSpecArr(:,i) = magSpec;
    i=i+1;
end    

plot(magSpecArr);
function [magSpec] = magAndPhase(shortTimeFrame)
frame_length = 320;
    ham = hamming(frame_length);
    ham_res = ham.*shortTimeFrame;
    dft = fft(ham_res);
    magSpecFull = abs(dft);
    magSpec = magSpecFull(1:(length(magSpecFull)/2));  
end







