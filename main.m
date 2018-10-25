



%reading audio file and resampling it
fs_target = 16000;
[audion_in,fs_old] = audioread('audio_out_training\WAV - 2018-10-16 13-24-22.wav');
x = y(:,1); %taking chanel 1 only
x_res = resample(x,fs_target,fs_old); %resampling onto 16khz


%applying hamming window
frame_length = 320;
samples_num = length(x_res);
frame_num = floor(samples_num/frame_length);
for frame = 1:frame_num
    sample1 = (frame * frame_length) - (frame_length - 1);
    sample2 = (frame * frame_length);
    ham = hamming(frame_length);
    tf = x_res(sample1:sample2);
    ham_res = ham.*tf;
    dft = fft(ham_res);
    plot(dft, 'r');
end    






