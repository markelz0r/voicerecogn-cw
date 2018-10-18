y1 = audioread('audio_out_training\WAV - 2018-10-16 13-24-22.wav');
x = y1(:,1);

n = randn(2000,1);
n2 = n / 10^4;

x(1:2000) = n2;


audiowrite('audio_out_training_m\WAV - 2018-10-16 13-24-22.wav',y,Fs);