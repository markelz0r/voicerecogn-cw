y = audioread('audio_out_training\WAV - 2018-10-16 13-24-22.wav');
x = y(:,1);

n = randn(2000,1);
n2 = 