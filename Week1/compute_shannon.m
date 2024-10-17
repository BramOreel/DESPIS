
fs = 16000; %(sampling rate)
duration = 1;
N = duration*fs;
windowS = 1024;
Noverlap = windowS/2;

distance_vec = zeros(10,1);
l =1;
p = input('type d when done',"s");
while(p ~= 'd')
display(l)



%We will start by calculating the power of the noise
white_n = zeros(N,1);
white_n(1) = 0.0001;
[ simin,nbsecs,fs] = initparams(white_n,fs,3);
sim('recplay');
% Retrieve recorded output
out=simout.signals.values(:,1);
%Calculate PSD
[s1,f1,t1,P] = spectrogram(white_n,windowS,Noverlap,windowS,fs);
PSD_rec_noise = sum(P,2);



%We will start by calculating the power of the signal
t = linspace(0,duration,N);


sinewave = sin(2*pi*1500*t)';
[ simin,nbsecs,fs] = initparams(sinewave,fs,3);
sim('recplay');
% Retrieve recorded output
out2=simout.signals.values(:,1);

%Trim the output so we only have the sine
%Find the max of the out signal
c2 = max(out2);
c_60 = 0.6*c2;
x = 0;
%Find the first x value that reaches this threshold

for i = 1:size(out2,1)
    if(out2(i) > c_60)
        x = i;
        break
    end    
end

yOnset = x; % Determine start of recorded signal [samples]

y = out2(yOnset-200:yOnset + N-201); % Extract the relevant output signal

%Calculate PSD
[s2,f2,t2,P2] = spectrogram(y,windowS,Noverlap,windowS,fs);
PSD_rec_sig = sum(P2,2);



PSD_rec_sig = PSD_rec_sig - PSD_rec_noise;
for k = 1:size(PSD_rec_sig)
    if PSD_rec_sig(k) < 0
        PSD_rec_sig(k) = 0;
    end
end




%Calculate the channel capacity;
C_channel = 0;
for k = 1:windowS/2

    C_channel = C_channel + log2(1 + PSD_rec_sig(k)/PSD_rec_noise(k));
end
C_channel = C_channel*fs/windowS;
display(C_channel)





distance_vec(l) = C_channel;
%Test the dependancy between distance & capacity
l = l+1; 
p = input('type d when done',"s");
end


display(distance_vec)
