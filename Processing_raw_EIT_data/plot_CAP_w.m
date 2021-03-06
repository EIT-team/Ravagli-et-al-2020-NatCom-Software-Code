
% Plot compound action potentials (CAPs)

EEG=pop_loadbv;

Fs = EEG.srate;

map_p=[1:29]; % Montage
map_=[1:18,20:30]; %REF = 19

t=cell2mat({EEG.event.latency})';

jj=zeros(length(t),1);
for i=1:length(t)
    if strcmp([EEG.event(1,i).type], ['S  2']);
        jj(i)=1;
    end
end

T_trig=t(jj==1);
T_trig = T_trig(5:end-5)

tau=100;
size_bin=floor(tau*Fs/1000);

Data= double(EEG.data)';

[b,a] = butter(5,40000/(Fs/2),'low');
Data = filtfilt(b,a,Data);

T=[1:size_bin]*1000/Fs-tau/2;

EP=zeros(length(T_trig)-2,size_bin,size(Data,2));
for i=2:length(T_trig)-1
    EP(i-1,:,:)=Data([T_trig(i)-floor(size_bin/2):T_trig(i)+floor(size_bin/2)-1],:);
end

EP_avg=detrend(squeeze(mean(EP,1)));

figure
plot(T,EP_avg);

