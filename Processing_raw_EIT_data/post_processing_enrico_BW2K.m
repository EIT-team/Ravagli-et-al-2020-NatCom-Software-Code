clc;
close all;
clear all;

files=dir('*EIT*.vhdr');

files={files.name};

for ffil = 1:length(files)
    
    EIT_fname = files{ffil};
    
    log_f=dir([EIT_fname(1:end-5) '_log*.mat']);
    log_f=log_f.name;
    load(log_f);
    
    % Electrode numbering on Actichamp
    map_p=[1:29];
    map_=[2:30];
    
    % Parameters for post-processing
    EP_cutoff = 2000;  % cutoff frequency for EPs (low-pass freq.)
    dZ_BW = 2000;      % Bandwidth for demodulation
    N_butter_EP = 5; % Butterworth order for filtering EPs
    N_butter_dZ = 5;  % Butterworth order for filtering dZs
    N_butter_notch = 3; % Butterworth order for notch filtering
    Noise_thres = 1000;  % Noise threshold (in uV)
    Filter_50Hz = true; % 50 Hz notch filter on EPs
    Plot_chan = [1:length(map_)];
    T_window = 0.05; % time to take around trigger
    
    %% Find all protocol switches
    
    HDR= sopen([EIT_fname(1:end-5) '.eeg']);
    Fs=HDR.SampleRate;
    
    SWITCH=HDR.EVENT.POS;
    
    for i=1:length(SWITCH)
        if ~(strcmp([HDR.EVENT.Desc{i}], ['S  3']) || strcmp([HDR.EVENT.Desc{i}], ['S  1']))
            SWITCH(i)=0;
        end
    end
    
    SWITCH(SWITCH==0)=[];
    SWITCH=SWITCH(2:end);
    
    SWITCH = [SWITCH;2*SWITCH(end)-SWITCH(end-1)];
    
    Prt = ExpSetup.Protocol;
    Prt_size = size(Prt,1);
    
    if (Prt_size+1==size(SWITCH,1))
        disp('All data seemed to be there');
    elseif (Prt_size+1>size(SWITCH,1))
        disp('Protocol size is bigger than the number of switches in the file');
        Prt_size=size(SWITCH,1)-1;
    else
        disp('Protocol size is smaller than the number of switches in the file!'); %%, assuming repeated protocol');
        %%write stuff to handle this
    end
    
    EEG = pop_loadbv('',EIT_fname,[SWITCH(1)+Fs SWITCH(1)+10*Fs]);
    
    EEG.data=EEG.data(map_p,:);    
    
    %% Not clear why it looks for injection on map_==Prt(1,2) but it is only for carrier detection so don't really care - Enrico
    
    V_inj = detrend(double(EEG.data(map_==Prt(1,1),:)'),'constant');
    NFFT = 2^nextpow2(length(V_inj));           % Next power of 2 from length of y
    Y = fft(V_inj,NFFT)/length(V_inj);
    f = Fs/2*linspace(0,1,NFFT/2+1);
    w_inj=2*abs(Y(1:NFFT/2+1));
    
    [~,maxw] = max((w_inj));
    Fc = f(maxw);
    disp(sprintf('****** Detected carrier frequency: Fc = %i Hz ******',Fc));    
    
    %% Processing on each pair
    
    EIT = [];
    hhplots = figure('Position',[10,50,1900,950],'PaperPositionMode','auto');
    
    for iPair = 1:Prt_size
        
        disp(sprintf('Processing protocol line %02i / %02i',iPair,Prt_size));
        
        EEG = pop_loadbv( '',EIT_fname,[SWITCH(iPair) SWITCH(iPair+1)]);
        
        EEG.data=EEG.data(map_p,:);
        
        T_trig=cell2mat({EEG.event.latency})';
               
        for i=1:length(T_trig)
            if ~strcmp([EEG.event(1,i).type], ['S  2']);
                T_trig(i)=0;
            end
        end
        
        T_trig(T_trig==0)=[];
        T_trig=T_trig(4:end-3);
        
        T_stim = mean(T_trig(2:end)-T_trig(1:end-1))/Fs;
        disp(['Stimulation every ' num2str(round(T_stim*1000)) ' ms']);
        
        Data = double(EEG.data');
        
        % Low-pass filter to retrieve EP's
        [b,a] = butter(N_butter_EP,EP_cutoff/(Fs/2),'low');
        X_ep = filtfilt(b,a,Data);
        
        if Filter_50Hz
            [b,a] = iirnotch(50/(Fs/2),(50/(Fs/2))/35);
            X_ep = filtfilt(b,a,X_ep);
        end
        
        % Band-pass filter and demodulate
        [b,a] = butter(N_butter_dZ,(Fc+dZ_BW*[-1,1])/(Fs/2));
        X_dz = filtfilt(b,a,Data);
        X_dz = hilbert(X_dz);
        A_dz   = abs(X_dz);              % amplitude
        
        %%%%%%%%%%%% Code for removal of motion artifact - Enrico %%%%%%%%%%%%%%%%%
        
      A_dz_mean=mean(A_dz);
%       [b,a] = butter(3,[100]/(Fs/2),'high');
%       A_dz_highpass = filtfilt(b,a,A_dz);
%       A_dz=A_dz_mean+A_dz_highpass;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Number of channels, triggers, bins
        N_chan = size(Data,2);
        N_trig = length(T_trig);
        N_bin = round(T_window*Fs);
        w = (1:N_bin)-round(N_bin*0.5);    % window around trigger
        T = 1e3*w/Fs;
        t0 = find(T>-5 & T<-1);
        
        % Get results of each channel in matrices
        EP   = cell(1,N_chan);
        DC   = cell(1,N_chan);
        dZ   = cell(1,N_chan);
        BV = cell(1,N_chan);
        for iChan = 1:N_chan
            EP{iChan} = zeros(N_bin,N_trig);
            dZ{iChan} = zeros(N_bin,N_trig);
            BV{iChan} = zeros(1,N_trig);
            
            % iterate over triggers
            for jTrig = 1:N_trig
                ival = T_trig(jTrig)+w;
                EP{iChan}(:,jTrig) = detrend(X_ep(ival,iChan),'constant');
                DC{iChan}(:,jTrig) = mean(X_ep(ival,iChan));
                dZ{iChan}(:,jTrig) = detrend(A_dz(ival,iChan),'constant');
                BV{iChan}(jTrig) = mean(A_dz(ival,iChan));
            end
        end
        
        % Compute average / std. error for EP, dZ, dPhi
        avg_EP = zeros(N_bin,N_chan);
        avg_dZ_rel = zeros(N_bin,N_chan);
        avg_dZ_abs = zeros(N_bin,N_chan);
        BV0=zeros(length(map_p),1);
        DC_avg=zeros(length(map_p),1);
              
        for iChan = 1:N_chan
            ep = EP{iChan};
            dz = dZ{iChan};
            bv = repmat(BV{iChan},N_bin,1);            
            avg_EP(:,iChan) = mean(ep,2)-mean(reshape(ep(t0,:),[],1));
            avg_dZ_abs(:,iChan) = mean(dz,2)-mean(reshape(dz(t0,:),[],1));
            avg_dZ_rel(:,iChan) = 100*mean(dz./bv,2);
            avg_dZ_rel(:,iChan) = avg_dZ_rel(:,iChan)-mean(avg_dZ_rel(t0,iChan));
            BV0(iChan)=mean(mean(bv,2));   
            DC_avg(iChan)=mean(DC{iChan});
        end
              
        EIT{iPair}.Pair = Prt(iPair,:);
        EIT{iPair}.EP_avg=avg_EP;
        EIT{iPair}.dZ_avg=avg_dZ_abs;
        EIT{iPair}.BV0=BV0;
        EIT{iPair}.DC_avg=DC_avg;
        EIT{iPair}.ch_max_abs=max(abs(Data));
        
        subplot(round(sqrt(Prt_size))+1,round(Prt_size/round(sqrt(Prt_size))),iPair);
        
        % Original plot (selective)
        Plot_chan = find(std(avg_dZ_abs(T<-1,:))<100000);
        if(nnz(Plot_chan))>0
            plot(T,avg_dZ_abs(:,Plot_chan));
            title(['Pair=' num2str(iPair)]);
            ylim([-20,20])
            xlim([-5,25])
            grid on
            ylabel('dZ (uV)')
        end
        
        strtt = EIT_fname(1:end-4);
        strtt(strtt=='_')= ' ';
        drawnow;
        
    end    
%     save([EIT_fname(1:end-5) '_lite_BW2K_hpfilt100.mat'],'EIT','T','Fs','Fc','-v7.3');
    save([EIT_fname(1:end-5) '_lite_BW2K.mat'],'EIT','T','Fs','Fc','-v7.3');
    
end







