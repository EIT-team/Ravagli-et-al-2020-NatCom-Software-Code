
clc;
clear all;
close all;

files=dir('*_lite*.mat'); 

files={files.name}; files';

n=1;

load(files{n});files{n}

%% Concatenation of data

dZ=[];
Prt_0=[];
BV=[];
DC=[];
EP=[];
max_abs=[];

for i = 1:length(EIT)
    Prt_0=[Prt_0; [repmat(EIT{i}.Pair,29,1),[2:30]',ones(29,1)]];
    dZ = [dZ; EIT{i}.dZ_avg'];
    BV=[BV; EIT{i}.BV0 ];
    DC=[DC; EIT{i}.DC_avg];
    EP=[EP; EIT{i}.EP_avg'];
    max_abs=[max_abs; EIT{i}.ch_max_abs'];
end

sprintf('Rows: %d - Original',size(Prt_0,1))

%% Remove recordings from reference electrode (Ch19 on ActiChamp)

ind = ~any(Prt_0'==19);
Prt_0 = Prt_0(ind,:);
dZ = dZ(ind,:);
BV = BV(ind,:);
DC = DC(ind,:);
EP = EP(ind,:);
max_abs=max_abs(ind,:);

sprintf('Rows: %d - Ref19',size(Prt_0,1))

%% Remap only Prt matrix according to numbering on forward model

% Forward model numbering in ActiChamp HW positions - Example: Acti=4, forward=18, so forward_on_acti(4)=18;

forward_on_acti= [ 29	9	3	18	24	14	8	2	19	25	13	7	4	1	17	20	23	26	29	12	6	15	21	27	11	5	16	22	28	10];
Prt_0=forward_on_acti(Prt_0);

%% Plot original

figure(1); hold on;
subplot(231);
plot(T,dZ); grid on;
title('Original');
subplot(234);
plot(T,EP);grid on;

%% Plot BV0s for check of saturated electrodes

% figure(2); hold on;
% bar(max_abs);
% ind=max_abs>400e3;
% Prt_0(ind,:)

%% Removal of pairs with DC saturation (keep DC level <400mV)

ind = abs(DC)<400e3;

Prt_0 = Prt_0(ind,:);
dZ = dZ(ind,:);
BV = BV(ind,:);
DC = DC(ind,:);
EP = EP(ind,:);
max_abs=max_abs(ind,:);

sprintf('Rows: %d - max_abs_channel>400mV',size(Prt_0,1))

%% Keep only data from ring1 (EIT ring used on cuff with two rings)

ind = any(Prt_0(:,1)'==[1:14]') & any(Prt_0(:,2)'==[1:14]') & any(Prt_0(:,3)'==[1:14]') ; 
Prt_0=Prt_0(ind,:);
dZ = dZ(ind,:);
BV = BV(ind,:);
DC = DC(ind,:);
EP = EP(ind,:);
max_abs=max_abs(ind,:);

sprintf('Rows: %d - Ring1',size(Prt_0,1))

%% How many injection pairs are left?

ind=(Prt_0(:,3)==Prt_0(:,1))  |  (Prt_0(:,3)==Prt_0(:,2));
sprintf('Left %d injection pairs ',nnz(ind))

%% Plot clean #1

figure(1); hold on;
subplot(232);
plot(T,dZ); grid on;
title('Clean #1');
subplot(235);
plot(T,EP);grid on;

%% Plot pre-stim noise std histogram to choose threshold for later

% figure(3); hold on;
% histogram(std(dZ(:,T<-3),0,2)',50);

%% Clean by pre-stim noise std

ind = std(dZ(:,T<-5),0,2)<5;
Prt_0 = Prt_0(ind,:);
dZ = dZ(ind,:);
BV = BV(ind,:);
DC = DC(ind,:);
EP = EP(ind,:);
max_abs=max_abs(ind,:);

sprintf('Rows: %d - Prestim noise std',size(Prt_0,1))

%% How many injection pairs are left?

ind=(Prt_0(:,3)==Prt_0(:,1))  |  (Prt_0(:,3)==Prt_0(:,2));
sprintf('Left %d injection pairs ',nnz(ind))

%% How many injection pairs are left?

% ind=(Prt_0(:,3)==Prt_0(:,1))  |  (Prt_0(:,3)==Prt_0(:,2));
% nnz(ind)

%% Plot clean #2

figure(1); hold on;
subplot(233);
plot(T,dZ); grid on;
% xlim([-2 2]); ylim([-20 20]);
title('Clean #1');
subplot(236);
plot(T,EP);grid on;
% xlim([-2 2]);

%% dZ flip for recon format

dZ=dZ';

%% Save results

save(['RECON_' files{n}(1:end-4) '.mat' ],'BV','Prt_0','T','dZ','Fc','Fs');






























