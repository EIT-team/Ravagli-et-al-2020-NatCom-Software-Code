
clc;
clear all;

load('test_RL_pp');

acti_on_forward=[ 14 8	3	13	26	21	12	7	2	30	25	20	11	6	22	27	15	4	9	16	23	28	17	5	10	18	24	29 ] -1;

N_plots=length(EIT)/7;

figure;

for i=1:14
       
    subplot(N_plots,7,i);    
    
    BV=EIT{i}.BV0(acti_on_forward);
    
    bar(BV(1:14));   

    
end



