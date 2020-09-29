
clc;
clear all;
% close all;

load('test_LL_sref.mat');
BV1=BV_exp;

load('test_LL_dref.mat');
BV2=BV_exp;

%% Check only specific parts of the protocol

% BV1=BV1(1:392); BV2=BV2(1:392);     % Ring1

%%

scatter(abs(BV1),abs(BV2));
R=corrcoef(abs(BV1),abs(BV2))























