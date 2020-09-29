
clc;
clear all;
close all;

%% EXP
load('test_LL_sref.mat');

%% EXP - Remove damaged electrodes

% ind = any(Prt_exp'==12 ) | any(Prt_exp'==12 )| any(Prt_exp'==12 );
% ind = ~ind;
% 
% Prt_exp=Prt_exp(ind,:);
% BV_exp=BV_exp(ind,:);

%% MOD

load('BV0_forward_SingleRef_full.mat');
% load('BV0_forward_DoubleRef_full.mat');

BV_mod=BV0;
Prt_mod=Prt_Forward;

%% MOD - Only ring 1

 ind = any(Prt_mod(:,1)'==[1:14]') & any(Prt_mod(:,2)'==[1:14]') & any(Prt_mod(:,3)'==[1:14]') ; 
 Prt_mod = Prt_mod(ind,:);
 BV_mod = BV_mod(ind,:);

%%  Intersect

[C,ind_exp,ind_mod] = intersect(Prt_exp(:,1:3),Prt_mod(:,1:3),'rows');

Prt_exp=Prt_exp(ind_exp,:);
BV_exp=BV_exp(ind_exp,:);

Prt_mod=Prt_mod(ind_mod,:);
BV_mod=BV_mod(ind_mod,:);

%%

scatter(abs(BV_mod),abs(BV_exp*1e-6));
R=corrcoef(abs(BV_mod),abs(BV_exp))
















