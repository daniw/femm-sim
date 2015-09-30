i = 20;

load('bn.mat');load('bt.mat');load('hn.mat');load('ht.mat');

plot(1:length(hn(:,i)),hn(:,i), 1:length(ht(:,i)),ht(:,i));
