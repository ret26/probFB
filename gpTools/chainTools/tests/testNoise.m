close all
clear all
clc
ntypes = {'babble','factory1','volvo','buccaneer1','leopard'};
ndir = './datasets/audio/noise/preMIRS/';
for k = 1:length(ntypes)
    rawfile = [ndir ntypes{k} '.raw'];
    fid = fopen(rawfile,'r');
    data = fread(fid,Inf,'int16');
    fs = 8000;
    figure, plot(data);
    keyboard
    fclose(fid);
end
