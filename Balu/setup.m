clear all
close all
clc

BaluPath = pwd;

io = 'InputOutput';

addpath(fullfile(BaluPath));
addpath(fullfile(BaluPath,io));

I = imresize(imread('LogoBalu.png'),0.5);
imshow(I)
title('Installing Balu 4.0 ...');



d = dir;

n = length(d);

ft = Bio_statusbar('Installing Balu 4.0');

for i=1:n
    ft = Bio_statusbar(i/n,ft);
    st = d(i).name;
    if and(exist(st,'dir'),length(st)>2)
        pause(0.2)
        fprintf('Adding directory %s...\n',st);
        if strcmp(io,st)~=1
            addpath(fullfile(BaluPath,st));
        end
    end
end
delete(ft)
savepath

title('Balu 4.0 installed');

fprintf('Balu 4.0 installed succefully!\n');

pause(2)
close all

