clear all; close all; clc;

buffer=[1 2 3 4 5 6 7 8 9]

while ~isempty(buffer)
    k=buffer(1);
    buffer(1)=[];
    if(isprime(k))
        fprintf('****');
%         buffer=[buffer 4];
        buffer(end+1)=4;
    end
    fprintf('k=%d\n',k);
end
