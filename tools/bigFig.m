function [ h1 ] = bigfig()
%BIGFIG - makes a figure window that covers the entire screen.

scrz=get(0,'ScreenSize');
h1=figure('Position',scrz);


end

