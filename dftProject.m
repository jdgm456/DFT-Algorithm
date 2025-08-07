%Algorithm for plotting Discrete Fourier Transform given a signal equation,
%duration, and number of padded zeroes
clc; %clear command window
clear; %clear variables
close all; %close all figures
syms t;




%Function will discretize the signal on the basis of fc and fs
%we wish to traverse from 0 to 1/fc - 1/fs samples with each 
%sample being evenly spaced by the sampling period 1/fs, 
%we always stop 1 sample before reaching 1/fc when discretizing 
%to count a cycle
function xn = discretize(fc,fs, sinusoid)
syms t;
   
    %discrete time steps evenly spaced by 1/fs, going from 0 to 1/fc
    n = 0:1/fs:1/fc-(1/fs);
    %change of variables
    sinusoid = subs(sinusoid,t,n);
    %casting xn as a double significantly improves computational time
    xn = double(sinusoid);
    
end


% Function will repeat the signal on the basis of however many times
%cycles are input, as well as padd zeros after repeating
    function xn_padded = cycleszeroes(cycles, z, xn)
    xn_padded = []; %pass in empty array that will be populated
    %we wish to repeat xn cycles's amount of times
    for i=1:cycles
        xn_padded = [xn_padded,xn];
    end
    %pad zeros after
    xn_padded = [xn_padded, zeros(1,z)];
    end

    %keeping track of zeros padded and sampling frequency
    %we can DFT padded xn, to preserve amplitude we normalize 
    %via length of signal, excluding padded zeros
    %frequency is normalized via length of signal plus
    %padded zeros to get the increased resolution effect
    function [xdft,f] = DiscreteFT(xn_padded,z,Fs)
        % w is a constant term, can be computed once and used in every
        % iteration
        w = exp(-1i*2*pi/(length(xn_padded)));
        %size of xdft is based on size of full time domain signal including
        %zeros
        xdft = zeros(1,length(xn_padded));
        
        %we wish to match frequency index length with xn padded
        for k=1:length(xn_padded)
            xdftsum = 0;
            
            %we also need to execute amount of xn padded iterations
            for n=1:length(xn_padded)
                %k-1 and n-1 so that we are still starting from 0, 
                %but matlab starts array index at 1
                xdftsum = xdftsum + xn_padded(n) * w^((k-1)*(n-1));
            end
            %Conversion factor to retrieve original amplitude, 
             %z is subtracted out to look at just non zero samples
            xdft(k) = abs(xdftsum*2/(length(xn_padded)- z ));
        end
        %Conversion factor to map frequency index onto hertz, divide by
        %length of non zero samples to get increased resolution!!
        f = (0:length(xn_padded)-1)*(Fs/(length(xn_padded) ));
       
    end

    %Test case a)
    xn1 = discretize(200,2000,47*cos(2*pi*200*t));
    xn1_padded = cycleszeroes(1,20,xn1);
    [xdft1,f1] = DiscreteFT(xn1_padded,20,2000);
    figure(1);
    hold on
    stem(f1,xdft1);
    %Fs/2 = 1000
    xlabel('frequency [hertz]');
    ylabel('Amplitude');
    title('DFT of a)');
    xlim([0 1000]);

    %Test case b)
    xn2 = discretize(200,2000,47*cos(2*pi*200*t));
    xn2_padded = cycleszeroes(1,200,xn1);
    [xdft2,f2] = DiscreteFT(xn2_padded,200,2000);
    figure(2);
    hold on
    stem(f2,xdft2);
    %Fs/2 = 1000
    xlabel('frequency [hertz]');
    ylabel('Amplitude');
    title('DFT of b)');
    xlim([0 1000]);

    %Test case c)
    xn3 = discretize(200,2000,47*cos(2*pi*200*t));
    xn3_padded = cycleszeroes(100,200,xn1);
    [xdft3,f3] = DiscreteFT(xn3_padded,200,2000);
    figure(3);
    hold on
    stem(f3,xdft3);
    %Fs/2 = 1000
    xlabel('frequency [hertz]');
    ylabel('Amplitude');
    title('DFT of c)');
    xlim([0 1000]);

    %Test case d)
    xn4 = discretize(400,2000,50*cos(2*pi*200*t)*sin(2*pi*200*t));
    xn4_padded = cycleszeroes(400,4000,xn4);
    [xdft4,f4] = DiscreteFT(xn4_padded,4000,2000);
    figure(4);
    hold on
    stem(f4,xdft4);
    %Fs/2 = 1000
    xlabel('frequency [hertz]');
    ylabel('Amplitude');
    title('DFT of d)');
    xlim([0 1000]);


    sinusoide = 47*cos(2*pi*300*t) + 50*sin(2*pi*400*t) + 74*cos(2*pi*1500*t) + 5*cos(2*pi*1600*t);
    %Test case e)
    xn5 = discretize(100,4000,sinusoide);
    xn5_padded = cycleszeroes(50,8000,xn5);
    [xdft5,f5] = DiscreteFT(xn5_padded,8000,4000);
    figure(5);
    hold on
    stem(f5,xdft5);
    %Fs/2 = 1000
    xlabel('frequency [hertz]');
    ylabel('Amplitude');
    title('DFT of e)');
    xlim([0 2000]);

    %Test case e)
    xn6 = discretize(100,1200,sinusoide);
    xn6_padded = cycleszeroes(50,8000,xn6);
    [xdft6,f6] = DiscreteFT(xn6_padded,8000,1200);
    figure(6);
    hold on
    stem(f6,xdft6);
    %Fs/2 = 1000
    xlabel('frequency [hertz]');
    ylabel('Amplitude');
    title('DFT of f)');
    xlim([0 600]);



    
                

