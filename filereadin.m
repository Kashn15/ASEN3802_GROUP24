clear;
clc;
close all;
dataFolder = 'ASEN3802_HeatConduction_FA25';
a = dir(fullfile(dataFolder, '*mA'));
A = (pi/4) * (0.0254^2);
for i = 1:length(a)
    %% Parse filename
    
    b = strsplit(a(i).name,'_');
material = b{1};

if strcmp(material,'Aluminum')
    k = 130;
elseif strcmp(material,'Brass')
    k = 115;
elseif strcmp(material,'Steel')
    k = 16.2;
end

    v = strsplit(b{2},'V');
    ampval = strsplit(b{3},'mA');
    volts(i) = str2double(v{1});
    amps(i)  = str2double(ampval{1});
    power = volts(i) * (amps(i))*1E-3; %convert to mA

    %% Import data file

    data = readmatrix(fullfile(dataFolder, a(i).name));
    % Delete first row
    data(1,:) = [];
    time = data(:,1);
    Dat = data(:,2:9);   % CH1–CH8

    %% steady state
    N = length(time);
    steady_index = round(0.9*N):N-2;   % last 10% of data assumed steady (we took the last 2 off which were NaN)

    Dat_steady = mean(Dat(steady_index,:),1);

    %% Thermocouple positions (meters)
x1= 1+3/8;
x2 = x1+0.5;
x3 = x2+0.5;
x4 = x3+0.5;
x5 = x4+0.5;
x6 = x5+0.5;
x7 = x6+0.5;
x8 = x7+0.5;
    x = 0.0254.*[x1 x2 x3 x4 x5 x6 x7 x8]; %METERS

    %% Linear fit for experimental

    p = polyfit(x, Dat_steady, 1);

    H_exp = p(1);
    T_0 = polyval(p,0);

    %% Analytical
    H_an = power/(k*A);

    %% Tables 
T_0_table(i) = T_0;
H_exp_table(i) = H_exp;
H_an_table(i) = H_an;
Table_Task1(i,:) = [T_0_table(i), H_exp_table(i), H_an_table(i)];
    %% Plot
    figure;
    plot(x, x.*H_an+T_0);
    hold on
    plot(x, polyval(p,x))
    errorbar(x,Dat_steady, 2*ones(size(Dat_steady)))
    xlabel('Position (m)')
    ylabel('Temperature (°C)')
    title([a(i).name ' Steady State'])
    grid on
hold off

%% Task 2
initial_index = 1:round(0.1*N);   % first 10% of data
T_init = mean(Dat(initial_index,:),1);

p_init = polyfit(x, T_init, 1);
M_exp = p_init(1);
M_0 = polyval(p_init,0);

figure;
scatter(x, T_init, 50)
hold on
plot(x, polyval(p_init,x))
xlabel('Position (m)')
ylabel('Temperature (°C)')
title([a(i).name ' Initial Temperature Distribution'])
legend('Initial Data','Linear Fit')
grid on

    %% Tables
    Table_Task2(i) = M_0;
    
end

