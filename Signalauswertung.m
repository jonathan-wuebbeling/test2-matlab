clear all; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%       script for showing different effect when measuring a vibartion
%       with different Sensors an using of differentiation and integartion
%
%       To get the raw data you have to perform the Lab No.3 of the
%       Vibrationmeasurement and Analysis lecture
%
%       This script contains several plots with methods of signal
%       processing. To the the different effect  it is recommended to watch
%       the different plot one after another using breakpoints. For every
%       plot the are some tasks. These Tasks are described in the
%       description in front of every plot section.
%
%
%       The evaluated *.csv file comes from an Picoscope
%
%       Channel A : acceleration sensor
%       Channel B : velocity sensor
%       Channel C : displacement sensor
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Settings
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%filename = '20180118-0005_25Hz_T50ms.csv';      % the CSV file has to be in the same folder as the script itself
%filename = '20190619-0001\20190619-0001_1.csv';      % 
filename = '20190619-0003.csv';
Measurement.SignalFrequency = 20.0;             % frequency  of signal (as configured in the frequency  generator)
SensitivityDisplacementSensor = 10;            % Displacement  in   V/mm
SensitivityVelocitySensor = 0.0233;             % Velocity      in   V/mm/s
SensitivityAccelerationSensor = 10.16*10^-6;    % Acceleration  in   V/mm/s²










%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   read *.csv file from Picoscope Measurement
%   the first have was crated automaticly with a matlab tool to read the
%   *.csv from the Picoscope
%   the second half calculates the sensorreadings from V to the
%   correponding physical units
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delimiter = ',';
startRow = 3;
formatSpec = '%f%f%f%f%[^\n\r]';
fileID = fopen(filename,'r');
textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false, 'EndOfLine', '\r\n');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'ReturnOnError', false);
fclose(fileID);
Array = importdata(filename);
%the measurement is rearenaged into simple vectors and conversion from
%Volts to physical units with Sensitivity is done
Measurement.Time = Array.data(:,1)/1000.;                                       % time [ms] --> [s]
Measurement.Acceleration = (Array.data(:,2)./1000.0)./SensitivityAccelerationSensor;    % acceleration [mV] --> [V]
Measurement.Velocity = Array.data(:,3)./SensitivityVelocitySensor;                  % velocity [V]
Measurement.Displacement = Array.data(:,4)./SensitivityDisplacementSensor;              % displacement [V]

clear('Array','dataArray','formatSpec','startRow','filename','fileID','delimiter','SensitivityDisplacementSensor','SensitivityVelocitySensor','SensitivityAccelerationSensor','startRow','ans') % delete obsolete variables from Workspace

%getting some general informations from Measurement which are needed later on
Measurement.TimeStep = (Measurement.Time(end)-Measurement.Time(1)) / length(Measurement.Time);        %calculation of timestep with average
Measurement.SamplingFrequency = 1./Measurement.TimeStep; 

% Analytische Funktionen ggf. zu Demonstrationszwecken 
%x = sin(2*pi*20*t);
%v = 2*pi*20*cos(2*pi*20*t+pi/3.);
%a = -4*pi^2*20*20*sin(2*pi*20*t+pi/3.);










%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Plot No. 1
%
%   this plot shows the raw date of the different Signal
%
%   ! have a look at the magnitude and the phase differences !
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
%plot all three Measurements in in graph
subplot(4,1,1)
plot(Measurement.Time,Measurement.Displacement,'g-',Measurement.Time,Measurement.Velocity,'r-',Measurement.Time,Measurement.Acceleration,'b-') 
title('Signals')
legend('Displacement','Velocity', 'Acceleration')
xlabel('time [s]') % x-axis label#
ylabel('$mm \ $,$\ \displaystyle\frac{mm}{s} \ $,$\ \displaystyle\frac{mm}{s^2}$,','interpreter','latex') 

subplot(4,1,2)
plot(Measurement.Time,Measurement.Displacement,'g-') 

legend('Displacement')
xlabel('time [s]') % x-axis label#
ylabel('$mm$','interpreter','latex') 

subplot(4,1,3)
plot(Measurement.Time,Measurement.Velocity,'r-') 

legend('Velocity')
xlabel('time [s]') % x-axis label#
ylabel('$\displaystyle\frac{mm}{s}$','interpreter','latex') 

subplot(4,1,4)
plot(Measurement.Time,Measurement.Acceleration,'b-') 

legend('Acceleration')
xlabel('time [s]') % x-axis label#
ylabel('$\displaystyle\frac{mm}{s^2}$','interpreter','latex') 








%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Plot No. 2
%
%   in the first section the Acceleration is calculated from the Signal of
%   the velocity sensor by differentiating
%   in the second section the calculated Velocity is plotted and
%   compared to the measured acceleration
%   
%   ! zoom in !
%   ? How do you explain the differences ?
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%In the next section the acceleration is calculated from velocity by differentiation 
Calculation.AccelerationFromVelocity = diff(Measurement.Velocity)./Measurement.TimeStep;              %calculation of acceleration as derivative of velocity
Calculation.Time = Measurement.Time(1:length(Measurement.Time)-1);                                    %getting a new time vector for the derivative


%plot of the calculated 
figure                             
plot(Calculation.Time,Calculation.AccelerationFromVelocity,'Magenta-',Measurement.Time,Measurement.Acceleration,'b-') 
title('Acceleration')
legend('Derivation of Velocity','Acceleration Sensor Measurement')
xlabel('time [s]') % x-axis label
ylabel('$\displaystyle\frac{mm}{s^2}$','interpreter','latex') % y-axis label





%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Plot No. 3
%
%   In Plot No. 2 there is a lot of noise pruduced by differentiation. This
%   can be avoided by appling a Low pass filter prior to the
%   differentiation
%   
%   ! use a Calculation.LowPassMultiplier of 1 !
%   ? how does filtering change the velocity Measurement , zoom in ?
%   ! change Calculation.LowPassMultiplier to a smaller number e.g. 10 !
%   ? whats happening an why ?
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Lowpass Filter using the Matlab Butter filter with a Factor between
%Sgnal Frequency an Cut off frequency of 1
Calculation.LowPassCutOffFrequencyFactor1 = 1*Measurement.SignalFrequency;         % calculation of Cut off frequency from Signalfrequency and a prior defined factor.
[b,c] = butter(1.,Calculation.LowPassCutOffFrequencyFactor1/(Measurement.SamplingFrequency/2.),'low');         % generate Low pass Butter filter of first order from the matlab Toolbox
Calculation.VelocityWithLowPassFactor1 = filter(b,c,Measurement.Velocity);                                     % apply filter to measured signal


%Factor of 10
Calculation.LowPassCutOffFrequencyFactor10 = 10*Measurement.SignalFrequency;         % calculation of Cut off frequency from Signalfrequency and a prior defined factor.
[b,c] = butter(1.,Calculation.LowPassCutOffFrequencyFactor10/(Measurement.SamplingFrequency/2.),'low');         % generate Low pass Butter filter of first order from the matlab Toolbox
Calculation.VelocityWithLowPassFactor10 = filter(b,c,Measurement.Velocity);                                     % apply filter to measured signal



figure                              
plot(Measurement.Time,Measurement.Velocity,'r-',Measurement.Time,Calculation.VelocityWithLowPassFactor1,'black-',Measurement.Time,Calculation.VelocityWithLowPassFactor10,'g') 
title('Velocity')
legend('Velocity Measurement',['Filtered Velocity Measurement f_c = ' num2str(Calculation.LowPassCutOffFrequencyFactor1)],['Filtered Velocity Measurement f_c = ' num2str(Calculation.LowPassCutOffFrequencyFactor10)])
xlabel('time [s]') % x-axis label
ylabel('$\displaystyle\frac{mm}{s}$','interpreter','latex') % y-axis label


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Plot No. 4
%
%   in this plot the filtered velocity is derived and compared to the
%   measured velocity
%
%   ! compare to plot No. 2 !
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Calculate acceleration from filtered velocity by differentiation as in plot No. 2
Calculation.AccelerationFromVelocity = diff(Calculation.VelocityWithLowPassFactor10)./Measurement.TimeStep;
Calculation.Time = Measurement.Time(1:length(Measurement.Time)-1);

figure
plot(Calculation.Time,Calculation.AccelerationFromVelocity,'Magenta-',Measurement.Time,Measurement.Acceleration,'b-') 
title('Acceleration')
legend('Derivation of filtered Velocity','Acceleration Sensor Measurement')
xlabel('time [s]') % x-axis label
ylabel('$\displaystyle\frac{mm}{s^2}$','interpreter','latex') % y-axis label
ylim([2*min(Measurement.Acceleration) 2*max(Measurement.Acceleration)])






%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Plot No. 5
%
%   This plot shows the effect of integration. Therefor the measured
%   acceleration is integrated two times and compared to velocity
%   and Displacement measurements
%
%   ? whate effect do you see ?
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Calculate velocity from acceleration by integration
Calculation.VelocityFromAcceleration = Measurement.TimeStep*cumtrapz(Measurement.Acceleration);

figure
subplot(2,1,1)
plot(Measurement.Time,Calculation.VelocityFromAcceleration,'Magenta-',Measurement.Time,Measurement.Velocity,'b-')
title('Velocity')
legend('Integral of Acceleration','Velocity Measurement')
xlabel('time [s]') % x-axis label
ylabel('$\displaystyle\frac{mm}{s}$','interpreter','latex') % y-axis label

% Calculate displacement from allready integrated acceleration by integration
Calculation.DisplacementFromAcceleration = Measurement.TimeStep*cumtrapz(Calculation.VelocityFromAcceleration);


subplot(2,1,2)
plot(Measurement.Time,Calculation.DisplacementFromAcceleration,'Magenta-',Measurement.Time,Measurement.Displacement,'b-')
title('Displacement')
legend('Second Integral of Acceleration','Displacement Measurement')
xlabel('time [s]') % x-axis label
ylabel('$mm$','interpreter','latex') % y-axis label





%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Plot No. 6
%
%   one way to coorect the drift is to substract a polynomial fit from the
%   integrals
%
%   ! compare plot 6 to plot 7 !
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%correction of velocity = integral of acceleration by polyfir of first
%order
Calculation.VelocityFromAccelerationPolyFit1 = polyval(polyfit(Measurement.Time,Calculation.VelocityFromAcceleration,1),Measurement.Time); %Linearen Anteil aus Integration abziehen 
Calculation.VelocityFromAccelerationPolyFitCorrected = Calculation.VelocityFromAcceleration - Calculation.VelocityFromAccelerationPolyFit1;

%correction of displacement = second integral of acceleration by polyfir of
%second order
Calculation.DisplacementFromAccelerationPolyFit2 = polyval(polyfit(Measurement.Time,Calculation.DisplacementFromAcceleration,2),Measurement.Time); %Linearen Anteil aus Integration abziehen 
Calculation.DisplacementFromAccelerationPolyFitCorrected = Calculation.DisplacementFromAcceleration - Calculation.DisplacementFromAccelerationPolyFit2;



figure
subplot(2,1,1)                 
plot(Measurement.Time,Calculation.VelocityFromAcceleration,'Magenta-',Measurement.Time,Calculation.VelocityFromAccelerationPolyFitCorrected,'b-')
title('Velocity')
legend('Integral of Acceleration','Integral od Acceleration corrected with Polynomial Fit of first Order')
xlabel('time [s]') % x-axis label
ylabel('$\displaystyle\frac{mm}{s}$','interpreter','latex') % y-axis label

subplot(2,1,2)
plot(Measurement.Time,Calculation.DisplacementFromAcceleration,'Magenta-',Measurement.Time,Calculation.DisplacementFromAccelerationPolyFitCorrected,'b-')
title('Displacement')
legend('Second Integral of Acceleration','Second Integral od Acceleration corrected with Polynomial Fit of second Order');
xlabel('time [s]') % x-axis label
ylabel('$mm$','interpreter','latex') % y-axis label

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Plot No. 7
%
%   still to be done : use a low pass Filter to correct the drift.
%   substract mean Value before integration
%
%   ! task !
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


help butter
