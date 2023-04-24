close all;
clear;
clc;

%% enables

enable_LOS = 1;
enable_ground_reflections = 1;
enable_diffraction = 1;

enable_one_reflection = 1;
enable_left_wall_reflection = 1;
enable_right_wall_reflection = 1;
enable_bottom_wall_reflection = 1;
eanble_rue_du_grand_hospice = 1;
enable_rue_du_rouleau = 1;
enable_rue_du_peuplier = 1;

enable_double_reflections = 1;
enable_left_and_right_wall = 1;
enable_left_and_bottom_wall = 1;
enable_left_and_streets = 1;
enable_right_and_left_wall = 1;
enable_right_and_bottom_wall = 1;
enable_bottom_and_left_wall = 1;
enable_bottom_and_right_wall = 1;
enable_bottom_and_streets = 1;
enable_streets_and_streets = 1;


%% variables

Fc = 27e9; %Hz
EIRP = 1; %W Gtx*Ptx
h_bs = 2; %m
h_ue = 2; %m
t_SNR = 5; %dB
r_noise_figure = 15; %dB
BW = 200e6; %Hz
d_basestation = 10; %m
permitivity = 4; %between 3-5
c = 3e8; %m/s
Ra = 73; %Ohm
time_resolution = 1/BW; %s

ray_traced_position = [15.5 ; 128.5];
impulse_position = [15.5; 128.5]; %point where impulse response will be plotted


%% draw image

[x, y] = drawLines();

figure
hold on
lines = line(x, y, 'Color', 'black');
axis('equal')
xlabel('width [m]')
ylabel('height [m]')
title('Positions where received power is calculated')

%% transmitter position

tx_x = 20;
tx_y = 300;

%% raytracing for heatmap

% step 1: move rx_x over possible squares (create array for this)
% step 2: check for LOS (includes ground reflection)
% step 3: raytracing

%% step 1

all_rx_xy = getReceiverPositions(d_basestation, tx_x, tx_y);


plot(all_rx_xy(1, :), all_rx_xy(2, :), '.');

%% check if wanted_ray_traced_position is valid


if(isempty(find(and(all_rx_xy(1,:) == ray_traced_position(1), all_rx_xy(2,:) == ray_traced_position(2)))))
    disp('Not a valid raytracing position')
    return
end
if(isempty(find(and(all_rx_xy(1,:) == impulse_position(1), all_rx_xy(2,:) == impulse_position(2)))))
    disp('Not a valid impulse response position')
    return
end

%% step 2

no_LOS = [];
for i = 1:length(all_rx_xy(1,:))
    %turn path from rx to tx in line

    intersec = intersectionCalculator(x, y, all_rx_xy(:,i), tx_x, tx_y, 0);

    if (intersec == 1)
        no_LOS = [no_LOS [all_rx_xy(:,i); i]];
    end
end

figure
hold on
lines = line(x, y, 'Color', 'black');
axis('equal')
xlabel('width [m]')
ylabel('height [m]')
title('Positions where there is no line of sight')
plot(no_LOS(1, :), no_LOS(2, :), '.');


%% step 3

% step 3.1 calculate LOS
% step 3.2 calculate ground reflection
% step 3.3 calculate single reflection
% step 3.4 calculate double reflection
% step 3.5 calculate diffraction

%% step 3.1
% IF HEIGHT OF BS AND UA ARE DIFFERENT -> THIS RAYTRACING IS NOT CORRECT!!!!!!!!!!!!!!

received_power = zeros(size(all_rx_xy(1,:)));
effective_h = -(c/Fc)/pi*cos(pi/2*cos(pi/2))/sin(pi/2)^2; %same for every spot except ground reflection
%diffraction_points = [40, 280; 40, 190; 40, 80; 40, 10; 0, 10];
beta = 2*pi*Fc/c;
[received_power, rice_factor, delay_spread, impulse_response] = receivedPowerCalc(Ra, c, permitivity, x, y, all_rx_xy, tx_x, tx_y, effective_h, beta, h_ue, h_bs, EIRP, Fc, no_LOS, enable_LOS, enable_ground_reflections, ...
enable_diffraction, enable_one_reflection, enable_left_wall_reflection, enable_right_wall_reflection, enable_bottom_wall_reflection, eanble_rue_du_grand_hospice, enable_rue_du_rouleau, ...
enable_rue_du_peuplier, enable_double_reflections, enable_left_and_right_wall, enable_left_and_bottom_wall, enable_left_and_streets, enable_right_and_left_wall, enable_right_and_bottom_wall, ...
enable_bottom_and_left_wall, enable_bottom_and_right_wall, enable_bottom_and_streets, enable_streets_and_streets);

%calculate dBm of received power
received_power_dBm = 10*log10(received_power) + 30; 

%plot 3D structure
figure
plot3(all_rx_xy(1, :), all_rx_xy(2, :), received_power_dBm, '.')
daspect([1 1 0.5])
xlabel('width [m]')
ylabel('height [m]')
zlabel('received power (dBm)')
title('Received power as a 3D plot')

%plot heatmap
figure
scatter(all_rx_xy(1, :), all_rx_xy(2, :), [], received_power_dBm, 'filled')
hold on
lines = line(x, y, 'Color', 'black');
xlabel('width [m]')
ylabel('height [m]')
axis equal
c = colorbar;
c.Label.String = "Received power [dBm]";
title('Received power as a heatmap')

%plot BS to bottom wall
line_index = find(and(all_rx_xy(1,:) == 20.5, all_rx_xy(2,:) <290)); %find elements on line between BS and bottom wall
figure
plot(log10(all_rx_xy(2,line_index) - 300), received_power_dBm(line_index))
%set(gca, 'XDir','reverse')
xlabel('Distance from basestation [log(m)]')
ylabel('Received power [dBm]')
title('Received power between basestation and bottom wall')

%% calculate SNR

SNR  = zeros(size(received_power));
thermal_noise_power = 10*log10(1.379*10^-23*293*BW); %dB
SNR = received_power_dBm - 30 - r_noise_figure - thermal_noise_power;

figure
scatter(all_rx_xy(1, :), all_rx_xy(2, :), [], SNR, 'filled')
hold on
lines = line(x, y, 'Color', 'black');
xlabel('width [m]')
ylabel('height [m]')
axis equal
c = colorbar;
c.Label.String = "Signal to Noise Ration [dB]";
title('SNR as a heatmap')

%plot BS to bottom wall
figure
plot(log10(all_rx_xy(2,line_index) - 300), SNR(line_index))
%set(gca, 'XDir','reverse')
xlabel('Distance from basestation [log(m)]')
ylabel('Signal to Noise Ratio [dB]')
title('SNR between basestation and bottom wall')

%% plot rice factor

rice_factor_dBm = 10*log10(rice_factor) + 30;

figure
scatter(all_rx_xy(1, :), all_rx_xy(2, :), [], rice_factor_dBm, 'filled')
hold on
lines = line(x, y, 'Color', 'black');
xlabel('width [m]')
ylabel('height [m]')
axis equal
c = colorbar;
c.Label.String = "Rice Factor [dBm]";
title('Rice Factor as a heatmap')

%plot BS to bottom wall
figure
plot(log10(all_rx_xy(2,line_index) - 300), rice_factor_dBm(line_index))
%set(gca, 'XDir','reverse')
xlabel('Distance from basestation [log(m)]')
ylabel('Rice Factor [dBm]')
title('Rice Factor between basestation and bottom wall')

%% plotting delays

figure
scatter(all_rx_xy(1, :), all_rx_xy(2, :), [], delay_spread, 'filled')
hold on
lines = line(x, y, 'Color', 'black');
xlabel('width [m]')
ylabel('height [m]')
axis equal
c = colorbar;
c.Label.String = "time [s]";
title('Delay Spread as a heatmap')

%plot BS to bottom wall
figure
plot(log10(all_rx_xy(2,line_index) - 300), delay_spread(line_index))
%set(gca, 'XDir','reverse')
xlabel('Distance from basestation [log(m)]')
ylabel('time [s]')
title('Delay Spread between basestation and bottom wall')

%% raytracing for wanted position with lines

figure
hold on
lines = line(x, y, 'Color', 'black');
axis('equal')
xlabel('width [m]')
ylabel('height [m]')
title('Ray tracing for one position')
plot([tx_x ray_traced_position(1)], [tx_y ray_traced_position(2)] , '*')

ray_tracing_index = find(and(all_rx_xy(1, :) == ray_traced_position(1), all_rx_xy(2, :) == ray_traced_position(2)));
drawRayTracingLines(x, y, all_rx_xy, tx_x, tx_y, no_LOS, ray_traced_position, impulse_response, ray_tracing_index);

%% plot impulse response

impulse_index = find(and(all_rx_xy(1, :) == impulse_position(1), all_rx_xy(2, :) == impulse_position(2)));

single_impulse_response = nonzeros((impulse_response(impulse_index, :, :)));
single_impulse_response = reshape(single_impulse_response, 2, length(single_impulse_response)/2);

%fysical impulse response
figure
hold on
stem(single_impulse_response(2, :), abs(single_impulse_response(1, :)))
xlabel('delay [s]')
ylabel('Amplitude [V]')
title('Channel impulse response')

%tapped delay line (assume unccorrelated scattering)
amount_of_taps = ceil(max(single_impulse_response(2,:))/time_resolution);
taps = (1:amount_of_taps)*time_resolution;
if(delay_spread(impulse_index) < time_resolution/10)
    disp('narrowband')
    figure
    hold on
    stem(taps(find(mean(single_impulse_response(2, :)) <= taps, 1, 'first')), abs(sum(single_impulse_response(1, :))))
    xlabel('delay [s]')
    ylabel('Amplitude [V]')
    title('Channel impulse response narrow band (Tapped delay line)')
else
    disp('wideband')
    taps_amplitude = zeros(size(taps));
    for i = 1:length(single_impulse_response(1,:))
        j = 1;
        while(single_impulse_response(2,i) > taps(j))
            j = j + 1;
        end
        taps_amplitude(j) = taps_amplitude(j) + single_impulse_response(1,i);
    end

    figure
    hold on
    stem(taps,abs(taps_amplitude))
    xlabel('delay [s]')
    ylabel('Amplitude [V]')
    title('Channel impulse response wideband (Tapped delay line)')
end









    
     


