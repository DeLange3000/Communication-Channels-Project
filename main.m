close all;
clear;
clc;

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
for i = 1:length(all_rx_xy(1,:))
    Voc = 0;
    %check if there is a LOS
    if(isempty(find(no_LOS(3,:) == i))) %true is there is LOS
        r_los = sqrt((all_rx_xy(1,i) - tx_x)^2 + (all_rx_xy(2,i) - tx_y)^2); %since same height
        E_los = sqrt(60*EIRP)/r_los*exp(-1i*2*pi*Fc/c*r_los);
        Voc = Voc + E_los*effective_h;

        %if LOS then also ground reflection
        r_gr = sqrt(r_los^2 + (2*h_bs)^2);
        theta = acos(h_bs/(r_gr/2));
        effective_h_gr = -(c/Fc)/pi*cos(pi/2*cos(theta))/sin(theta)^2;
        lamba_parallel = (cos(theta) - sqrt(1/permitivity)*sqrt(1-1/permitivity*sin(theta)^2))/(cos(theta) + sqrt(1/permitivity)*sqrt(1-1/permitivity*sin(theta)^2));
        E_gr = lamba_parallel*sqrt(60*EIRP)/r_gr*exp(-1i*2*pi*Fc/c*r_gr);
        Voc = Voc + E_gr*effective_h_gr;
    else

        %diffraction!!
    end

    % one reflection
    % only code for this specific geometry to reduce complexity
    % how: mirror transmitter to left and right and bottom and for every street. Draw line from BS to UE. if there is
    % only one intersection with a wall -> calc reflection

    % mirror BS to left

    tx_y_mirror = tx_y;
    tx_x_mirror = tx_x + 40;
    [amount_of_intersec, intersections] = intersectionCalculator(x, y, all_rx_xy(:,i), tx_x_mirror, tx_y_mirror, 1);
    

    %visualize intersections
    if (amount_of_intersec == 1)
            if(intersections(1) == 40)
            figure
            hold on
            lines = line(x, y, 'Color', 'black');
            axis('equal')
            xlabel('width [m]')
            ylabel('height [m]')
            title('Positions where there is no line of sight')
            plot([tx_x intersections(1) all_rx_xy(1,i)], [tx_y intersections(2), all_rx_xy(2,i)]);
            end
    end




    % mirror BS to right


    % mirror BS to bottom


    % mirror BS to Rue du grand hospice


    % mirror BS to Rue du rouleau


    % mirror BS to Rue du Peuplier



    %check if wall is horizontal or vertical to get theta (same for 2
    %reflections?)
    










    % turn E_field into received power using effective height of dipole to get
    % Voc and then Ra (see exercise slides 25 -> 19)
    received_power(i) = 1/(8*Ra)*abs(Voc)^2;
end

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
xlabel('width [m]')
ylabel('height [m]')
axis equal
c = colorbar;
c.Label.String = "Received power [dBm]";
title('Received power as a heatmap')

%plot BS to bottom wall
line_index = find(all_rx_xy(1,:) == 20.5); %find elements on line between BS and bottom wall
figure
plot(log10(all_rx_xy(2,line_index) - 300), received_power_dBm(line_index), '*')
%set(gca, 'XDir','reverse')
xlabel('Distance from basestation [log(m)]')
ylabel('Received power [dBm]')
title('Received power between basestation and bottom wall')


%% raytracing for wanted position with lines


    
     


