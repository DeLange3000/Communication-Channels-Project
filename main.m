close all;
clear;
clc;

%% enables

enable_LOS = 0;
enable_ground_reflections = 0;
enable_diffraction = 0;
enable_one_reflection = 0;
enable_left_wall_reflection = 0;
enable_right_wall_reflection = 0;
enable_bottom_wall_reflection = 0;
eanble_rue_du_grand_hospice = 0;
enable_rue_du_rouleau = 0;
enable_rue_du_peuplier = 0;
enable_double_reflections = 1;
enable_left_and_right_wall = 0;
enable_left_and_bottom_wall = 0;
enable_left_and_streets = 0;
enable_right_and_left_wall = 0;
enable_right_and_bottom_wall = 0;
enable_bottom_and_left_wall = 0;
enable_bottom_and_right_wall = 1;
enable_bottom_and_streets = 0;


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
        if(enable_LOS == 1)
            r_los = sqrt((all_rx_xy(1,i) - tx_x)^2 + (all_rx_xy(2,i) - tx_y)^2); %since same height
            E_los = sqrt(60*EIRP)/r_los*exp(-1i*2*pi*Fc/c*r_los);
            Voc = Voc + E_los*effective_h;
        end

        %if LOS then also ground reflection
        if(enable_ground_reflections == 1)
            r_gr = sqrt(r_los^2 + (2*h_bs)^2);
            theta = acos(h_bs/(r_gr/2));
            effective_h_gr = -(c/Fc)/pi*cos(pi/2*cos(theta))/sin(theta)^2;
            lamba_parallel = (cos(theta) - sqrt(1/permitivity)*sqrt(1-1/permitivity*sin(theta)^2))/(cos(theta) + sqrt(1/permitivity)*sqrt(1-1/permitivity*sin(theta)^2));
            E_gr = lamba_parallel*sqrt(60*EIRP)/r_gr*exp(-1i*2*pi*Fc/c*r_gr);
            Voc = Voc + E_gr*effective_h_gr;
        end
    else
        %diffraction!!
        if(enable_diffraction == 1)

        end
    end

    % one reflection
    % only code for this specific geometry to reduce complexity
    % how: mirror transmitter to left and right and bottom and for every street. Draw line from BS to UE. if there is
    % only one intersection with a wall -> calc reflection

    if(enable_one_reflection == 1)
        if(enable_left_wall_reflection == 1)
            % mirror BS to left
            tx_y_mirror = tx_y;
            tx_x_mirror = -tx_x;
            [amount_of_intersec, intersections] = intersectionCalculator(x, y, all_rx_xy(:,i), tx_x_mirror, tx_y_mirror, 1);
    
            %remove any reflections from the left street
            j = 1;
            while j <= amount_of_intersec
                if(intersections(1,j) < 0)
                    intersections(:,j) = [];
                    amount_of_intersec = amount_of_intersec - 1;
                else
                    j = j + 1;
                end
            end
            if(amount_of_intersec == 1)
    
    %             figure
    %             hold on
    %             lines = line(x, y, 'Color', 'black');
    %             axis('equal')
    %             xlabel('width [m]')
    %             ylabel('height [m]')
    %             title('Positions where there is no line of sight')
    %             plot([tx_x intersections(1) all_rx_xy(1,i)], [tx_y intersections(2), all_rx_xy(2,i)]);
    
                r_refl = sqrt((all_rx_xy(1,i) - tx_x_mirror)^2 + (all_rx_xy(2,i) - tx_y_mirror)^2);
                theta = atan((tx_y_mirror - intersections(2))/(intersections(1) - tx_x_mirror));
                lamba_ortho = (cos(theta) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta)^2))/(cos(theta) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta)^2));
                E_refl = lamba_ortho*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                Voc = Voc + E_refl*effective_h;
            end
        end
    
        if(enable_right_wall_reflection == 1)    

            % mirror BS to right
            tx_y_mirror = tx_y;
            tx_x_mirror = tx_x + 40;
            [amount_of_intersec, intersections] = intersectionCalculator(x, y, all_rx_xy(:,i), tx_x_mirror, tx_y_mirror, 1);
            %exclude intersections with vertical right streets
            j = 1;
            while j <= amount_of_intersec
                if(intersections(1,j) > 40)
                    intersections(:,j) = [];
                    amount_of_intersec = amount_of_intersec - 1;
                else
                    j = j + 1;
                end
            end
    
            if(amount_of_intersec == 1)
                r_refl = sqrt((all_rx_xy(1,i) - tx_x_mirror)^2 + (all_rx_xy(2,i) - tx_y_mirror)^2);
                theta = atan((tx_y_mirror - intersections(2))/(tx_x_mirror - intersections(1)));
                lamba_ortho = (cos(theta) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta)^2))/(cos(theta) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta)^2));
                E_refl = lamba_ortho*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                Voc = Voc + E_refl*effective_h;
            end
        end
    
        if(enable_bottom_wall_reflection == 1)
        % mirror BS to bottom
    
            tx_y_mirror = -tx_y;
            tx_x_mirror = tx_x;
            [amount_of_intersec, intersections] = intersectionCalculator(x, y, all_rx_xy(:,i), tx_x_mirror, tx_y_mirror, 1);
    
            if(amount_of_intersec == 1)
                %determine if reflection can reach point (streets will
                %otherwise receive signal while they cannot)

                [a, b] = intersectionCalculator(x, y, intersections, tx_x, tx_y, 1);
                if(a <= 1)
                    r_refl = sqrt((all_rx_xy(1,i) - tx_x_mirror)^2 + (all_rx_xy(2,i) - tx_y_mirror)^2);
                    theta = atan(abs(abs(tx_y_mirror) - abs(intersections(2)))/abs(tx_x_mirror - abs(intersections(1))));
                    lamba_ortho = (cos(theta) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta)^2))/(cos(theta) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta)^2));
                    E_refl = lamba_ortho*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                    Voc = Voc + E_refl*effective_h;
                end
            end
        end
    
        % mirror BS to Rue du grand hospice
        if(eanble_rue_du_grand_hospice == 1)
            tx_x_mirror = tx_x;
            tx_y_mirror = tx_y - 60;
            [amount_of_intersec, intersections] = intersectionCalculator(x, y, all_rx_xy(:,i), tx_x_mirror, tx_y_mirror, 1);
            if(amount_of_intersec == 2 && intersections(1,1) == 40 && intersections(2,2) == 270)
                %determine if reflection can reach point (streets will
                %otherwise receive signal while they cannot)
                [a, b] = intersectionCalculator(x, y, intersections(:,2), tx_x, tx_y, 1);
                if(a <= 1)
                    r_refl = sqrt((all_rx_xy(1,i) - tx_x_mirror)^2 + (all_rx_xy(2,i) - tx_y_mirror)^2);
                    theta = atan(abs(abs(tx_y_mirror) - abs(intersections(2)))/abs(tx_x_mirror - abs(intersections(1))));
                    lamba_ortho = (cos(theta) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta)^2))/(cos(theta) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta)^2));
                    E_refl = lamba_ortho*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                    Voc = Voc + E_refl*effective_h;
                end

            end

        end
    
        % mirror BS to Rue du rouleau
        if(enable_rue_du_rouleau == 1)
            tx_x_mirror = tx_x;
            tx_y_mirror = tx_y - 240;
            [amount_of_intersec, intersections] = intersectionCalculator(x, y, all_rx_xy(:,i), tx_x_mirror, tx_y_mirror, 1);
            if(amount_of_intersec == 2 && intersections(1,1) == 40 && intersections(2,2) == 180)
                %determine if reflection can reach point (streets will
                %otherwise receive signal while they cannot)
                [a, b] = intersectionCalculator(x, y, intersections(:,2), tx_x, tx_y, 1);
                if(a <= 1)
                    r_refl = sqrt((all_rx_xy(1,i) - tx_x_mirror)^2 + (all_rx_xy(2,i) - tx_y_mirror)^2);
                    theta = atan(abs(abs(tx_y_mirror) - abs(intersections(2)))/abs(tx_x_mirror - abs(intersections(1))));
                    lamba_ortho = (cos(theta) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta)^2))/(cos(theta) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta)^2));
                    E_refl = lamba_ortho*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                    Voc = Voc + E_refl*effective_h;
                end

            end

        end
    
    
        % mirror BS to Rue du Peuplier
        if(enable_rue_du_peuplier == 1)
            tx_x_mirror = tx_x;
            tx_y_mirror = tx_y - 460;
            [amount_of_intersec, intersections] = intersectionCalculator(x, y, all_rx_xy(:,i), tx_x_mirror, tx_y_mirror, 1);
            if(amount_of_intersec == 3 && intersections(2,1) == 0 && intersections(1,2) == 40 && intersections(2,3) == 70) %also intersects with bottom wall
                %determine if reflection can reach point (streets will
                %otherwise receive signal while they cannot)
                [a, b] = intersectionCalculator(x, y, intersections(:,3), tx_x, tx_y, 1);
                if(a <= 1)
                    r_refl = sqrt((all_rx_xy(1,i) - tx_x_mirror)^2 + (all_rx_xy(2,i) - tx_y_mirror)^2);
                    theta = atan(abs(abs(tx_y_mirror) - abs(intersections(2)))/abs(tx_x_mirror - abs(intersections(1))));
                    lamba_ortho = (cos(theta) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta)^2))/(cos(theta) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta)^2));
                    E_refl = lamba_ortho*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                    Voc = Voc + E_refl*effective_h;
                end

            end

        end
    end


        % double reflections
        if(enable_double_reflections == 1)

            %bottom wall + right wall
            if(enable_left_and_right_wall == 1)
                tx_x_mirror1 = tx_x - 20;
                tx_y_mirror1 = tx_y;
                tx_x_mirror2 = tx_x + 60;
                tx_y_mirror2 = tx_y;
                [amount_of_intersec, intersections] = intersectionCalculator(x, y, all_rx_xy(:,i), tx_x_mirror2, tx_y_mirror2, 1);
                %exclude intersections with vertical right streets
                j = 1;
                while j <= amount_of_intersec
                    if(intersections(1,j) > 40)
                        intersections(:,j) = [];
                        amount_of_intersec = amount_of_intersec - 1;
                    else
                        j = j + 1;
                    end
                end
        
                if(amount_of_intersec == 1)
                    r_refl = sqrt((all_rx_xy(1,i) - tx_x_mirror2)^2 + (all_rx_xy(2,i) - tx_y_mirror2)^2);
                    theta2 = atan((tx_y_mirror2 - intersections(2))/(tx_x_mirror2 - intersections(1)));
                    theta1 = theta2;
                    lamba_ortho1 = (cos(theta1) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2))/(cos(theta1) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2));
                    lamba_ortho2 = (cos(theta2) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2))/(cos(theta2) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2));
                    E_refl = lamba_ortho1*lamba_ortho2*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                    Voc = Voc + E_refl*effective_h;
                end
            end

            %left and bottom wall
            if(enable_left_and_bottom_wall == 1)
                tx_x_mirror1 = tx_x - 40;
                tx_y_mirror1 = tx_y;
                tx_x_mirror2 = tx_x_mirror1;
                tx_y_mirror2 = -tx_y_mirror1;
                [amount_of_intersec, intersections] = intersectionCalculator(x, y, all_rx_xy(:,i), tx_x_mirror2, tx_y_mirror2, 1);
              
                if(amount_of_intersec == 1 && intersections(2,1) == 0)
                    [amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections, tx_x_mirror1, tx_y_mirror1, 1);
                    [a, b] = intersectionCalculator(x, y, intersections, tx_x, tx_y, 1);
                    if(a <= 1 && isempty(intersections_mirror1(intersections_mirror1 == 10)))
                        r_refl = sqrt((all_rx_xy(1,i) - tx_x_mirror2)^2 + (all_rx_xy(2,i) - tx_y_mirror2)^2);
                        theta2 = atan((tx_y_mirror2 - intersections(2))/(tx_x_mirror2 - intersections(1)));
                        theta1 = pi/2 - theta2;
                        lamba_ortho1 = (cos(theta1) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2))/(cos(theta1) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2));
                        lamba_ortho2 = (cos(theta2) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2))/(cos(theta2) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2));
                        E_refl = lamba_ortho1*lamba_ortho2*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                        Voc = Voc + E_refl*effective_h;
                    end
                end
            end

            if(enable_left_and_streets == 1)

                %rue du grand hospice
                tx_x_mirror1 = tx_x - 40;
                tx_y_mirror1 = tx_y;
                tx_x_mirror2 = tx_x_mirror1;
                tx_y_mirror2 = tx_y - 60;
                [amount_of_intersec, intersections] = intersectionCalculator(x, y, all_rx_xy(:,i), tx_x_mirror2, tx_y_mirror2, 1);
              
                if(amount_of_intersec == 3 && intersections(1,2) == 40 && intersections(2,3) == 270)
                    [amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections(:,3), tx_x_mirror1, tx_y_mirror1, 1);
                    [a, b] = intersectionCalculator(x, y, intersections, tx_x, tx_y, 1);
                    if(a <= 1 && amount_of_intersec <= 2)
                        r_refl = sqrt((all_rx_xy(1,i) - tx_x_mirror2)^2 + (all_rx_xy(2,i) - tx_y_mirror2)^2);
                        theta2 = atan((tx_y_mirror2 - intersections(2))/(tx_x_mirror2 - intersections(1)));
                        theta1 = pi/2 - theta2;
                        lamba_ortho1 = (cos(theta1) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2))/(cos(theta1) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2));
                        lamba_ortho2 = (cos(theta2) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2))/(cos(theta2) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2));
                        E_refl = lamba_ortho1*lamba_ortho2*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                        Voc = Voc + E_refl*effective_h;
                    end
                end

                %rue du rouleau
                tx_x_mirror1 = tx_x - 40;
                tx_y_mirror1 = tx_y;
                tx_x_mirror2 = tx_x_mirror1;
                tx_y_mirror2 = tx_y - 240;
                [amount_of_intersec, intersections] = intersectionCalculator(x, y, all_rx_xy(:,i), tx_x_mirror2, tx_y_mirror2, 1);
              
                if(amount_of_intersec == 3 && intersections(1,2) == 40 && intersections(2,3) == 180)
                    [amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections(:,3), tx_x_mirror1, tx_y_mirror1, 1);
                    [a, b] = intersectionCalculator(x, y, intersections, tx_x, tx_y, 1);
                    if(a <= 1 && amount_of_intersec <= 2)
                        r_refl = sqrt((all_rx_xy(1,i) - tx_x_mirror2)^2 + (all_rx_xy(2,i) - tx_y_mirror2)^2);
                        theta2 = atan((tx_y_mirror2 - intersections(2))/(tx_x_mirror2 - intersections(1)));
                        theta1 = pi/2 - theta2;
                        lamba_ortho1 = (cos(theta1) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2))/(cos(theta1) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2));
                        lamba_ortho2 = (cos(theta2) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2))/(cos(theta2) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2));
                        E_refl = lamba_ortho1*lamba_ortho2*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                        Voc = Voc + E_refl*effective_h;
                    end
                end

                %rue du peuplier
                tx_x_mirror1 = tx_x - 40;
                tx_y_mirror1 = tx_y;
                tx_x_mirror2 = tx_x_mirror1;
                tx_y_mirror2 = tx_y - 460;
                [amount_of_intersec, intersections] = intersectionCalculator(x, y, all_rx_xy(:,i), tx_x_mirror2, tx_y_mirror2, 1);
              
                if(amount_of_intersec == 3 && intersections(1,2) == 40 && intersections(2,3) == 70)
                    [amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections(:,3), tx_x_mirror1, tx_y_mirror1, 1);
                    [a, b] = intersectionCalculator(x, y, intersections, tx_x, tx_y, 1);
                    if(a <= 1 && amount_of_intersec <= 2)
                        r_refl = sqrt((all_rx_xy(1,i) - tx_x_mirror2)^2 + (all_rx_xy(2,i) - tx_y_mirror2)^2);
                        theta2 = atan((tx_y_mirror2 - intersections(2))/(tx_x_mirror2 - intersections(1)));
                        theta1 = pi/2 - theta2;
                        lamba_ortho1 = (cos(theta1) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2))/(cos(theta1) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2));
                        lamba_ortho2 = (cos(theta2) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2))/(cos(theta2) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2));
                        E_refl = lamba_ortho1*lamba_ortho2*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                        Voc = Voc + E_refl*effective_h;
                    end
                end
            end

            %right wall and left wall
            if(enable_right_and_left_wall == 1)
                tx_x_mirror1 = tx_x + 40;
                tx_y_mirror1 = tx_y;
                tx_x_mirror2 = tx_x - 60;
                tx_y_mirror2 = tx_y;
                [amount_of_intersec, intersections] = intersectionCalculator(x, y, all_rx_xy(:,i), tx_x_mirror2, tx_y_mirror2, 1);
                %exclude intersections with vertical right streets
                j = 1;
                while j <= amount_of_intersec
                    if(intersections(1,j) > 40 || intersections(1,j) < 0)
                        intersections(:,j) = [];
                        amount_of_intersec = amount_of_intersec - 1;
                    else
                        j = j + 1;
                    end
                end
        
                if(amount_of_intersec == 1)
                    r_refl = sqrt((all_rx_xy(1,i) - tx_x_mirror2)^2 + (all_rx_xy(2,i) - tx_y_mirror2)^2);
                    theta2 = atan((tx_y_mirror2 - intersections(2))/(tx_x_mirror2 - intersections(1)));
                    theta1 = theta2;
                    lamba_ortho1 = (cos(theta1) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2))/(cos(theta1) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2));
                    lamba_ortho2 = (cos(theta2) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2))/(cos(theta2) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2));
                    E_refl = lamba_ortho1*lamba_ortho2*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                    Voc = Voc + E_refl*effective_h;
                end
            end
                
            % right and bottom wall

            if(enable_right_and_bottom_wall == 1)
                tx_x_mirror1 = tx_x + 40;
                tx_y_mirror1 = tx_y;
                tx_x_mirror2 = tx_x_mirror1;
                tx_y_mirror2 = -tx_y_mirror1;
                [amount_of_intersec, intersections] = intersectionCalculator(x, y, all_rx_xy(:,i), tx_x_mirror2, tx_y_mirror2, 1);
              
                if(amount_of_intersec == 1 && intersections(2,1) == 0)
                    [amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections, tx_x_mirror1, tx_y_mirror1, 1);
                    [a, b] = intersectionCalculator(x, y, intersections, tx_x, tx_y, 1);
                    if(a <= 1 && isempty(intersections_mirror1(intersections_mirror1 == 10)))
                        r_refl = sqrt((all_rx_xy(1,i) - tx_x_mirror2)^2 + (all_rx_xy(2,i) - tx_y_mirror2)^2);
                        theta2 = atan((tx_y_mirror2 - intersections(2))/(tx_x_mirror2 - intersections(1)));
                        theta1 = pi/2 - theta2;
                        lamba_ortho1 = (cos(theta1) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2))/(cos(theta1) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2));
                        lamba_ortho2 = (cos(theta2) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2))/(cos(theta2) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2));
                        E_refl = lamba_ortho1*lamba_ortho2*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                        Voc = Voc + E_refl*effective_h;
                    end
                end
            end

            %bottom wall to left wall
            if(enable_bottom_and_left_wall == 1)
                tx_x_mirror1 = tx_x;
                tx_y_mirror1 = - tx_y;
                tx_x_mirror2 = tx_x_mirror1 - 40;
                tx_y_mirror2 = tx_y_mirror1;
                [amount_of_intersec, intersections] = intersectionCalculator(x, y, all_rx_xy(:,i), tx_x_mirror2, tx_y_mirror2, 1);
              
                if(amount_of_intersec == 3 && intersections(2,1) == 0 && intersections(1, 2) == 0)
                    [amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections, tx_x_mirror1, tx_y_mirror1, 1);
                    r_refl = sqrt((all_rx_xy(1,i) - tx_x_mirror2)^2 + (all_rx_xy(2,i) - tx_y_mirror2)^2);

                    theta2 = atan((tx_y_mirror2 - intersections(2))/(tx_x_mirror2 - intersections(1)));
                    theta1 = pi/2 - theta2;
                    lamba_ortho1 = (cos(theta1) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2))/(cos(theta1) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2));
                    lamba_ortho2 = (cos(theta2) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2))/(cos(theta2) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2));
                    E_refl = lamba_ortho1*lamba_ortho2*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                    Voc = Voc + E_refl*effective_h;
                end
            end

            %bottom wall right wall
            if(enable_bottom_and_right_wall == 1)
                

            end
        








        end



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
hold on
lines = line(x, y, 'Color', 'black');
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


    
     


