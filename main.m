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

ray_traced_position = [42.5 ; 188.5];



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


if(isempty(all_rx_xy(all_rx_xy == ray_traced_position)))
    disp('Not a valid raytracing position')
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
            r_los = sqrt((all_rx_xy(1,i) - tx_x)^2 + (all_rx_xy(2,i) - tx_y)^2); %since same height
            r_gr = sqrt(r_los^2 + (2*h_bs)^2);
            theta = acos(h_bs/(r_gr/2));
            effective_h_gr = -(c/Fc)/pi*cos(pi/2*cos(theta))/sin(theta)^2;
            lamba_parallel = (cos(theta) - sqrt(1/permitivity)*sqrt(1-1/permitivity*sin(theta)^2))/(cos(theta) + sqrt(1/permitivity)*sqrt(1-1/permitivity*sin(theta)^2));
            E_gr = lamba_parallel*sqrt(60*EIRP)/r_gr*exp(-1i*2*pi*Fc/c*r_gr);
            Voc = Voc + E_gr*effective_h_gr;
        end
    else
        %diffraction!!
        %get right diffraction point
        if(enable_diffraction == 1)
            if(all_rx_xy(2,i) > 270)
                diff_point = [40, 280];
            elseif(all_rx_xy(2,i) > 180)
                diff_point = [40, 190];
            elseif(all_rx_xy(2,i) > 70)
                diff_point = [40, 80];
            elseif(all_rx_xy(1,i) > 40)
                diff_point = [40, 10];
            else
                diff_point = [0, 10];
            end
            
            r_los = sqrt((all_rx_xy(1,i) - tx_x)^2 + (all_rx_xy(2,i) - tx_y)^2); %since same height
            r_diff = sqrt((all_rx_xy(1,i) - diff_point(1))^2 + (all_rx_xy(2,i) - diff_point(2))^2) + sqrt((tx_x - diff_point(1))^2 + (tx_y - diff_point(2))^2);
            v = sqrt(2/pi*beta*(r_diff - r_los));
            ampl_F = sqrt(10^((-6.9 -20*log10(sqrt((v - 0.1)^2 + 1) + v -0.1))/10));
            angle_F = -pi/4 - pi/2*v^2;

            E_diff = ampl_F*sqrt(60*EIRP)/r_los*exp(-1i*2*pi*Fc/c*r_los + angle_F*1i);
            Voc = Voc + E_diff*effective_h;

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

            %left wall + right wall
            if(enable_left_and_right_wall == 1)
                tx_x_mirror1 = tx_x - 40;
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

                    theta2 = atan((tx_y_mirror2 - intersections(2, 3))/(tx_x_mirror2 - intersections(1, 3)));
                    theta1 = pi/2 - theta2;
                    lamba_ortho1 = (cos(theta1) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2))/(cos(theta1) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2));
                    lamba_ortho2 = (cos(theta2) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2))/(cos(theta2) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2));
                    E_refl = lamba_ortho1*lamba_ortho2*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                    Voc = Voc + E_refl*effective_h;
                end
            end

            %bottom wall right wall
            if(enable_bottom_and_right_wall == 1)
                tx_x_mirror1 = tx_x;
                tx_y_mirror1 = - tx_y;
                tx_x_mirror2 = tx_x_mirror1 + 40;
                tx_y_mirror2 = tx_y_mirror1;
                [amount_of_intersec, intersections] = intersectionCalculator(x, y, all_rx_xy(:,i), tx_x_mirror2, tx_y_mirror2, 1);

                %exclude intersections with vertical right streets
                j = 1;
                while j <= amount_of_intersec
                    if(intersections(1,j) > 40 || intersections(2, j) == 0)
                        intersections(:,j) = [];
                        amount_of_intersec = amount_of_intersec - 1;
                    else
                        j = j + 1;
                    end
                end
              
                if(amount_of_intersec == 1)
                    %[amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections(:,2), tx_x_mirror1, tx_y_mirror1, 1);
                    r_refl = sqrt((all_rx_xy(1,i) - tx_x_mirror2)^2 + (all_rx_xy(2,i) - tx_y_mirror2)^2);

                    theta2 = atan((tx_y_mirror2 - intersections(2))/(tx_x_mirror2 - intersections(1)));
                    theta1 = pi/2 - theta2;
                    lamba_ortho1 = (cos(theta1) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2))/(cos(theta1) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2));
                    lamba_ortho2 = (cos(theta2) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2))/(cos(theta2) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2));
                    E_refl = lamba_ortho1*lamba_ortho2*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                    Voc = Voc + E_refl*effective_h;
                end
            end

            %bottom wall and street
            if(enable_bottom_and_streets == 1)

                % bottom wall and top of sint katelijneplijn left and right
                tx_x_mirror1 = tx_x;
                tx_y_mirror1 = - tx_y;
                tx_x_mirror2 = tx_x_mirror1;
                tx_y_mirror2 = - tx_y_mirror1 + 20;
                [amount_of_intersec, intersections] = intersectionCalculator(x, y, all_rx_xy(:,i), tx_x_mirror2, tx_y_mirror2, 1);

                %exclude intersections with vertical right streets
                j = 1;
                while j <= amount_of_intersec
                    if(intersections(2, j) > 10)
                        intersections(:,j) = [];
                        amount_of_intersec = amount_of_intersec - 1;
                    else
                        j = j + 1;
                    end
                end
              
                if(amount_of_intersec == 1)
                    %[amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections(:,2), tx_x_mirror1, tx_y_mirror1, 1);
                    r_refl = sqrt((all_rx_xy(1,i) - tx_x_mirror2)^2 + (all_rx_xy(2,i) - tx_y_mirror2)^2);
                    [amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections, tx_x_mirror1, tx_y_mirror1, 1);
                    [amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections_mirror1, tx_x, tx_y, 1);
                    if(amount_of_intersec <= 1)

                        theta2 = atan((tx_y_mirror2 - intersections(2))/(tx_x_mirror2 - intersections(1)));
                        theta1 = theta2;
                        lamba_ortho1 = (cos(theta1) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2))/(cos(theta1) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2));
                        lamba_ortho2 = (cos(theta2) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2))/(cos(theta2) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2));
                        E_refl = lamba_ortho1*lamba_ortho2*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                        Voc = Voc + E_refl*effective_h;
                    end
                end

                % bottom wall and top of rue du peuplier
                tx_x_mirror1 = tx_x;
                tx_y_mirror1 = - tx_y;
                tx_x_mirror2 = tx_x_mirror1;
                tx_y_mirror2 = - tx_y_mirror1 + 160;
                [amount_of_intersec, intersections] = intersectionCalculator(x, y, all_rx_xy(:,i), tx_x_mirror2, tx_y_mirror2, 1);

                %exclude intersections with vertical right streets
                j = 1;
                while j <= amount_of_intersec
                    if(intersections(2, j) > 80)
                        intersections(:,j) = [];
                        amount_of_intersec = amount_of_intersec - 1;
                    else
                        j = j + 1;
                    end
                end
              
                if(amount_of_intersec == 1)
                    %[amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections(:,2), tx_x_mirror1, tx_y_mirror1, 1);
                    r_refl = sqrt((all_rx_xy(1,i) - tx_x_mirror2)^2 + (all_rx_xy(2,i) - tx_y_mirror2)^2);
                    [amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections, tx_x_mirror1, tx_y_mirror1, 1);
                    %only keep intersections with street and right wall
                    j = 1;
                    while j <= amount_of_intersec
                        if((intersections_mirror1(2, j) > 60 || intersections_mirror1(2, j) <= 0) && intersections_mirror1(1, j) > 0 )
                            intersections_mirror1(:,j) = [];
                            amount_of_intersec = amount_of_intersec - 1;
                        else
                            j = j + 1;
                        end
                    end
                    if(amount_of_intersec == 0)
                        theta2 = atan((tx_y_mirror2 - intersections(2))/(tx_x_mirror2 - intersections(1)));
                        theta1 = theta2;
                        lamba_ortho1 = (cos(theta1) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2))/(cos(theta1) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2));
                        lamba_ortho2 = (cos(theta2) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2))/(cos(theta2) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2));
                        E_refl = lamba_ortho1*lamba_ortho2*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                        Voc = Voc + E_refl*effective_h;
                    end
                end

                % bottom wall and top of rue du rouleau
                % bottom wall and top of rue du grand hospice
                % both have no reflections from bottom wall since for one
                % reflection the bottom wall does not reach these streets
            end

            %double reflection inside streets (excluding bottom streets
            %since they are included in bottom wall + streets double reflections)
            if(enable_streets_and_streets == 1)

                %rue du peuplier
                tx_x_mirror1 = tx_x;
                tx_y_mirror1 = - tx_y + 140;
                tx_x_mirror2 = tx_x_mirror1;
                tx_y_mirror2 = tx_y + 20;
                [amount_of_intersec, intersections] = intersectionCalculator(x, y, all_rx_xy(:,i), tx_x_mirror2, tx_y_mirror2, 1);

                %exclude intersections with vertical right streets
                j = 1;
                while j <= amount_of_intersec
                    if(intersections(2, j) > 80)
                        intersections(:,j) = [];
                        amount_of_intersec = amount_of_intersec - 1;
                    else
                        j = j + 1;
                    end
                end
              
                if(amount_of_intersec == 1)
                    %[amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections(:,2), tx_x_mirror1, tx_y_mirror1, 1);
                    r_refl = sqrt((all_rx_xy(1,i) - tx_x_mirror2)^2 + (all_rx_xy(2,i) - tx_y_mirror2)^2);
                    [amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections, tx_x_mirror1, tx_y_mirror1, 1);
                    %only keep intersections with street and right wall
                    j = 1;
                    while j <= amount_of_intersec
                        if(intersections_mirror1(2, j) ~= 70)
                            intersections_mirror1(:,j) = [];
                            amount_of_intersec = amount_of_intersec - 1;
                        else
                            j = j + 1;
                        end
                    end
                    if(amount_of_intersec == 1)
                        [amount_of_intersec, intersections_LOS] = intersectionCalculator(x, y, intersections_mirror1, tx_x, tx_y, 1);
                        if(amount_of_intersec <= 1)
                            theta2 = atan((tx_y_mirror2 - intersections(2))/(tx_x_mirror2 - intersections(1)));
                            theta1 = theta2;
                            lamba_ortho1 = (cos(theta1) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2))/(cos(theta1) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2));
                            lamba_ortho2 = (cos(theta2) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2))/(cos(theta2) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2));
                            E_refl = lamba_ortho1*lamba_ortho2*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                            Voc = Voc + E_refl*effective_h;
                        end
                     end
                end

                %rue du rouleau
                tx_x_mirror1 = tx_x;
                tx_y_mirror1 = tx_y - 240;
                tx_x_mirror2 = tx_x_mirror1;
                tx_y_mirror2 = tx_y + 20;
                [amount_of_intersec, intersections] = intersectionCalculator(x, y, all_rx_xy(:,i), tx_x_mirror2, tx_y_mirror2, 1);

                %exclude intersections with vertical right streets
                j = 1;
                while j <= amount_of_intersec
                    if(intersections(2, j) > 190)
                        intersections(:,j) = [];
                        amount_of_intersec = amount_of_intersec - 1;
                    else
                        j = j + 1;
                    end
                end
              
                if(amount_of_intersec == 1)
                    %[amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections(:,2), tx_x_mirror1, tx_y_mirror1, 1);
                    r_refl = sqrt((all_rx_xy(1,i) - tx_x_mirror2)^2 + (all_rx_xy(2,i) - tx_y_mirror2)^2);
                    [amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections, tx_x_mirror1, tx_y_mirror1, 1);
                    %only keep intersections with street and right wall
                    j = 1;
                    while j <= amount_of_intersec
                        if(intersections_mirror1(2, j) ~= 180)
                            intersections_mirror1(:,j) = [];
                            amount_of_intersec = amount_of_intersec - 1;
                        else
                            j = j + 1;
                        end
                    end
                    if(amount_of_intersec == 1)
                        [amount_of_intersec, intersections_LOS] = intersectionCalculator(x, y, intersections_mirror1, tx_x, tx_y, 1);
                        if(amount_of_intersec <= 1)
                            theta2 = atan((tx_y_mirror2 - intersections(2))/(tx_x_mirror2 - intersections(1)));
                            theta1 = theta2;
                            lamba_ortho1 = (cos(theta1) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2))/(cos(theta1) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2));
                            lamba_ortho2 = (cos(theta2) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2))/(cos(theta2) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2));
                            E_refl = lamba_ortho1*lamba_ortho2*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                            Voc = Voc + E_refl*effective_h;
                        end
                     end
                end
                
                %rue du grand hospice
                tx_x_mirror1 = tx_x;
                tx_y_mirror1 = tx_y - 60;
                tx_x_mirror2 = tx_x_mirror1;
                tx_y_mirror2 = tx_y + 20;
                [amount_of_intersec, intersections] = intersectionCalculator(x, y, all_rx_xy(:,i), tx_x_mirror2, tx_y_mirror2, 1);

                %exclude intersections with vertical right streets
                j = 1;
                while j <= amount_of_intersec
                    if(intersections(2, j) > 280)
                        intersections(:,j) = [];
                        amount_of_intersec = amount_of_intersec - 1;
                    else
                        j = j + 1;
                    end
                end
              
                if(amount_of_intersec == 1)
                    %[amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections(:,2), tx_x_mirror1, tx_y_mirror1, 1);
                    r_refl = sqrt((all_rx_xy(1,i) - tx_x_mirror2)^2 + (all_rx_xy(2,i) - tx_y_mirror2)^2);
                    [amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections, tx_x_mirror1, tx_y_mirror1, 1);
                    %only keep intersections with street and right wall
                    j = 1;
                    while j <= amount_of_intersec
                        if(intersections_mirror1(2, j) ~= 270)
                            intersections_mirror1(:,j) = [];
                            amount_of_intersec = amount_of_intersec - 1;
                        else
                            j = j + 1;
                        end
                    end
                    if(amount_of_intersec == 1)
                        [amount_of_intersec, intersections_LOS] = intersectionCalculator(x, y, intersections_mirror1, tx_x, tx_y, 1);
                        if(amount_of_intersec <= 1)
                            theta2 = atan((tx_y_mirror2 - intersections(2))/(tx_x_mirror2 - intersections(1)));
                            theta1 = theta2;
                            lamba_ortho1 = (cos(theta1) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2))/(cos(theta1) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2));
                            lamba_ortho2 = (cos(theta2) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2))/(cos(theta2) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2));
                            E_refl = lamba_ortho1*lamba_ortho2*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                            Voc = Voc + E_refl*effective_h;
                        end
                     end
                end
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

figure
hold on
lines = line(x, y, 'Color', 'black');
axis('equal')
xlabel('width [m]')
ylabel('height [m]')
title('Ray tracing for one position')
plot([tx_x ray_traced_position(1)], [tx_y ray_traced_position(2)] , '*')


%check if there is a LOS
if(isempty(find(no_LOS([1 2],:) == ray_traced_position))) %true is there is LOS
plot([tx_x ray_traced_position(1)], [tx_y ray_traced_position(2)]);

else
%diffraction!!
%get right diffraction point
    if(ray_traced_position(2) > 270)
        diff_point = [40, 280];
    elseif(ray_traced_position(2) > 180)
        diff_point = [40, 190];
    elseif(ray_traced_position(2) > 70)
        diff_point = [40, 80];
    elseif(ray_traced_position(1) > 40)
        diff_point = [40, 10];
    else
        diff_point = [0, 10];
    end
    plot([tx_x diff_point(1) ray_traced_position(1)], [tx_y diff_point(2) ray_traced_position(2)]);

end

% mirror BS to left
tx_y_mirror = tx_y;
tx_x_mirror = -tx_x;
[amount_of_intersec, intersections] = intersectionCalculator(x, y, ray_traced_position, tx_x_mirror, tx_y_mirror, 1);

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
    plot([tx_x intersections(1) ray_traced_position(1)], [tx_y intersections(2), ray_traced_position(2)]);
end

% mirror BS to right
tx_y_mirror = tx_y;
tx_x_mirror = tx_x + 40;
[amount_of_intersec, intersections] = intersectionCalculator(x, y, ray_traced_position, tx_x_mirror, tx_y_mirror, 1);
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
    plot([tx_x intersections(1) ray_traced_position(1)], [tx_y intersections(2), ray_traced_position(2)]);
end

% mirror BS to bottom
tx_y_mirror = -tx_y;
tx_x_mirror = tx_x;
[amount_of_intersec, intersections] = intersectionCalculator(x, y, ray_traced_position, tx_x_mirror, tx_y_mirror, 1);

if(amount_of_intersec == 1)
    [a, b] = intersectionCalculator(x, y, intersections, tx_x, tx_y, 1);
    if(a <= 1)
        plot([tx_x intersections(1) ray_traced_position(1)], [tx_y intersections(2), ray_traced_position(2)]);
    end
end

% mirror BS to Rue du grand hospice
tx_x_mirror = tx_x;
tx_y_mirror = tx_y - 60;
[amount_of_intersec, intersections] = intersectionCalculator(x, y, ray_traced_position, tx_x_mirror, tx_y_mirror, 1);
if(amount_of_intersec == 2 && intersections(1,1) == 40 && intersections(2,2) == 270)
    [a, b] = intersectionCalculator(x, y, intersections(:,2), tx_x, tx_y, 1);
    if(a <= 1)
        plot([tx_x intersections(1, 2) ray_traced_position(1)], [tx_y intersections(2, 2), ray_traced_position(2)]);
    end
end

% mirror BS to Rue du rouleau
tx_x_mirror = tx_x;
tx_y_mirror = tx_y - 240;
[amount_of_intersec, intersections] = intersectionCalculator(x, y, ray_traced_position, tx_x_mirror, tx_y_mirror, 1);
if(amount_of_intersec == 2 && intersections(1,1) == 40 && intersections(2,2) == 180)
    [a, b] = intersectionCalculator(x, y, intersections(:,2), tx_x, tx_y, 1);
    if(a <= 1)
        plot([tx_x intersections(1, 2) ray_traced_position(1)], [tx_y intersections(2, 2), ray_traced_position(2)]);
    end
end


% mirror BS to Rue du Peuplier
tx_x_mirror = tx_x;
tx_y_mirror = tx_y - 460;
[amount_of_intersec, intersections] = intersectionCalculator(x, y, ray_traced_position, tx_x_mirror, tx_y_mirror, 1);
if(amount_of_intersec == 3 && intersections(2,1) == 0 && intersections(1,2) == 40 && intersections(2,3) == 70) %also intersects with bottom wall
    [a, b] = intersectionCalculator(x, y, intersections(:,3), tx_x, tx_y, 1);
    if(a <= 1)
        plot([tx_x intersections(1, 2) ray_traced_position(1)], [tx_y intersections(2, 2), ray_traced_position(2)]);
    end
end


% double reflections

%left wall + right wall
tx_x_mirror1 = tx_x - 40;
tx_y_mirror1 = tx_y;
tx_x_mirror2 = tx_x + 60;
tx_y_mirror2 = tx_y;
[amount_of_intersec, intersections] = intersectionCalculator(x, y, ray_traced_position, tx_x_mirror2, tx_y_mirror2, 1);
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
    [amount_of_intersec, intersections_double_refl] = intersectionCalculator(x, y, intersections, tx_x_mirror1, tx_y_mirror1, 1);
    j = 1;
    while j <= amount_of_intersec
        if(intersections_double_refl(1,j) ~= 0)
            intersections_double_refl(:,j) = [];
            amount_of_intersec = amount_of_intersec - 1;
        else
            j = j + 1;
        end
    end
    plot([tx_x intersections_double_refl(1) intersections(1) ray_traced_position(1)], [tx_y intersections_double_refl(2) intersections(2), ray_traced_position(2)]);
end

%left and bottom wall
tx_x_mirror1 = tx_x - 40;
tx_y_mirror1 = tx_y;
tx_x_mirror2 = tx_x_mirror1;
tx_y_mirror2 = -tx_y_mirror1;
[amount_of_intersec, intersections] = intersectionCalculator(x, y, ray_traced_position, tx_x_mirror2, tx_y_mirror2, 1);

if(amount_of_intersec == 1 && intersections(2,1) == 0)
    [amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections, tx_x_mirror1, tx_y_mirror1, 1);

    [a, b] = intersectionCalculator(x, y, intersections, tx_x, tx_y, 1);
    if(a <= 1 && isempty(intersections_mirror1(intersections_mirror1 == 10)))
        j = 1;
        while j <= amount_of_intersec
            if(intersections_mirror1(1,j) ~= 0)
                intersections_mirror1(:,j) = [];
                amount_of_intersec = amount_of_intersec - 1;
            else
                j = j + 1;
            end
        end
        plot([tx_x intersections_mirror1(1) intersections(1) ray_traced_position(1)], [tx_y intersections_mirror1(2) intersections(2), ray_traced_position(2)]);
    end
end

%left and rue du grand hospice
tx_x_mirror1 = tx_x - 40;
tx_y_mirror1 = tx_y;
tx_x_mirror2 = tx_x_mirror1;
tx_y_mirror2 = tx_y - 60;
[amount_of_intersec, intersections] = intersectionCalculator(x, y, ray_traced_position, tx_x_mirror2, tx_y_mirror2, 1);

if(amount_of_intersec == 3 && intersections(1,2) == 40 && intersections(2,3) == 270)
    [amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections(:,3), tx_x_mirror1, tx_y_mirror1, 1);

    [a, b] = intersectionCalculator(x, y, intersections, tx_x, tx_y, 1);
    if(a <= 1 && amount_of_intersec <= 2)
        j = 1;
        while j <= amount_of_intersec
            if(intersections_mirror1(1,j) ~= 0)
                intersections_mirror1(:,j) = [];
                amount_of_intersec = amount_of_intersec - 1;
            else
                j = j + 1;
            end
        end
        plot([tx_x intersections_mirror1(1) intersections(1) ray_traced_position(1)], [tx_y intersections_mirror1(2) intersections(2), ray_traced_position(2)]);
    end
end

%rue du rouleau
tx_x_mirror1 = tx_x - 40;
tx_y_mirror1 = tx_y;
tx_x_mirror2 = tx_x_mirror1;
tx_y_mirror2 = tx_y - 240;
[amount_of_intersec, intersections] = intersectionCalculator(x, y, ray_traced_position, tx_x_mirror2, tx_y_mirror2, 1);

if(amount_of_intersec == 3 && intersections(1,2) == 40 && intersections(2,3) == 180)
    [amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections(:,3), tx_x_mirror1, tx_y_mirror1, 1);

    [a, b] = intersectionCalculator(x, y, intersections, tx_x, tx_y, 1);
    if(a <= 1 && amount_of_intersec <= 2)
        j = 1;
        while j <= amount_of_intersec
            if(intersections_mirror1(1,j) ~= 0)
                intersections_mirror1(:,j) = [];
                amount_of_intersec = amount_of_intersec - 1;
            else
                j = j + 1;
            end
        end
        plot([tx_x intersections_mirror1(1) intersections(1) ray_traced_position(1)], [tx_y intersections_mirror1(2) intersections(2), ray_traced_position(2)]);
    end
end

%left and rue du peuplier
tx_x_mirror1 = tx_x - 40;
tx_y_mirror1 = tx_y;
tx_x_mirror2 = tx_x_mirror1;
tx_y_mirror2 = tx_y - 460;
[amount_of_intersec, intersections] = intersectionCalculator(x, y, ray_traced_position, tx_x_mirror2, tx_y_mirror2, 1);

if(amount_of_intersec == 3 && intersections(1,2) == 40 && intersections(2,3) == 70)
    [amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections(:,3), tx_x_mirror1, tx_y_mirror1, 1);
    [a, b] = intersectionCalculator(x, y, intersections, tx_x, tx_y, 1);
    if(a <= 1 && amount_of_intersec <= 2)
        j = 1;
        while j <= amount_of_intersec
            if(intersections_mirror1(1,j) ~= 0)
                intersections_mirror1(:,j) = [];
                amount_of_intersec = amount_of_intersec - 1;
            else
                j = j + 1;
            end
        end
        plot([tx_x intersections_mirror1(1) intersections(1) ray_traced_position(1)], [tx_y intersections_mirror1(2) intersections(2), ray_traced_position(2)]);
    end
end


%right wall and left wall
tx_x_mirror1 = tx_x + 40;
tx_y_mirror1 = tx_y;
tx_x_mirror2 = tx_x - 60;
tx_y_mirror2 = tx_y;
[amount_of_intersec, intersections] = intersectionCalculator(x, y, ray_traced_position, tx_x_mirror2, tx_y_mirror2, 1);
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
    [amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections, tx_x_mirror1, tx_y_mirror1, 1);
    j = 1;
    while j <= amount_of_intersec
        if(intersections_mirror1(1,j) ~= 40)
            intersections_mirror1(:,j) = [];
            amount_of_intersec = amount_of_intersec - 1;
        else
            j = j + 1;
        end
    end
    plot([tx_x intersections_mirror1(1) intersections(1) ray_traced_position(1)], [tx_y intersections_mirror1(2) intersections(2), ray_traced_position(2)]);
end


% right and bottom wall
tx_x_mirror1 = tx_x + 40;
tx_y_mirror1 = tx_y;
tx_x_mirror2 = tx_x_mirror1;
tx_y_mirror2 = -tx_y_mirror1;
[amount_of_intersec, intersections] = intersectionCalculator(x, y, ray_traced_position, tx_x_mirror2, tx_y_mirror2, 1);

if(amount_of_intersec == 1 && intersections(2,1) == 0)
    [amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections, tx_x_mirror1, tx_y_mirror1, 1);
    j = 1;
    while j <= amount_of_intersec
        if(intersections_mirror1(1,j) ~= 40)
            intersections_mirror1(:,j) = [];
            amount_of_intersec = amount_of_intersec - 1;
        else
            j = j + 1;
        end
    end
    [a, b] = intersectionCalculator(x, y, intersections, tx_x, tx_y, 1);
    if(a <= 1 && isempty(intersections_mirror1(intersections_mirror1 == 10)))
        plot([tx_x intersections_mirror1(1) intersections(1) ray_traced_position(1)], [tx_y intersections_mirror1(2) intersections(2), ray_traced_position(2)]);  
    end
end

%bottom wall to left wall
tx_x_mirror1 = tx_x;
tx_y_mirror1 = - tx_y;
tx_x_mirror2 = tx_x_mirror1 - 40;
tx_y_mirror2 = tx_y_mirror1;
[amount_of_intersec, intersections] = intersectionCalculator(x, y, ray_traced_position, tx_x_mirror2, tx_y_mirror2, 1);

if(amount_of_intersec == 3 && intersections(2,1) == 0 && intersections(1, 2) == 0)
    [amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections, tx_x_mirror1, tx_y_mirror1, 1);
    j = 1;
    while j <= amount_of_intersec
        if(intersections_mirror1(2,j) ~= 0)
            intersections_mirror1(:,j) = [];
            amount_of_intersec = amount_of_intersec - 1;
        else
            j = j + 1;
        end
    end
    plot([tx_x intersections_mirror1(1) intersections(1) ray_traced_position(1)], [tx_y intersections_mirror1(2) intersections(2), ray_traced_position(2)]);  
end

%bottom wall right wall
tx_x_mirror1 = tx_x;
tx_y_mirror1 = - tx_y;
tx_x_mirror2 = tx_x_mirror1 + 40;
tx_y_mirror2 = tx_y_mirror1;
[amount_of_intersec, intersections] = intersectionCalculator(x, y, ray_traced_position, tx_x_mirror2, tx_y_mirror2, 1);

%exclude intersections with vertical right streets
j = 1;
while j <= amount_of_intersec
    if(intersections(1,j) > 40 || intersections(2, j) == 0)
        intersections(:,j) = [];
        amount_of_intersec = amount_of_intersec - 1;
    else
        j = j + 1;
    end
end

if(amount_of_intersec == 1)
    [amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections(:,2), tx_x_mirror1, tx_y_mirror1, 1);
    j = 1;
    while j <= amount_of_intersec
        if(intersections_mirror1(2,j) ~= 0)
            intersections_mirror1(:,j) = [];
            amount_of_intersec = amount_of_intersec - 1;
        else
            j = j + 1;
        end
    end
    plot([tx_x intersections_mirror1(1) intersections(1) ray_traced_position(1)], [tx_y intersections_mirror1(2) intersections(2), ray_traced_position(2)]);  
end


%bottom wall and street

% bottom wall and top of sint katelijneplijn left and right
tx_x_mirror1 = tx_x;
tx_y_mirror1 = - tx_y;
tx_x_mirror2 = tx_x_mirror1;
tx_y_mirror2 = - tx_y_mirror1 + 20;
[amount_of_intersec, intersections] = intersectionCalculator(x, y, ray_traced_position, tx_x_mirror2, tx_y_mirror2, 1);

%exclude intersections with vertical right streets
j = 1;
while j <= amount_of_intersec
    if(intersections(2, j) > 10)
        intersections(:,j) = [];
        amount_of_intersec = amount_of_intersec - 1;
    else
        j = j + 1;
    end
end

if(amount_of_intersec == 1)
    %[amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections(:,2), tx_x_mirror1, tx_y_mirror1, 1);
    [amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections, tx_x_mirror1, tx_y_mirror1, 1);
    [amount_of_intersec, intersections_double_refl] = intersectionCalculator(x, y, intersections_mirror1, tx_x, tx_y, 1);
    if(amount_of_intersec <= 1)
        j = 1;
        while j <= amount_of_intersec
            if(intersections_mirror1(2,j) ~= 0)
                intersections_mirror1(:,j) = [];
                amount_of_intersec = amount_of_intersec - 1;
            else
                j = j + 1;
            end
        end
        plot([tx_x intersections_mirror1(1) intersections(1) ray_traced_position(1)], [tx_y intersections_mirror1(2) intersections(2), ray_traced_position(2)]);  
    end
end


% bottom wall and top of rue du peuplier
tx_x_mirror1 = tx_x;
tx_y_mirror1 = - tx_y;
tx_x_mirror2 = tx_x_mirror1;
tx_y_mirror2 = - tx_y_mirror1 + 160;
[amount_of_intersec, intersections] = intersectionCalculator(x, y, ray_traced_position, tx_x_mirror2, tx_y_mirror2, 1);

%exclude intersections with vertical right streets
j = 1;
while j <= amount_of_intersec
    if(intersections(2, j) > 80)
        intersections(:,j) = [];
        amount_of_intersec = amount_of_intersec - 1;
    else
        j = j + 1;
    end
end

if(amount_of_intersec == 1)
    %[amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections(:,2), tx_x_mirror1, tx_y_mirror1, 1);
    [amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections, tx_x_mirror1, tx_y_mirror1, 1);
    %only keep intersections with street and right wall
    j = 1;
    while j <= amount_of_intersec
        if((intersections_mirror1(2, j) > 60 || intersections_mirror1(2, j) < 0) && intersections_mirror1(1, j) > 0 )
            intersections_mirror1(:,j) = [];
            amount_of_intersec = amount_of_intersec - 1;
        else
            j = j + 1;
        end
    end
    if(amount_of_intersec == 1)
        plot([tx_x intersections_mirror1(1) intersections(1) ray_traced_position(1)], [tx_y intersections_mirror1(2) intersections(2), ray_traced_position(2)]);  
    end
end

% bottom wall and top of rue du rouleau
% bottom wall and top of rue du grand hospice
% both have no reflections from bottom wall since for one
% reflection the bottom wall does not reach these streets


%double reflection inside streets (excluding bottom streets
%since they are included in bottom wall + streets double reflections)

%rue du peuplier
tx_x_mirror1 = tx_x;
tx_y_mirror1 = - tx_y + 140;
tx_x_mirror2 = tx_x_mirror1;
tx_y_mirror2 = tx_y + 20;
[amount_of_intersec, intersections] = intersectionCalculator(x, y, ray_traced_position, tx_x_mirror2, tx_y_mirror2, 1);

%exclude intersections with vertical right streets
j = 1;
while j <= amount_of_intersec
    if(intersections(2, j) > 80)
        intersections(:,j) = [];
        amount_of_intersec = amount_of_intersec - 1;
    else
        j = j + 1;
    end
end

if(amount_of_intersec == 1)
    %[amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections(:,2), tx_x_mirror1, tx_y_mirror1, 1);
    r_refl = sqrt((all_rx_xy(1,i) - tx_x_mirror2)^2 + (all_rx_xy(2,i) - tx_y_mirror2)^2);
    [amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections, tx_x_mirror1, tx_y_mirror1, 1);
    %only keep intersections with street and right wall
    j = 1;
    while j <= amount_of_intersec
        if(intersections_mirror1(2, j) ~= 70)
            intersections_mirror1(:,j) = [];
            amount_of_intersec = amount_of_intersec - 1;
        else
            j = j + 1;
        end
    end
    if(amount_of_intersec == 1)
        [amount_of_intersec, intersections_LOS] = intersectionCalculator(x, y, intersections_mirror1, tx_x, tx_y, 1);
        if(amount_of_intersec <= 1)
                    plot([tx_x intersections_mirror1(1) intersections(1) ray_traced_position(1)], [tx_y intersections_mirror1(2) intersections(2), ray_traced_position(2)]); 
        end
     end
end

%rue du rouleau
tx_x_mirror1 = tx_x;
tx_y_mirror1 = tx_y - 240;
tx_x_mirror2 = tx_x_mirror1;
tx_y_mirror2 = tx_y + 20;
[amount_of_intersec, intersections] = intersectionCalculator(x, y, ray_traced_position, tx_x_mirror2, tx_y_mirror2, 1);

%exclude intersections with vertical right streets
j = 1;
while j <= amount_of_intersec
    if(intersections(2, j) > 190)
        intersections(:,j) = [];
        amount_of_intersec = amount_of_intersec - 1;
    else
        j = j + 1;
    end
end

if(amount_of_intersec == 1)
    %[amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections(:,2), tx_x_mirror1, tx_y_mirror1, 1);
    r_refl = sqrt((all_rx_xy(1,i) - tx_x_mirror2)^2 + (all_rx_xy(2,i) - tx_y_mirror2)^2);
    [amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections, tx_x_mirror1, tx_y_mirror1, 1);
    %only keep intersections with street and right wall
    j = 1;
    while j <= amount_of_intersec
        if(intersections_mirror1(2, j) ~= 180)
            intersections_mirror1(:,j) = [];
            amount_of_intersec = amount_of_intersec - 1;
        else
            j = j + 1;
        end
    end
    if(amount_of_intersec == 1)
        [amount_of_intersec, intersections_LOS] = intersectionCalculator(x, y, intersections_mirror1, tx_x, tx_y, 1);
        if(amount_of_intersec <= 1)
           plot([tx_x intersections_mirror1(1) intersections(1) ray_traced_position(1)], [tx_y intersections_mirror1(2) intersections(2), ray_traced_position(2)]); 
        end
     end
end
        
%rue du grand hospice
tx_x_mirror1 = tx_x;
tx_y_mirror1 = tx_y - 60;
tx_x_mirror2 = tx_x_mirror1;
tx_y_mirror2 = tx_y + 20;
[amount_of_intersec, intersections] = intersectionCalculator(x, y, ray_traced_position, tx_x_mirror2, tx_y_mirror2, 1);

%exclude intersections with vertical right streets
j = 1;
while j <= amount_of_intersec
    if(intersections(2, j) > 280)
        intersections(:,j) = [];
        amount_of_intersec = amount_of_intersec - 1;
    else
        j = j + 1;
    end
end

if(amount_of_intersec == 1)
    %[amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections(:,2), tx_x_mirror1, tx_y_mirror1, 1);
    [amount_of_intersec, intersections_mirror1] = intersectionCalculator(x, y, intersections, tx_x_mirror1, tx_y_mirror1, 1);
    %only keep intersections with street and right wall
    j = 1;
    while j <= amount_of_intersec
        if(intersections_mirror1(2, j) ~= 270)
            intersections_mirror1(:,j) = [];
            amount_of_intersec = amount_of_intersec - 1;
        else
            j = j + 1;
        end
    end
    if(amount_of_intersec == 1)
        [amount_of_intersec, intersections_LOS] = intersectionCalculator(x, y, intersections_mirror1, tx_x, tx_y, 1);
        if(amount_of_intersec <= 1)
            plot([tx_x intersections_mirror1(1) intersections(1) ray_traced_position(1)], [tx_y intersections_mirror1(2) intersections(2), ray_traced_position(2)]); 
        end
     end
end


    
     


