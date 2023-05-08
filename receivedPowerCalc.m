function [received_power, rice_factor, delay_spread, impulse_response] = receivedPowerCalc(Ra, c, permitivity, x, y, all_rx_xy, tx_x, tx_y, effective_h, beta, h_ue, h_bs, EIRP, Fc, no_LOS, enable_LOS, enable_ground_reflections, ...
enable_diffraction, enable_one_reflection, enable_left_wall_reflection, enable_right_wall_reflection, enable_bottom_wall_reflection, eanble_rue_du_grand_hospice, enable_rue_du_rouleau, ...
enable_rue_du_peuplier, enable_double_reflections, enable_left_and_right_wall, enable_left_and_bottom_wall, enable_left_and_streets, enable_right_and_left_wall, enable_right_and_bottom_wall, ... 
enable_bottom_and_left_wall, enable_bottom_and_right_wall, enable_bottom_and_streets, enable_streets_and_streets)

%variables
received_power = zeros(size(all_rx_xy(1,:)));
rice_factor = zeros(size(received_power));
delay_spread = zeros(size(received_power));
impulse_response = zeros(length(received_power), 2, 23);

for i = 1:length(all_rx_xy(1,:))
    Voc = 0;
    a0 = 0;
    an = 0;
    delays = [];
    %check if there is a LOS
    if(isempty(find(no_LOS(3,:) == i))) %true is there is LOS
        if(enable_LOS == 1)
            r_los = sqrt((all_rx_xy(1,i) - tx_x)^2 + (all_rx_xy(2,i) - tx_y)^2); %since same height
            E_los = sqrt(60*EIRP)/r_los*exp(-1i*2*pi*Fc/c*r_los);
            Voc = Voc + E_los*effective_h;
            a0 = abs(E_los*effective_h)^2;
            delays = [delays r_los/c];
            impulse_response(i, :, 1) = [E_los*effective_h, r_los/c];
        end

        %if LOS then also ground reflection
        if(enable_ground_reflections == 1)
            r_los = sqrt((all_rx_xy(1,i) - tx_x)^2 + (all_rx_xy(2,i) - tx_y)^2); %since same height
            r_gr = sqrt(r_los^2 + (2*h_bs)^2);
            theta = acos(h_bs/(r_gr/2));
            effective_h_gr = -(c/Fc)/pi*cos(pi/2*cos(theta))/sin(theta)^2;
            lamba_parallel = (cos(theta) - sqrt(1/permitivity)*sqrt(1-1/permitivity*sin(theta)^2))/(cos(theta) + sqrt(1/permitivity)*sqrt(1-1/permitivity*sin(theta)^2));
            E_gr = lamba_parallel*sqrt(60*EIRP)/r_gr*exp(-1i*2*pi*Fc/c*r_gr);
            Voc = Voc + E_gr*effective_h_gr*sin(theta);
            an = an + abs(E_gr*effective_h_gr*sin(theta))^2;
            delays = [delays r_gr/c];
            impulse_response(i, :, 2) = [E_gr*effective_h_gr*sin(theta), r_gr/c];
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
            an = an + abs(E_diff*effective_h)^2;
            delays = [delays r_diff/c];
            impulse_response(i, :, 23) = [E_diff*effective_h, r_diff/c];
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
                an = an + abs(E_refl*effective_h)^2;
                delays = [delays r_refl/c];
                impulse_response(i, :, 3) = [E_refl*effective_h, r_refl/c];
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
                an = an + abs(E_refl*effective_h)^2;
                delays = [delays r_refl/c];
                impulse_response(i, :, 4) = [E_refl*effective_h, r_refl/c];
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
                    an = an + abs(E_refl*effective_h)^2;
                    delays = [delays r_refl/c];
                    impulse_response(i, :, 5) = [E_refl*effective_h, r_refl/c];
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
                    an = an + abs(E_refl*effective_h)^2;
                    delays = [delays r_refl/c];
                    impulse_response(i, :, 6) = [E_refl*effective_h, r_refl/c];
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
                [a, ~] = intersectionCalculator(x, y, intersections(:,2), tx_x, tx_y, 1);
                if(a <= 1)
                    r_refl = sqrt((all_rx_xy(1,i) - tx_x_mirror)^2 + (all_rx_xy(2,i) - tx_y_mirror)^2);
                    theta = atan(abs(abs(tx_y_mirror) - abs(intersections(2)))/abs(tx_x_mirror - abs(intersections(1))));
                    lamba_ortho = (cos(theta) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta)^2))/(cos(theta) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta)^2));
                    E_refl = lamba_ortho*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                    Voc = Voc + E_refl*effective_h;
                    an = an + abs(E_refl*effective_h)^2;
                    delays = [delays r_refl/c];
                    impulse_response(i, :, 7) = [E_refl*effective_h, r_refl/c];
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
                    an = an + abs(E_refl*effective_h)^2;
                    delays = [delays r_refl/c];
                    impulse_response(i, :, 8) = [E_refl*effective_h, r_refl/c];
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
                    theta2 = atan(abs(tx_y_mirror2 - all_rx_xy(2,i))/abs(tx_x_mirror2 - all_rx_xy(1,i)));
                    theta1 = theta2;
                    lamba_ortho1 = (cos(theta1) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2))/(cos(theta1) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2));
                    lamba_ortho2 = (cos(theta2) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2))/(cos(theta2) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2));
                    E_refl = lamba_ortho1*lamba_ortho2*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                    Voc = Voc + E_refl*effective_h;
                    an = an + abs(E_refl*effective_h)^2;
                    delays = [delays r_refl/c];
                    impulse_response(i, :, 9) = [E_refl*effective_h, r_refl/c];
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
                        theta2 = atan(abs(tx_x_mirror2 - all_rx_xy(1,i))/abs(tx_y_mirror2 - all_rx_xy(2,i)));
                        theta1 = pi/2 - theta2;
                        lamba_ortho1 = (cos(theta1) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2))/(cos(theta1) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2));
                        lamba_ortho2 = (cos(theta2) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2))/(cos(theta2) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2));
                        E_refl = lamba_ortho1*lamba_ortho2*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                        Voc = Voc + E_refl*effective_h;
                        an = an + abs(E_refl*effective_h)^2;
                        delays = [delays r_refl/c];
                        impulse_response(i, :, 10) = [E_refl*effective_h, r_refl/c];
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
                        theta2 = atan(abs(tx_x_mirror2 - all_rx_xy(1,i))/abs(tx_y_mirror2 - all_rx_xy(2,i)));
                        theta1 = pi/2 - theta2;
                        lamba_ortho1 = (cos(theta1) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2))/(cos(theta1) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2));
                        lamba_ortho2 = (cos(theta2) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2))/(cos(theta2) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2));
                        E_refl = lamba_ortho1*lamba_ortho2*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                        Voc = Voc + E_refl*effective_h;
                        an = an + abs(E_refl*effective_h)^2; 
                        delays = [delays r_refl/c];
                        impulse_response(i, :, 11) = [E_refl*effective_h, r_refl/c];
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
                        theta2 = atan(abs(tx_x_mirror2 - all_rx_xy(1,i))/abs(tx_y_mirror2 - all_rx_xy(2,i)));
                        theta1 = pi/2 - theta2;
                        lamba_ortho1 = (cos(theta1) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2))/(cos(theta1) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2));
                        lamba_ortho2 = (cos(theta2) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2))/(cos(theta2) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2));
                        E_refl = lamba_ortho1*lamba_ortho2*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                        Voc = Voc + E_refl*effective_h;
                        an = an + abs(E_refl*effective_h)^2;
                        delays = [delays r_refl/c];
                        impulse_response(i, :, 12) = [E_refl*effective_h, r_refl/c];
                    end
                end

                %rue du peuplier
                tx_x_mirror1 = tx_x - 40;
                tx_y_mirror1 = tx_y;
                tx_x_mirror2 = tx_x_mirror1;
                tx_y_mirror2 = tx_y - 460;
                [amount_of_intersec, intersections] = intersectionCalculator(x, y, all_rx_xy(:,i), tx_x_mirror2, tx_y_mirror2, 1);
              
                if(amount_of_intersec == 3 && intersections(1,2) == 40 && intersections(2,3) == 70)
                    [amount_of_intersec, ~] = intersectionCalculator(x, y, intersections(:,3), tx_x_mirror1, tx_y_mirror1, 1);
                    [a, b] = intersectionCalculator(x, y, intersections, tx_x, tx_y, 1);
                    if(a <= 1 && amount_of_intersec <= 2)
                        r_refl = sqrt((all_rx_xy(1,i) - tx_x_mirror2)^2 + (all_rx_xy(2,i) - tx_y_mirror2)^2);
                        theta2 = atan(abs(tx_x_mirror2 - all_rx_xy(1,i))/abs(tx_y_mirror2 - all_rx_xy(2,i)));
                        theta1 = pi/2 - theta2;
                        lamba_ortho1 = (cos(theta1) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2))/(cos(theta1) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2));
                        lamba_ortho2 = (cos(theta2) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2))/(cos(theta2) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2));
                        E_refl = lamba_ortho1*lamba_ortho2*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                        Voc = Voc + E_refl*effective_h;
                        an = an + abs(E_refl*effective_h)^2;
                        delays = [delays r_refl/c];
                        impulse_response(i, :, 13) = [E_refl*effective_h, r_refl/c];
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
                    if(amount_of_intersec >= 1)
                        r_refl = sqrt((all_rx_xy(1,i) - tx_x_mirror2)^2 + (all_rx_xy(2,i) - tx_y_mirror2)^2);
                        theta2 = atan(abs(tx_y_mirror2 - all_rx_xy(2,i))/abs(tx_x_mirror2 - all_rx_xy(1,i)));
                        theta1 = theta2;
                        lamba_ortho1 = (cos(theta1) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2))/(cos(theta1) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2));
                        lamba_ortho2 = (cos(theta2) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2))/(cos(theta2) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2));
                        E_refl = lamba_ortho1*lamba_ortho2*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                        Voc = Voc + E_refl*effective_h;
                        an = an + abs(E_refl*effective_h)^2;
                        delays = [delays r_refl/c];
                        impulse_response(i, :, 14) = [E_refl*effective_h, r_refl/c];
                    end
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
                        theta2 = atan(abs((tx_x_mirror2 - all_rx_xy(1,i)))/(tx_y_mirror2 - all_rx_xy(2,i)));
                        theta1 = pi/2 - theta2;
                        lamba_ortho1 = (cos(theta1) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2))/(cos(theta1) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2));
                        lamba_ortho2 = (cos(theta2) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2))/(cos(theta2) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2));
                        E_refl = lamba_ortho1*lamba_ortho2*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                        Voc = Voc + E_refl*effective_h;
                        an = an + abs(E_refl*effective_h)^2;
                        delays = [delays r_refl/c];
                        impulse_response(i, :, 15) = [E_refl*effective_h, r_refl/c];
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
                    [~, intersections_mirror1] = intersectionCalculator(x, y, intersections, tx_x_mirror1, tx_y_mirror1, 1);
                    r_refl = sqrt((all_rx_xy(1,i) - tx_x_mirror2)^2 + (all_rx_xy(2,i) - tx_y_mirror2)^2);

                    theta2 = atan(abs(tx_y_mirror2 - all_rx_xy(2,i))/abs(tx_x_mirror2 - all_rx_xy(1,i)));
                    theta1 = pi/2 - theta2;
                    lamba_ortho1 = (cos(theta1) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2))/(cos(theta1) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2));
                    lamba_ortho2 = (cos(theta2) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2))/(cos(theta2) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2));
                    E_refl = lamba_ortho1*lamba_ortho2*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                    Voc = Voc + E_refl*effective_h;
                    an = an + abs(E_refl*effective_h)^2;
                    delays = [delays r_refl/c];
                    impulse_response(i, :, 16) = [E_refl*effective_h, r_refl/c];
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

                    theta2 = atan(abs(tx_y_mirror2 - all_rx_xy(2,i))/abs(tx_x_mirror2 - all_rx_xy(1,i)));
                    theta1 = pi/2 - theta2;
                    lamba_ortho1 = (cos(theta1) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2))/(cos(theta1) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2));
                    lamba_ortho2 = (cos(theta2) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2))/(cos(theta2) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2));
                    E_refl = lamba_ortho1*lamba_ortho2*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                    Voc = Voc + E_refl*effective_h;
                    an = an + abs(E_refl*effective_h)^2;
                    delays = [delays r_refl/c];
                    impulse_response(i, :, 17) = [E_refl*effective_h, r_refl/c];
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

                        theta2 = atan(abs(tx_x_mirror2 - all_rx_xy(1,i))/abs(tx_y_mirror2 - all_rx_xy(2,i)));
                        theta1 = theta2;
                        lamba_ortho1 = (cos(theta1) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2))/(cos(theta1) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2));
                        lamba_ortho2 = (cos(theta2) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2))/(cos(theta2) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2));
                        E_refl = lamba_ortho1*lamba_ortho2*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                        Voc = Voc + E_refl*effective_h;
                        an = an + abs(E_refl*effective_h)^2;
                        delays = [delays r_refl/c];
                        impulse_response(i, :, 18) = [E_refl*effective_h, r_refl/c];
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
                        theta2 = atan(abs(tx_x_mirror2 - all_rx_xy(1,i))/abs(tx_y_mirror2 - all_rx_xy(2,i)));
                        theta1 = theta2;
                        lamba_ortho1 = (cos(theta1) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2))/(cos(theta1) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2));
                        lamba_ortho2 = (cos(theta2) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2))/(cos(theta2) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2));
                        E_refl = lamba_ortho1*lamba_ortho2*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                        Voc = Voc + E_refl*effective_h;
                        an = an + abs(E_refl*effective_h)^2;
                        delays = [delays r_refl/c];
                        impulse_response(i, :, 19) = [E_refl*effective_h, r_refl/c];
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
                            theta2 = atan(abs(tx_x_mirror2 - all_rx_xy(1,i))/abs(tx_y_mirror2 - all_rx_xy(2,i)));
                            theta1 = theta2;
                            lamba_ortho1 = (cos(theta1) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2))/(cos(theta1) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2));
                            lamba_ortho2 = (cos(theta2) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2))/(cos(theta2) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2));
                            E_refl = lamba_ortho1*lamba_ortho2*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                            Voc = Voc + E_refl*effective_h;
                            an = an + abs(E_refl*effective_h)^2;
                            delays = [delays r_refl/c];
                            impulse_response(i, :, 20) = [E_refl*effective_h, r_refl/c];
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
                            theta2 = atan(abs(tx_x_mirror2 - all_rx_xy(1,i))/abs(tx_y_mirror2 - all_rx_xy(2,i)));
                            theta1 = theta2;
                            lamba_ortho1 = (cos(theta1) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2))/(cos(theta1) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2));
                            lamba_ortho2 = (cos(theta2) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2))/(cos(theta2) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2));
                            E_refl = lamba_ortho1*lamba_ortho2*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                            Voc = Voc + E_refl*effective_h;
                            an = an + abs(E_refl*effective_h)^2;
                            delays = [delays r_refl/c];
                            impulse_response(i, :, 21) = [E_refl*effective_h, r_refl/c];
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
                            theta2 = atan(abs(tx_x_mirror2 - all_rx_xy(1,i))/abs(tx_y_mirror2 - all_rx_xy(2,i)));
                            theta1 = theta2;
                            lamba_ortho1 = (cos(theta1) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2))/(cos(theta1) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta1)^2));
                            lamba_ortho2 = (cos(theta2) - sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2))/(cos(theta2) + sqrt(permitivity)*sqrt(1-1/permitivity*sin(theta2)^2));
                            E_refl = lamba_ortho1*lamba_ortho2*sqrt(60*EIRP)/r_refl*exp(-1i*2*pi*Fc/c*r_refl);
                            Voc = Voc + E_refl*effective_h;
                            an = an + abs(E_refl*effective_h)^2;
                            delays = [delays r_refl/c];
                            impulse_response(i, :, 22) = [E_refl*effective_h, r_refl/c];
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
    rice_factor(i) = a0/an;
    impulse_response(i,1,:) = impulse_response(i,1,:)./(Ra*sqrt(2*received_power(i)/Ra));

    if(length(delays) == 1)
        delay_spread(i) = 0;
    else
        delay_spread(i) =  max(delays) - min(delays);
    end
end
end