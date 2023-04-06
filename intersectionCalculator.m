function [return_value, intersections] = intersectionCalculator(walls_x, walls_y, rx_xy, tx_x, tx_y, a)

A = [-(rx_xy(2) - tx_y)/(rx_xy(1) - tx_x) 1; 0 0];
B = [tx_y + A(1,1)*tx_x; 0];

return_value = 0;
amount_of_intersections = 0;
intersections = [];
for j = 1:length(walls_x(1,:))
    %turn walls into lines
    A(2,:) = [-(walls_y(2,j) - walls_y(1,j))/(walls_x(2,j) - walls_x(1,j)) 1];
    B(2) = [walls_y(1,j) + A(2,1)*walls_x(1,j)];
    %calculate intersection of LOS line with walls
    intersec_y = (B(2) - (A(2,1)*B(1))/A(1,1))/(A(2,2) - (A(2,1)*A(1,2))/A(1,1));
    intersec_x = (B(1) - A(1,2)*intersec_y)/A(1,1);

    %check if intersection point falls within wall and LOS line limits
    %(they have finite length)
    if(intersec_x > max([tx_x, rx_xy(1)]) || intersec_x < min([tx_x, rx_xy(1)]))
        intersec_x = NaN;
    elseif(intersec_x > max([walls_x(1,j), walls_x(2,j)]) || intersec_x < min([walls_x(1,j), walls_x(2,j)]))
        intersec_x = NaN;
    end
    if (intersec_y > max([tx_y, rx_xy(2)]) || intersec_y < min([tx_y, rx_xy(2)]))
        intersec_y = NaN;
    elseif(intersec_y > max([walls_y(1,j), walls_y(2,j)]) || intersec_y < min([walls_y(1,j), walls_y(2,j)]))
        intersec_y = NaN;
    end

    % add non LOS points to an array
    if (not(isnan(intersec_x) || isnan(intersec_y)))
        if(a == 0)
        return_value = 1;
        return
        end

        if(a == 1)
            amount_of_intersections = amount_of_intersections + 1;
            intersections = [intersections [intersec_x; intersec_y]];
        end
       
    end

    return_value = amount_of_intersections;


    end
end