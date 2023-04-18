function [] = drawRayTracingLines(x, y, all_rx_xy, tx_x, tx_y, no_LOS, ray_traced_position)
%check if there is a LOS
if(isempty(find(no_LOS([1 2],:) == ray_traced_position))) %true is there is LOS
plot([tx_x ray_traced_position(1)], [tx_y ray_traced_position(2)], 'blue');

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
    plot([tx_x diff_point(1) ray_traced_position(1)], [tx_y diff_point(2) ray_traced_position(2)], 'blue');

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
    plot([tx_x intersections(1) ray_traced_position(1)], [tx_y intersections(2), ray_traced_position(2)], 'blue');
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
    plot([tx_x intersections(1) ray_traced_position(1)], [tx_y intersections(2), ray_traced_position(2)], 'blue');
end

% mirror BS to bottom
tx_y_mirror = -tx_y;
tx_x_mirror = tx_x;
[amount_of_intersec, intersections] = intersectionCalculator(x, y, ray_traced_position, tx_x_mirror, tx_y_mirror, 1);

if(amount_of_intersec == 1)
    [a, b] = intersectionCalculator(x, y, intersections, tx_x, tx_y, 1);
    if(a <= 1)
        plot([tx_x intersections(1) ray_traced_position(1)], [tx_y intersections(2), ray_traced_position(2)], 'blue');
    end
end

% mirror BS to Rue du grand hospice
tx_x_mirror = tx_x;
tx_y_mirror = tx_y - 60;
[amount_of_intersec, intersections] = intersectionCalculator(x, y, ray_traced_position, tx_x_mirror, tx_y_mirror, 1);
if(amount_of_intersec == 2 && intersections(1,1) == 40 && intersections(2,2) == 270)
    [a, b] = intersectionCalculator(x, y, intersections(:,2), tx_x, tx_y, 1);
    if(a <= 1)
        plot([tx_x intersections(1, 2) ray_traced_position(1)], [tx_y intersections(2, 2), ray_traced_position(2)], 'blue');
    end
end

% mirror BS to Rue du rouleau
tx_x_mirror = tx_x;
tx_y_mirror = tx_y - 240;
[amount_of_intersec, intersections] = intersectionCalculator(x, y, ray_traced_position, tx_x_mirror, tx_y_mirror, 1);
if(amount_of_intersec == 2 && intersections(1,1) == 40 && intersections(2,2) == 180)
    [a, b] = intersectionCalculator(x, y, intersections(:,2), tx_x, tx_y, 1);
    if(a <= 1)
        plot([tx_x intersections(1, 2) ray_traced_position(1)], [tx_y intersections(2, 2), ray_traced_position(2)], 'blue');
    end
end


% mirror BS to Rue du Peuplier
tx_x_mirror = tx_x;
tx_y_mirror = tx_y - 460;
[amount_of_intersec, intersections] = intersectionCalculator(x, y, ray_traced_position, tx_x_mirror, tx_y_mirror, 1);
if(amount_of_intersec == 3 && intersections(2,1) == 0 && intersections(1,2) == 40 && intersections(2,3) == 70) %also intersects with bottom wall
    [a, b] = intersectionCalculator(x, y, intersections(:,3), tx_x, tx_y, 1);
    if(a <= 1)
        plot([tx_x intersections(1, 2) ray_traced_position(1)], [tx_y intersections(2, 2), ray_traced_position(2)], 'blue');
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
    plot([tx_x intersections_double_refl(1) intersections(1) ray_traced_position(1)], [tx_y intersections_double_refl(2) intersections(2), ray_traced_position(2)], 'blue');
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
        plot([tx_x intersections_mirror1(1) intersections(1) ray_traced_position(1)], [tx_y intersections_mirror1(2) intersections(2), ray_traced_position(2)], 'blue');
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
        plot([tx_x intersections_mirror1(1) intersections(1) ray_traced_position(1)], [tx_y intersections_mirror1(2) intersections(2), ray_traced_position(2)], 'blue');
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
        plot([tx_x intersections_mirror1(1) intersections(1) ray_traced_position(1)], [tx_y intersections_mirror1(2) intersections(2), ray_traced_position(2)], 'blue');
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
        plot([tx_x intersections_mirror1(1) intersections(1) ray_traced_position(1)], [tx_y intersections_mirror1(2) intersections(2), ray_traced_position(2)], 'blue');
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
    plot([tx_x intersections_mirror1(1) intersections(1) ray_traced_position(1)], [tx_y intersections_mirror1(2) intersections(2), ray_traced_position(2)], 'blue');
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
        plot([tx_x intersections_mirror1(1) intersections(1) ray_traced_position(1)], [tx_y intersections_mirror1(2) intersections(2), ray_traced_position(2)], 'blue');  
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
    plot([tx_x intersections_mirror1(1) intersections(1) ray_traced_position(1)], [tx_y intersections_mirror1(2) intersections(2), ray_traced_position(2)], 'blue');  
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
    plot([tx_x intersections_mirror1(1) intersections(1) ray_traced_position(1)], [tx_y intersections_mirror1(2) intersections(2), ray_traced_position(2)], 'blue');  
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
        plot([tx_x intersections_mirror1(1) intersections(1) ray_traced_position(1)], [tx_y intersections_mirror1(2) intersections(2), ray_traced_position(2)], 'blue');  
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
        plot([tx_x intersections_mirror1(1) intersections(1) ray_traced_position(1)], [tx_y intersections_mirror1(2) intersections(2), ray_traced_position(2)], 'blue');  
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
                    plot([tx_x intersections_mirror1(1) intersections(1) ray_traced_position(1)], [tx_y intersections_mirror1(2) intersections(2), ray_traced_position(2)], 'blue'); 
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
           plot([tx_x intersections_mirror1(1) intersections(1) ray_traced_position(1)], [tx_y intersections_mirror1(2) intersections(2), ray_traced_position(2)], 'blue'); 
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
            plot([tx_x intersections_mirror1(1) intersections(1) ray_traced_position(1)], [tx_y intersections_mirror1(2) intersections(2), ray_traced_position(2)], 'blue'); 
        end
     end
end
end