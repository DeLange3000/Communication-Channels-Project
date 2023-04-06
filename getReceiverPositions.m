function [all_rx_xy] = getReceiverPositions(d_basestation, tx_x, tx_y)
%creates array with all possible receiver positions
width = 40;
height = 320;
width_streets = 20;
height_streets = 10;
all_rx_xy = [];
% fill in vismarkt center
for i = 1:width
    for j = 1:height
        all_rx_xy = [all_rx_xy [i-0.5; j-0.5]];
    end
end

%add side streets
offsets_xy = [-20, 0;  40, 0; 40, 70; 40, 180; 40, 270];
for a = 1:5
    for i = 1:width_streets
        for j = 1:height_streets
            all_rx_xy= [all_rx_xy [i-0.5 + offsets_xy(a,1); j - 0.5 + offsets_xy(a,2)]];
        end
    end
end



%removing points closer then minimal distance to base station
i = 1;
while i < length(all_rx_xy(1,:))
    if(sqrt((all_rx_xy(1,i) - tx_x)^2 + (all_rx_xy(2,i)-tx_y)^2) < d_basestation)
        all_rx_xy(:,i) = [];
    else
        i = i + 1;
    end
    
end

end