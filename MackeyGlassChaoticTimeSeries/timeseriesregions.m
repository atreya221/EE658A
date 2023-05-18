% We define a function to divide domain interval of x into
% 2N+1 fuzzy regions
function [fuzzy_regions] = timeseriesregions(N, x_dom)
    R_count = 2*N + 1;
    x_min = x_dom(1);
    x_max = x_dom(end);
    vertices = zeros(R_count,3);
    halfbase = (x_max-x_min)/(2*N);
    vertices(1,:) = [x_min-halfbase x_min x_min+halfbase];
    vertices(R_count,:) = [x_max-halfbase x_max x_max+halfbase];
    fuzzy_regions = zeros(size(x_dom,2),R_count);
    fuzzy_regions(:,1) = trapmf(x_dom, [vertices(1,1) vertices(1,1) vertices(1,2) vertices(1,3)]);
    fuzzy_regions(:,R_count) = trapmf(x_dom, [vertices(R_count,1) vertices(R_count,2) vertices(R_count,3) vertices(R_count,3)]);

    for i = 2:(R_count-1)
        vertices(i,:) = [x_min+(i-2)*halfbase x_min+(i-1)*halfbase x_min+i*halfbase];
        fuzzy_regions(:,i) = trimf(x_dom, [vertices(i,1) vertices(i,2) vertices(i,3)]);
        [max_val, max_indx] = max(fuzzy_regions(:,i));
        fuzzy_regions(max_indx,i) = ceil(max_val);
    end
end