%% Mackey-Glass Chaotic Time Series Prediction problem

clc;                % clear command window
close all;          % close all figures
clearvars;          % clear all workspace variables
%% 
% Generate and load Mackey-Glass chaotic time series

mackeyglass; % generate Mackey-Glass chaotic time series data
load mgchaotic.dat; % load mackey-glass time series data

% Plot the time series for first 1000 points
plot(mgchaotic(1:1000,1),mgchaotic(1:1000,2),'LineWidth',1);
xlabel('t');ylabel('x(t)');
title('Mackey-Glass Chaotic Time Series for {\tau} = 30');

% Fix some parameters as mentioned in paper
tau = 30;           % delay parameter
m = 9;              % no of points taken for prediction
l = 1;              % we will be predicting next point from m points
M = size(mgchaotic,1);  % size of dataset
%%
% Now let us form (M-m) input-output pairs as mentioned in the paper
input_data = zeros(M-m,m);
output_data = zeros(M-m,1);

for i = 1:(M-m)
    input_data(i,:) = mgchaotic(i:i+m-1,2)';
    output_data(i,:) = mgchaotic(i+m,2);
end
%%
% Now let us divide the data into train and test set
% First 700 data pairs will constitute our training data
% Next 300 data pairs will constitute our testing data

train_data = [input_data(1:700,:) output_data(1:700)];
test_data = [input_data(701:1000,:) output_data(701:1000)];
%%
N = [3 7 14];
x_min = min(mgchaotic(1:end,2));
x_max = max(mgchaotic(1:end,2));
disp([x_min x_max]);
% Since min is > 0.2 and max is < 1.4, lets take the domain
% interval [0.2, 1.4]
x_dom = 0.2:0.01:1.4;

% We define a function to complete the 5 steps of training
% and testing as mentioned in the paper
%%
%function pred_output = applyNumericFuzzy(train_data, test_data, x_dom, N)
for j=1:size(N,2)
    % Step 1
    R_count = 2*N(j) + 1; % no of regions into which we divide domain interval of x
    fuzzy_regions = timeseriesregions(N(j), x_dom);
    figure;
    hax = axes;
    for i=1:R_count
        y_loc = find(fuzzy_regions(:,i)==1);
        y_t(1,i) = x_dom(1,y_loc);
        plot(hax,x_dom, fuzzy_regions(:,i));
        title(sprintf('Membership Function for Chaotic Time Series Prediction N=%d',N(j)));
        xlabel('x(t)');
        ylabel('\mu(x)');
        hold on;
    end
    hold off;
    saveas(gcf,sprintf("chaotic_membership_%d.png",N(j)));

    % Step 2
    n = size(train_data,1);
    rules = zeros(n,m+1);
    degrees = zeros(n,m+1);
    % regional_degree = zeros(n,m+1,R_count);
    for k=1:n
        for l=1:(m+1)
            regional_degree = zeros(R_count,1);
            train_data(k,l) = round(train_data(k,l)*100)/100;
            for i = 1:R_count
                temp_matrix = abs(x_dom - repmat(train_data(k,l),size(x_dom)));
                dom_idx = find(temp_matrix < 0.001);
                if (dom_idx)
                    regional_degree(i,:) = fuzzy_regions(dom_idx,i);
                end
            end
            [degrees(k,l), rules(k,l)] = max(regional_degree);
         end
    end
%%
    % Step 3
    rules_degree = prod(degrees, 2);
    [temp,ia,ic] = unique(rules(:,1:m), "rows", "sorted");
    ctr = 1;
    for i=1:size(temp,1)
        dup_idx = find(ic==i);
        [~,ib] = max(rules_degree(dup_idx));
        final_rules(ctr,:) = rules(dup_idx(ib),:);
        ctr = ctr + 1;
    end
    % Step 4
    fuzzy_rule_base = final_rules;
    
%%
    %Step 5
    Degree_test = 0;
    in_mf_prod = 0;
    y_bar = zeros(1,1);
    test_output = 0;
    for k=1:size(test_data,1)
        sample = test_data(k,:);
        for l=1:size(fuzzy_rule_base,1)
            for r = 1:size(test_data,2)-1
                regional_degree = zeros(R_count,1);
                testdata = round(sample(r)*100)/100;
                for i = 1:R_count
                    temp_matrix = abs(x_dom - repmat(testdata,size(x_dom)));
                    dom_idx = find(temp_matrix < 0.001);
                    if (dom_idx)
                        regional_degree(i,:) = fuzzy_regions(dom_idx,i);
                    end
                end
                Degree_test(1,r) = regional_degree(fuzzy_rule_base(l,r));
            end
            in_mf_prod(l,:) = prod(Degree_test);
            y_bar(l,:) = y_t(fuzzy_rule_base(l,10));
        end
        test_output(k,1) = sum(in_mf_prod.*y_bar)/sum(in_mf_prod) ;
    end
    fprintf("mean square error=%f\n",norm(test_output-test_data(:,end)));
    test_outputs(j,:) = test_output;
end
%%
figure;
actual_output= test_data(:,end);
plot(1:300,actual_output,'LineWidth',2.1);hold on;
plot(1:300,test_outputs(1,:),'LineWidth',1.7);hold on;
plot(1:300,test_outputs(2,:),'LineWidth',1.5);hold on;
plot(1:300,test_outputs(3,:),'LineWidth',1.5);
title('Time Series Prediction with 700 train data');
xlabel('t','FontSize',12);
ylabel('x(t)','FontSize',12);
legend('Actual Value','Prediction at N=3','Prediction at N=7','Prediction at N=14');
legend('Location','southwest');
set(gcf,'position',[50 50 1000 600]);
saveas(gcf,"timeseries_700.png");