%% Example 3 - Truck Backer-Upper Problem (using Numerical Fuzzy approach)

clc;                % clear command window
close all;          % close all figures
clearvars;          % clear all workspace variables
% 
% Load data tables from csv files
% Now, we load first three data pairs of desired trajectories from 14 different 
% initial states.

dataDir = dir("data\*.csv");        % directory containing data tables
concatTable = [];                   % concatenate all data tables and store
tables = cell(1,14);                % store all tables in cell array
for i=1:length(dataDir)             % iterate over all tables
    file = dataDir(1).folder + "\table" + num2str(i) + ".csv";  % i-th data table
    tables{i} = readtable(file);                % read table from csv file
    truncTable = tables{i}(1:3,:);
    concatTable = [concatTable; truncTable];     % concatenate i-th table
end
concatTable = table2array(concatTable);         % convert table to double array
clear dataDir i file;
% 
% Generating fuzzy rules from Numerical Data
% For generating fuzzy rules from the given data, we follow the 5 step procedure 
% proposed in the paper.
%% 
% 
% Step 1 - Divide the Input and Output spaces into Fuzzy regions
% First, we identify the input-output data pairs. In our case the $i^{th}$ input-output 
% pair is $(\:x^{(i)}, \:\phi^{(i)}\:; \:\theta^{(i)}\:)$. Now, we fix the *domain 
% intervals* of the input variables ($x$ and $\phi$) and the output variable ($\theta$). 
% As provided in the paper, we fix the domain intervals as $x \in [0,20]$, $\phi 
% \in [-90\degree, 270\degree]$ and $\theta \in [-40\degree, 40\degree]$. Now 
% we divide the domain interval of each variable into $2N + 1$ regions and assign 
% each region a membership function.
% 
% 
% 
% For the input $x$, let us take $N = 2$. Hence we assign the 5 regions formed 
% (S2, S1, CE, B1, B2) triangular membership functions  as follows (except for 
% the end intervals S2 and B2, which are assigned trapezoidal membership functions).

x_dom = 0:0.1:20;                                   % domain interval of x
x_mf = zeros(201, 5);
x_mf(:,1) = trapmf(x_dom, [0 0 1.5 7]);             % S2
x_mf(:,2) = trimf(x_dom, [4, 7, 10]);               % S1
x_mf(:,3) = trimf(x_dom, [9, 10, 11]);              % CE
x_mf(:,4) = trimf(x_dom, [10, 13, 16]);             % B1
x_mf(:,5) = trapmf(x_dom, [13, 18.5, 20, 20]);      % B2
subplot(311);
plot(x_dom, x_mf, 'Linewidth', 1.5);
xlabel('x');ylabel('\mu(x)');
title('Fuzzy Membership Function for Input x');
legend('S2','S1','CE','B1','B2');
legend('Location','northeastoutside');
xlim([0 20]);
ylim([0 1]);
%% 
% For the input $\phi$, let us take $N = 3$. Hence we assign the 7 regions formed 
% (S3, S2, S1, CE, B1, B2, B3) triangular membership functions (The end regions 
% S3 and B3 extend beyond the domain interval).

phi_dom = -90:0.1:270;                              % domain interval of phi
phi_mf = zeros(3601, 7);
phi_mf(:,1) = trimf(phi_dom, [-115, -65, -15]);     % S3
phi_mf(:,2) = trimf(phi_dom, [-45, 0, 45]);         % S2
phi_mf(:,3) = trimf(phi_dom, [15, 52.5, 90]);       % S1
phi_mf(:,4) = trimf(phi_dom, [80, 90, 100]);        % CE
phi_mf(:,5) = trimf(phi_dom, [90, 127.5, 165]);     % B1
phi_mf(:,6) = trimf(phi_dom, [135, 180, 225]);      % B2
phi_mf(:,7) = trimf(phi_dom, [195, 245, 295]);      % B3
subplot(312);
plot(phi_dom, phi_mf, 'Linewidth', 1.5);
xlabel('\phi');ylabel('\mu(\phi)');
title('Fuzzy Membership Function for Input \phi');
legend('S3', 'S2','S1','CE','B1','B2', 'B3');
legend('Location','northeastoutside');
xlim([-90 270]);
ylim([0 1]);
%% 
% For the output $\theta$, let us take $N = 3$. Hence we assign the 7 regions 
% formed (S3, S2, S1, CE, B1, B2, B3) triangular membership functions (The end 
% regions S3 and B3 extend beyond the domain interval).

theta_dom = -40:0.1:40;                             % domain interval of theta
theta_mf = zeros(801, 7);
theta_mf(:,1) = trimf(theta_dom, [-40, -40, -20]);  % S3
theta_mf(:,2) = trimf(theta_dom, [-33, -20, -7]);   % S2
theta_mf(:,3) = trimf(theta_dom, [-14, -7, 0]);     % S1
theta_mf(:,4) = trimf(theta_dom, [-4, 0, 4]);       % CE
theta_mf(:,5) = trimf(theta_dom, [0, 7, 14]);       % B1
theta_mf(:,6) = trimf(theta_dom, [7, 20, 33]);      % B2
theta_mf(:,7) = trimf(theta_dom, [20, 40, 40]);     % B3
subplot(313);
plot(theta_dom, theta_mf, 'Linewidth', 1.5);
xlabel('\theta');ylabel('\mu(\theta)');
title('Fuzzy Membership Function for Output \theta');
legend('S3', 'S2','S1','CE','B1','B2', 'B3');
legend('Location','northeastoutside');
xlim([-40 40]);
ylim([0 1]);
% 
% Step 2 - Generate Fuzzy rules from given Data pairs
% To determine the fuzzy rules, we first determine the degrees of given $x^{(i)}$, 
% $\phi^{(i)}$ and $\theta^{(i)}$ in different regions that we had previously 
% defined. 
% 
% 
% 
% 1. Degree of input $x$ in different regions

x_deg_S2 = trapmf(concatTable(:,1), [0 0 1.5 7]);           % deg of x in S2
x_deg_S1 = trimf(concatTable(:,1), [4, 7, 10]);             % deg of x in S1
x_deg_CE = trimf(concatTable(:,1), [9, 10, 11]);            % deg of x in CE
x_deg_B1 = trimf(concatTable(:,1), [10, 13, 16]);           % deg of x in B1
x_deg_B2 = trapmf(concatTable(:,1), [13, 18.5, 20, 20]);    % deg of x in B2
%% 
% Now, we define the degree of $x$ as the maximum of the degrees of $x$ in different 
% regions. And then, we assign $x$ to the region with maximum degree.

[x_deg, x_reg] = max([x_deg_S2 x_deg_S1 x_deg_CE x_deg_B1 x_deg_B2], [], 2);
%% 
% 
% 
% 2. Degree of input $\phi$ in different regions

phi_deg_S3 = trimf(concatTable(:,2), [-115, -65, -15]);     % deg of phi in S3
phi_deg_S2 = trimf(concatTable(:,2), [-45, 0, 45]);         % deg of phi in S2
phi_deg_S1 = trimf(concatTable(:,2), [15, 52.5, 90]);       % deg of phi in S1
phi_deg_CE = trimf(concatTable(:,2), [80, 90, 100]);        % deg of phi in CE
phi_deg_B1 = trimf(concatTable(:,2), [90, 127.5, 165]);     % deg of phi in B1
phi_deg_B2 = trimf(concatTable(:,2), [135, 180, 225]);      % deg of phi in B2
phi_deg_B3 = trimf(concatTable(:,2), [195, 245, 295]);      % deg of phi in B3
%% 
% Similar to $x$, we define the degree of $\phi$ and assign $\phi$ to the region 
% with maximum degree.

[phi_deg, phi_reg] = max([phi_deg_S3, phi_deg_S2, phi_deg_S1, ...
    phi_deg_CE, phi_deg_B1, phi_deg_B2, phi_deg_B3], [], 2);
% 
% 3. Degree of output $\theta$ in different regions

theta_deg_S3 = trimf(concatTable(:,3), [-60, -40, -20]);  % deg of theta in S3
theta_deg_S2 = trimf(concatTable(:,3), [-33, -20, -7]);   % deg of theta in S2
theta_deg_S1 = trimf(concatTable(:,3), [-14, -7, 0]);     % deg of theta in S1
theta_deg_CE = trimf(concatTable(:,3), [-4, 0, 4]);       % deg of theta in CE
theta_deg_B1 = trimf(concatTable(:,3), [0, 7, 14]);       % deg of theta in B1
theta_deg_B2 = trimf(concatTable(:,3), [7, 20, 33]);      % deg of theta in B2
theta_deg_B3 = trimf(concatTable(:,3), [20, 40, 60]);     % deg of theta in B3
%% 
% Similar to $x$ and $\phi$, we define the degree of $\theta$ and assign $\theta$ 
% to the region with maximum degree.

[theta_deg, theta_reg] = max([theta_deg_S3, theta_deg_S2, theta_deg_S1, ...
    theta_deg_CE, theta_deg_B1, theta_deg_B2, theta_deg_B3], [], 2);
%% 
% 
% 
% Since we have assigned every input-output pair $(\:x^{(i)}, \:\phi^{(i)}\:; 
% \:\theta^{(i)}\:)$ to some region, lets say (B1, CE; S1). Then we define the 
% corresponding rule to be
% 
% *IF* $x$ *is B1 AND* $\phi$ *is CE, THEN* $\theta$ *is S1.*
% 
% Now, let us generate rules from all input-output pairs.

x_regions = ["S2", "S1", "CE", "B1", "B2"];
phi_regions = ["S3", "S2", "S1", "CE", "B1", "B2", "B3"];
theta_regions = ["S3", "S2", "S1", "CE", "B1", "B2", "B3"];

% We represent the above rule in code as 
% [index(B1) index(CE) index(S1)] = [4 4 3]
% where index represents the numbering of the 
% region (from left to right) in the variable's fuzzy regions

rules = [x_reg, phi_reg, theta_reg];

% 
% Step 3 - Assign a Degree to each rule
% Now, in the above generated rules, there will be many redundant and conflicting 
% rules.  We define a conflict group as all rules which have the same "IF" part. 
% We also define the *degree of a rule as the product of the degrees of all the 
% inputs and outputs*. As proposed in the paper, this problem can be solved *by 
% accepting only one rule from a conflict group, which has the* *maximum degree*. 
% 
% 

rules_deg = x_deg .* phi_deg .* theta_deg;    % calculate degree of rule
N = length(rules_deg);

for i=1:N
    for j=1:N
        if((rules(i,1) == rules(j,1)) && (rules(i,2) == rules(j,2)))
            if(rules_deg(i) >= rules_deg(j))
                rules_deg(j) = rules_deg(i);
                rules(j,3) = rules(i,3);
            else
                rules_deg(i) = rules_deg(j);
                rules(i,3) = rules(j,3);
            end
        end
    end
end

% Now, we select the unique rules in sorted order
reduced_rules_with_deg = unique([rules rules_deg], 'rows', 'sorted');

%Now, lets print all the rules along with their degrees.
for i=1:length(reduced_rules_with_deg)
    fprintf("IF x is %s AND phi is %s, THEN theta is %s\n", ...
        x_regions(reduced_rules_with_deg(i,1)), ...
        phi_regions(reduced_rules_with_deg(i,2)), ...
        phi_regions(reduced_rules_with_deg(i,3)));
end
%% 
% Now, let us print the obtained rules as a *fuzzy rule base table*, where row 
% headers represent the different $\phi$-regions, and column headers represent 
% the different $x$-regions.

n = length(phi_regions);
m = length(x_regions);
tablecell = cell(n, m);

for i=1:n
    for j=1:m
        tablecell(i,j) = cellstr("X");
    end
end

for l=1:length(reduced_rules_with_deg)
    tablecell(reduced_rules_with_deg(l,2), ...
        reduced_rules_with_deg(l,1)) = cellstr( ...
        theta_regions(reduced_rules_with_deg(l,3)));
end

resTable = cell2table(tablecell, ...
    "RowNames", phi_regions, ...
    "VariableNames", x_regions);

disp(resTable);
% 
% Step 4 - Create a Combined Fuzzy Rule Base
% This step is redundant in this problem, because in the paper (on page 8) it 
% was assumed that *there are no linguistic rules*. Hence the final fuzzy rule 
% base is same as the one obtained above.
% 
% Step 5 - Determine a Mapping based on the Combined Fuzzy Rule Base
% Now, having created the combined fuzzy rule base, our training part is complete. 
% Given any new $(x_1, x_2)$, we can determine the output control $y$ using the 
% following defuzzification strategy. Let us try the strategy on the inputs mentioned 
% in the paper:
% 
% $$(x, \phi\degree) = (3, -30)$$
% 
% $$(x, \phi\degree) = (10, 220)$$
% 
% $$(x, \phi\degree) = (13, 30)$$
% 
% Given these initial state of the truck, we use the combined fuzzy rule base 
% obtained in Step 4 to route the truck to the landing dock state $(x_f, \phi_f\degree) 
% = (10, 90)$. 
% 
% 
% 
% First, for each fuzzy rule, we combine the degrees of the antecedents of $x$ 
% and $\phi$ in the corresponding input regions to determine the degree of the 
% output control.
% 
% Then, we use the *centroid defuzzification formula,* mentioned in the paper, 
% to determine the output control $\theta$. Once we get the output control value, 
% we can use the approximate truck kinematics equations (13)-(15) to get the new 
% values of $x$ and $\phi$. We keep on repeating this procedure for some number 
% of iterations until the truck reaches the landing dock state.

x_test = [3 10 13];         % test inputs' initial x coordinates
phi_test = [-30, 220, 30];  % test inputs' initial phi coordinates
iter = 100;                  % maximum number of iterations
%% 
% Now, we define the centroids of all the output regions.

theta_centroids = [-40 -20 -7 0 7 20 40];   % centroids of output regions
theta_bar = theta_centroids(reduced_rules_with_deg(:,3)); % theta centroids for rules
%% 
% Though the $y$ coordinate is not considered as an input for the purpose of 
% driving the truck, however for the purpose of trajectory plotting, let us assume 
% the following values for initial$y$ coordinate of the 3 test inputs (values 
% taken approximately from figure in the paper).

y_test = [7 10 5];          % test inputs' initial y coordinates
b = 4;                      % length of truck
%% 
% Now let us add the following linguistic rules to the rule base, and repeat 
% Step 4 and Step 5.
% 
% IF x is CE and phi is S1, THEN theta is S2.
% 
% IF x is CE and phi is CE, THEN theta is CE.
% 
% IF x is CE and phi is B1, THEN theta is B2.
% 
% IF x is CE and phi is B2, THEN theta is B3.
% Step 4

linguistic_rules = [3 3 2 1; 3 4 4 1; 3 5 6 1; 3 6 7 1];
newrules = [reduced_rules_with_deg; linguistic_rules];
theta_bar = theta_centroids(newrules(:,3));
iter = 80;
% Step 5

figure;
axis([0 20 -10 100]);
trajectory_cols = ["Red", "Blue", "Green"];
for j=1:3
    x_final = x_test(j);
    phi_final = phi_test(j);
    y_final = y_test(j);
    plot(x_final,y_final,'.-','MarkerSize',9.0, ...
        'Color', trajectory_cols(j));
    title('Truck trajectories using numerical-fuzzy controller');
    xlim([0 20]);
    ylim([0 100]);
    xlabel('x (in meters)');
    ylabel('y (in meters)');
    xline = [x_final x_final - b*cosd(phi_final)];
    yline = [y_final y_final - b*sind(phi_final)];
    % disp([xline, yline]);
    % line(xline, yline, "Color", trajectory_cols(j));
    hold on;
    for itr=1:iter
        if(isnan(x_final) || isnan(phi_final))
            break
        end
        % Degrees of x in different regions
        x_test_deg_S2 = trapmf(x_final, [0 0 1.5 7]);       % deg of x in S2
        x_test_deg_S1 = trimf(x_final, [4, 7, 10]);         % deg of x in S1
        x_test_deg_CE = trimf(x_final, [9, 10, 11]);        % deg of x in CE
        x_test_deg_B1 = trimf(x_final, [10, 13, 16]);       % deg of x in B1
        x_test_deg_B2 = trapmf(x_final, [13, 18.5, 20, 20]);% deg of x in B2

        x_test_deg = [x_test_deg_S2 x_test_deg_S1 x_test_deg_CE x_test_deg_B1 x_test_deg_B2];

        % Degrees of phi in different regions
        phi_test_deg_S3 = trimf(phi_final, [-115, -65, -15]); % deg of phi in S3
        phi_test_deg_S2 = trimf(phi_final, [-45, 0, 45]);     % deg of phi in S2
        phi_test_deg_S1 = trimf(phi_final, [15, 52.5, 90]);   % deg of phi in S1
        phi_test_deg_CE = trimf(phi_final, [80, 90, 100]);    % deg of phi in CE
        phi_test_deg_B1 = trimf(phi_final, [90, 127.5, 165]); % deg of phi in B1
        phi_test_deg_B2 = trimf(phi_final, [135, 180, 225]);  % deg of phi in B2
        phi_test_deg_B3 = trimf(phi_final, [195, 245, 295]);  % deg of phi in B3

        phi_test_deg = [phi_test_deg_S3 phi_test_deg_S2 phi_test_deg_S1 phi_test_deg_CE phi_test_deg_B1 phi_test_deg_B2 phi_test_deg_B3];
        
        theta_test_deg = x_test_deg( ...
            newrules(:,1)).*phi_test_deg( ...
            newrules(:,2));
        theta_final = sum(theta_test_deg.*theta_bar)/sum(theta_test_deg);
        x_final = x_final + cosd(phi_final + theta_final) + sind( ...
            theta_final)*sind(phi_final);
        y_final = y_final + sind(phi_final + theta_final) - sind( ...
            theta_final)*cosd(phi_final);
        phi_final = phi_final - asind(2*sind(theta_final)/b);
        %disp([x_final y_final phi_final]);
        plots(:,j) = plot(x_final,y_final,'.-','MarkerSize',8.0, ...
            'Color', trajectory_cols(j), 'DisplayName', num2str(j));
        xline = [x_final x_final - b*cosd(phi_final)];
        yline = [y_final y_final - b*sind(phi_final)];
        % disp([xline, yline]);
        %line(xline, yline, "Color", trajectory_cols(j));
        xlim([0 20]);
        ylim([0 100]);
        xlabel('x (in meters)');
        ylabel('y (in meters)');
        drawnow;
        hold on;
    end
end
leglines = [plots(1,1) plots(1,2) plots(1,3)];
legend(leglines, 'Truck 1', 'Truck 2', 'Truck 3');
saveas(gcf,'trucktrajectory3.png');
%% 
% Now, we see that the new combined rule base is enough to drive the truck to 
% the loading dock position.
%% 
% 
% 
% 
% 
%