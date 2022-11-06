%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jere Knuutinen                                                   %
% Comparison of linear optimization control and rule based control %
% 20.7.2021                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%
clear all; clearvars; clc; close all
global t_start
global t_end
global t_start_rule_based
global t_end_rule_based
%% %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants that can be changed %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
E_bandwidth = 17;                % Max bandwidth [kWh/h]
E_charge = 11;                   % Max charging power [kWh/h]
E_need = 18;                     % Needed energy for the battery per day [kWh/d]
N = 24;                          % Optimization horizon (one day)
t_start = 16;                    % Starting time when electric vechile can be charged, for example at 16.00
t_end = 8;                       % Ending time when electric vechile can not anymore be charged, for example at 8.00
whole_day = 0;                   % If you want to simulate by assuming that electric vehicle is alwayes ready to be charged, 1 == car is not available at 8-16, 0 == car is alwayes availabel
Year = 2018;                     % Select 2017, 2018 or 2019
t_start_rule_based = 14;
t_end_rule_based = 18; 

[total_price_optimization, total_price_reference, prob] = main(E_bandwidth, E_charge, E_need, N, whole_day, Year)

function [total_price_optimization, total_price_reference, prob] = main(E_bandwidth, E_charge, E_need, N, whole_day, Year)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Constants that cannot be changed %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Year == 2017
        Tax = 2.7937/100;                 % Tax [€/kWh]
        Trans = 2.780/100;                % Transmission, energy [€/kWh]
        MarginSpot = 0.220/100;           % Energy, marginal to SPOT [€/kWh]
        MarginSellingSpot = 0.200/100;    % Selling, marginal to SPOT [€/kWh]
        SellingTrans = 0.0868/100;        % Selling, transmission [€/kWh]
        Alv = 1.24;                       % Value added tax
    elseif Year == 2018
        Tax = 2.7937/100;                 % Tax [€/kWh]
        Trans = 3.063/100;                % Transmission, energy [€/kWh]
        MarginSpot = 0.220/100;           % Energy, marginal to SPOT [€/kWh]
        MarginSellingSpot = 0.200/100;    % Selling, marginal to SPOT [€/kWh]
        SellingTrans = 0.0868/100;        % Selling, transmission [€/kWh]
        Alv = 1.24;                       % Value added tax   
    else
        Tax = 2.7937/100;                 % Tax [€/kWh]
        Trans = 3.640/100;                % Transmission, energy [€/kWh]
        MarginSpot = 0.220/100;           % Energy, marginal to SPOT [€/kWh]
        MarginSellingSpot = 0.200/100;    % Selling, marginal to SPOT [€/kWh]
        SellingTrans = 0.0868/100;        % Selling, transmission [€/kWh]
        Alv = 1.24;                       % Value added tax 
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Load, PV and Spot data and create timetable called 'information' %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    data = readtable('Datat.xlsx'); 
    [Information] = ReadFiles(data, Year);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Form grid price and PV price %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Information] = FormElectricityPrices(Information, Tax, Trans, MarginSpot,MarginSellingSpot, SellingTrans, Alv);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Editing table: We add boolean variables to the table. The boolean %
    % variable incates if the electic vechile can be charged.           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Information] = CreateBoolean(Information);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Simulate optimization and reference controls %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Information, prob] = simulate(Information, N, E_need, E_charge, E_bandwidth, whole_day);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate price of the control methods %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [total_price_optimization, total_price_reference] = calculate_prices(Information);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate self consumption %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [self_consumption_optimization, self_consumption_reference] = calculate_self_consumption(Information);

    %%%%%%%%%%%%%%%%
    % Plot figures %
    %%%%%%%%%%%%%%%%
    plot_figures(Information)
end

function [Information, prob] = simulate(Information, N, E_need, E_charge, E_bandwidth, whole_day)
    E_grid_load_vec = [];   % Store energy flow frid to load optimization control
    E_pv_load_vec = [];     % Store energy flow from PV to load optimization control
    E_pv_battery_vec = [];  % Store energy flow from PV to battery optimization control
    E_grid_battery_vec = [];% Store energy flow from grid to battery optimization control
    E_pv_grid_vec = [];     % Store energy flow from PV to grid optimization control

    E_grid_load_vec_ref = [];   % Store energy flow frid to load refernce control
    E_pv_load_vec_ref = [];     % Store energy flow from PV to load refernce control
    E_pv_battery_vec_ref = [];  % Store energy flow from PV to battery refernce control
    E_grid_battery_vec_ref = [];% Store energy flow from grid to battery refernce control
    E_pv_grid_vec_ref = [];     % Store energy flow from PV to grid refernce control
    
    % Year long simulation
    for i = N:N:365*N

        % Call optimization control function
        [E_grid_load, E_pv_load, E_pv_battery, E_grid_battery, E_pv_grid, prob] = optimization_control(N,Information.PV((i+1)-N:i),Information.BaseLoad((i+1)-N:i), ...
            Information.P_Grid_electricity_buy((i+1)-N:i),...
             Information.P_PV_electricity_buy((i+1)-N:i),Information.P_PV_electricity_sell((i+1)-N:i), ...
             E_bandwidth, E_charge, E_need, Information.boolean_flag((i+1)-N:i), whole_day);
        
         % Call reference control function
        [E_grid_load_ref, E_pv_load_ref, E_pv_battery_ref, E_grid_battery_ref, E_pv_grid_ref] = reference_control(Information.PV((i+1)-N:i), ...
            Information.BaseLoad((i+1)-N:i), E_charge, E_need, E_bandwidth, Information.boolean_flag2((i+1)-N:i));
        
        % Store daily data to the vectors
        E_grid_load_vec = [E_grid_load_vec; E_grid_load]; 
        E_pv_load_vec = [E_pv_load_vec; E_pv_load];
        E_pv_battery_vec = [E_pv_battery_vec; E_pv_battery];
        E_grid_battery_vec = [E_grid_battery_vec; E_grid_battery];
        E_pv_grid_vec = [E_pv_grid_vec; E_pv_grid];

        E_grid_load_vec_ref = [E_grid_load_vec_ref; E_grid_load_ref']; 
        E_pv_load_vec_ref = [E_pv_load_vec_ref; E_pv_load_ref'];
        E_pv_battery_vec_ref = [E_pv_battery_vec_ref; E_pv_battery_ref'];
        E_grid_battery_vec_ref = [E_grid_battery_vec_ref; E_grid_battery_ref'];
        E_pv_grid_vec_ref = [E_pv_grid_vec_ref; E_pv_grid_ref'];
    end
    
    
    Information.E_grid_load = E_grid_load_vec;
    Information.E_pv_load = E_pv_load_vec;
    Information.E_pv_battery = E_pv_battery_vec;
    Information.E_grid_battery = E_grid_battery_vec;
    Information.E_pv_grid = E_pv_grid_vec;

    Information.E_grid_load_ref = E_grid_load_vec_ref;
    Information.E_pv_load_ref = E_pv_load_vec_ref;
    Information.E_pv_battery_ref = E_pv_battery_vec_ref;
    Information.E_grid_battery_ref = E_grid_battery_vec_ref;
    Information.E_pv_grid_ref = E_pv_grid_vec_ref;
end

function [total_price_optimization, total_price_reference] = calculate_prices(Information)
    % price of cost optimization control
    price_PV = [];
    price_grid = [];
    price_PV_sell = [];
    for i = 1:height(Information)
        price_PV(i) = Information.E_pv_load(i)*Information.P_PV_electricity_buy(i) + Information.E_pv_battery(i)*Information.P_PV_electricity_buy(i);
        price_grid(i) = Information.E_grid_load(i)*Information.P_Grid_electricity_buy(i) + Information.E_grid_battery(i)*Information.P_Grid_electricity_buy(i);
        price_PV_sell(i) = Information.E_pv_grid(i)*Information.P_PV_electricity_sell(i);  
    end

    % price of cost reference control
    price_PV_ref = [];
    price_grid_ref = [];
    price_PV_sell_ref = [];
    for i = 1:height(Information)
        price_PV_ref(i) = Information.E_pv_load_ref(i)*Information.P_PV_electricity_buy(i) + Information.E_pv_battery_ref(i)*Information.P_PV_electricity_buy(i);
        price_grid_ref(i) = Information.E_grid_load_ref(i)*Information.P_Grid_electricity_buy(i) + Information.E_grid_battery_ref(i)*Information.P_Grid_electricity_buy(i);
        price_PV_sell_ref(i) = Information.E_pv_grid_ref(i)*Information.P_PV_electricity_sell(i);  
    end

    total_price_optimization = sum(price_grid)+ sum(price_PV) - sum(price_PV_sell);
    total_price_reference = sum(price_grid_ref)+ sum(price_PV_ref) - sum(price_PV_sell_ref);
end


function [self_consumption_optimization, self_consumption_reference] = calculate_self_consumption(Information) 
    self_consumption_optimization = (sum(Information.E_pv_battery) + sum(Information.E_pv_load))/(sum(Information.PV))
    self_consumption_reference = (sum(Information.E_pv_battery_ref) + sum(Information.E_pv_load_ref))/(sum(Information.PV))
end

function [E_grid_load, E_pv_load, E_pv_battery, E_grid_battery, E_pv_grid, prob] = optimization_control(N, PV, Load, P_Grid_electricity_buy, P_PV_electricity_buy, P_PV_electricity_sell, E_bandwidth, E_charge, E_need, boolean_flag, whole_day)
    
    % Minimize total cost of energy in household with PV forecast, electric vehicle that have to charged, load forecast and grid connection 
    % In optimization linear programming is used
    %   Inputs:
    %       N                       - Optimization horizon
    %       PV                      - PV daily forecast
    %       Load                    - Load daily forecast
    %       P_Grid_electricity_buy  - Electricity price when purchased from the grid
    %       P_PV_electricity_buy    - Electricity price when purchased from the PV system
    %       P_PV_electricity_sell   - Profit when own PV production sell to the grid
    %       E_bandwidth             - Max bandwidth
    %       E_charge                - Max charging energy
    %       E_need                  - Daily energy that battery needs
    %       Boolean_flag            - Indicates when charging of electric vehicle is possible
    
    %   Outputs:
    %       E_grid_load             - Hourly energy flow from grid to load
    %       E_pv_load               - Hourly energy flow from PV to load
    %       E_pv_battery            - Hourly energy flow from PV to battery
    %       E_grid_battery          - Hourly energy flow from grid to battery
    %       E_pv_grid               - Hourly energy flow from PV to grid
    
    % Use a problem based approach
    prob = optimproblem;
    
    % Decision variables
    E_pv_load_V = optimvar('E_pv_load_V',N);            % Energy flow from PV to load
    E_pv_battery_V = optimvar('E_pv_battery_V',N);      % Energy flow from PV to Battery
    E_grid_load_V = optimvar('E_grid_load_V',N);        % Energy flow from grid to load
    E_grid_battery_V = optimvar('E_grid_battery_V',N);  % Energy flow from grid to battery
    E_pv_grid_V = optimvar('E_pv_grid_V',N);            % Energy flow from PV to grid
    
    % Objective function
    prob.ObjectiveSense = 'minimize';
    prob.Objective = (E_pv_load_V' + E_pv_battery_V')*P_PV_electricity_buy + (E_grid_load_V' + E_grid_battery_V')*P_Grid_electricity_buy - E_pv_grid_V'*P_PV_electricity_sell;
    
    % constrains for the optimization
    prob.Constraints = optimconstr(N);
    
    % Energy from the grid can not exceed the bandwidth
    prob.Constraints.bandWidth = E_grid_load_V + E_grid_battery_V <= E_bandwidth;
    
    % Energy from the grid to the base load can not exceed bandwidth or be less than zero
    prob.Constraints.Grid_load_min = E_grid_load_V >= 0;
    prob.Constraints.Grid_load_max = E_grid_load_V <= E_bandwidth;
    
    if whole_day == 1
        % boolean flag gives constains when the vechile can be charged
        [time_vec] = find(boolean_flag==0);
        t1 = time_vec(1);
        t2 = time_vec(end);
        prob.Constraints.Boolean = optimconstr(N);
        prob.Constraints.Boolean = sum(E_grid_battery_V(t1:t2)) + sum(E_pv_battery_V(t1:t2)) == 0;
    end
    
    
    % Final energy of the battery have be defined
    prob.Constraints.Battery = optimconstr(N);
    prob.Constraints.Battery = sum(E_grid_battery_V(1:N)) + sum(E_pv_battery_V(1:N)) == E_need;
    
    % Energy from the grid to the battery can not exceed maximum charging limit or be less than zero
    prob.Constraints.Grid_battery_min = E_grid_battery_V >= 0;
    prob.Constraints.Grid_battery_max = E_grid_battery_V <= E_charge;
    
    % Energy from PV and grid must meet at least the baseload
    prob.Constraints.loadBalance = E_pv_load_V + E_grid_load_V == Load;
    
    % Min and max limits for the PV 
    prob.Constraints.PVmin = E_pv_load_V + E_pv_battery_V + E_pv_grid_V >= 0;
    prob.Constraints.PVmax = E_pv_load_V + E_pv_battery_V + E_pv_grid_V <= PV;
    
    % Limits from PV to the grid 
    prob.Constraints.PV_grid_min = E_pv_grid_V >= 0;
    prob.Constraints.PV_grid_max = E_pv_grid_V <= PV;
    
    % Limits from PV to the base load 
    prob.Constraints.PV_load_min = E_pv_load_V >= 0;
    prob.Constraints.PV_load_max = E_pv_load_V <= PV;
    
    % Limits from PV to the battery 
    prob.Constraints.PV_battery_min = E_pv_battery_V >= 0;
    prob.Constraints.PV_battery_max = E_pv_battery_V <= E_charge;
    
    %show(prob);
    
    % Solve the linear programming problem
    options = optimoptions(prob.optimoptions,'Display','none');
    [values,~,exitflag] = solve(prob,'Options',options);
  
    % Parse optmization results
    if exitflag <= 0
        E_grid_load = zeros(N,1);
        E_pv_load = zeros(N,1);
        E_pv_battery = zeros(N,1);
        E_grid_battery = zeros(N,1);
        E_pv_grid = zeros(N,1);
    else
        E_grid_load = values.E_grid_load_V;
        E_pv_load = values.E_pv_load_V;
        E_pv_battery = values.E_pv_battery_V;
        E_grid_battery = values.E_grid_battery_V;
        E_pv_grid = values.E_pv_grid_V;
    end
end

function [E_grid_load, E_pv_load, E_pv_battery, E_grid_battery, E_pv_grid] = reference_control(PV, Load, E_charge, E_need, E_bandwidth, boolean_flag)
%     global t_start_ruled_based
%     global t_end_ruled_based
%     boolean_flag = zeros(1,24);
%     boolean_flag(t_start_ruled_based+1:t_end_ruled_based) = 1;
    E_batt  = 0;
    E_grid_load = [];
    E_pv_load = [];
    E_pv_battery = [];
    E_grid_battery = [];
    E_pv_grid = [];
    
    for i = 1:length(PV)
       E_grid_load(i) = 0;      % Electricty purchased from the grid to load
       E_pv_load(i) = 0;        % PV electricity from PV to laod
       E_pv_battery(i) = 0;     % PV to battery
       E_grid_battery(i) = 0;   % Grid to battery
       E_pv_grid(i) = 0;        % PV to grid
       % Charging of vechile is possible
        if boolean_flag(i) == 1
           % Battery is full
            if E_batt == E_need
                % PV production is more than laod
                if PV(i) >= Load(i)
                    E_pv_load(i) = Load(i);
                    E_pv_grid(i) = PV(i) - E_pv_load(i);
                % Pv production is less than the load
                else
                    E_pv_load(i) = PV(i);
                    E_grid_load(i) = Load(i) - E_pv_load(i);
                end
           % Battery is not full yet
            else
               % Own PV production is higher than the baseload
                if PV(i) >= Load(i)
                    E_pv_load(i) = Load(i); % because PV production is more than load. 
                    % Remaining PV is more than max charging power
                    if (PV(i) - Load(i)) >= E_charge
                        % Battery's still needed energy is more than max
                        % charging capcity
                        if (E_need - E_batt) >= E_charge
                           E_pv_battery(i) = E_charge;                          % Charge with full capcity
                           E_pv_grid(i) = PV(i) - Load(i) - E_pv_battery(i);    % Sell rest to the grid
                           E_batt = E_batt + E_pv_battery(i);                 % Increment the energy of the battery
                        % Battery's still needed energy is less than max
                        % charging capcity
                        else
                            E_pv_battery(i) = E_need - E_batt;                  % Charge with still needed energy
                            E_pv_grid(i) = PV(i) - Load(i) - E_pv_battery(i);   % Sell rest to the grid
                            E_batt = E_batt + E_pv_battery(i);                  % Increment the eenrgy of the battery
                        end
                        
                    % Remaining PV is less than max charging power
                    else
                        
                        if (E_need - E_batt >= E_charge)
                            E_pv_battery(i) = PV(i) - Load(i);
                            E_grid_battery(i) = E_charge - E_pv_battery(i);
                            E_batt = E_batt + E_pv_battery(i) + E_grid_battery(i);

                        else
                                % remaining PV is more than still needed energy
                                if (PV(i) - Load(i)) >= (E_need - E_batt)
                                   E_pv_battery(i) = E_need - E_batt;                   % Still needed energy is less than remaining PV
                                   E_pv_grid(i) = PV(i) - Load(i) - E_pv_battery(i);    % Sell rest to the grid
                                   E_batt = E_batt + E_pv_battery(i);                   % Incremtn the enrgy of the battery
                                % Remaining PV is less than still needed energy
                                else
                                    E_pv_battery(i) = PV(i) - Load(i);                          % Use remaining PV to battery
                                    E_grid_battery(i) = (E_need - E_batt) - E_pv_battery(i);    % take rest from the grid
                                    E_batt = E_batt + E_pv_battery(i) + E_grid_battery(i);      % Increment the enrgy of the battery
                                end
                        end
                        
                    end
               
               % Own PV production is less than the baseload
                else
                    E_pv_load(i) = PV(i); % All PV goes to the load
                    E_grid_load(i) = Load(i) - E_pv_load(i); % Rest is taken from the grid
                    
                    % Bandwidth does not limit full charging
                    if (E_charge + E_grid_load(i)) <= E_bandwidth
                        % Remaining energy of the battery does not limit full charging
                        if (E_need - E_batt) >= E_charge
                            E_grid_battery(i) = E_charge;           % Battery charging with full capacity
                            E_batt = E_batt + E_grid_battery(i);    % Increment energy of the battery
                        % Remaining energy of the battery limits full charging
                        else
                            E_grid_battery(i) = E_need - E_batt;    % Battery charging with needed energy
                            E_batt = E_batt + E_grid_battery(i);    % Increment energy of the battery
                        end
                    
                    % Bandwidth limits full charging
                    else
                        % Battery's still needed energy is smaller than max
                        % charging capcity that is possible
                        if (E_need - E_batt) <= (E_bandwidth - E_grid_load)
                            E_grid_battery(i) = E_need - E_batt;    % Battery charging with still needed energy
                            E_batt = E_batt + E_grid_battery(i);    % Increment energy of the battery
                        % Battery's still needed energy is more than max
                        % charging capacity that is possible
                        else
                            E_grid_battery(i) = E_bandwidth - E_grid_load(i);   % Battery charging with the energy that is the largest possible
                            E_batt = E_batt + E_grid_battery(i);                % Increment energy of the battery
                        end
                    end
                end
            end
       % Charging of vechile is not possible
        else
           % PV production is more than laod
           if PV(i) > Load(i)
               E_pv_load(i) = Load(i);
               E_pv_grid(i) = PV(i) - E_pv_load(i);
           % Pv production is less than the load
           else
               E_pv_load(i) = PV(i);
               E_grid_load(i) = Load(i) - E_pv_load(i);
           end
        end          
    end
end

function [Information] = ReadFiles(data, Year)
    % This function convertes tables to the timetables and makes other
    % modifications too
    
    % You maybe have to change these if you read some different data file
    if Year == 2017
        Spot = data.SPOT;                   
        PV = data.Etel_ + data.It__l_nsi;
        BaseLoad = data.Kulutus;
    elseif Year == 2018
        Spot = data.SPOT_1;                   
        PV = data.Etel__1 + data.It__l_nsi_1;
        BaseLoad = data.Kulutus_1;
    else
        Spot = data.SPOT_2;                   
        PV = data.Etel__2 + data.It__l_nsi_2;
        BaseLoad = data.Kulutus_2;
    end
        

    % Read time verctor from excle (This can also maybe do in the more elegant way)
    Temp = readtable('Time.xlsx');
    Temp.Spot = Spot./1000; % Convert €/MWh to €/kWh
    Temp.BaseLoad = BaseLoad;
    Temp.PV = PV;
    
    % Convert table to timetable. Also for some reason temp table contains time
    % weird timestaps so we round them.
    temp_information = table2timetable(Temp);
   
    Information = retime(temp_information,'regular','linear','TimeStep', hours(1)); 
   
end

function Information = FormElectricityPrices(Information, Tax, Trans, MarginSpot,MarginSellingSpot, SellingTrans, Alv)
    
    Information.P_Grid_electricity_buy = Information.Spot*Alv + Tax + Trans + MarginSpot;       % Price of grid electricity
    Information.P_PV_electricity_buy = Information.Spot - SellingTrans - MarginSellingSpot;     % Price of PV electricity, This is the lost revenue when you do not sell PV electricity to grid
    Information.P_PV_electricity_sell = Information.Spot - SellingTrans - MarginSellingSpot;    % When you sell PV electricuty to grid  
    
end


function [Information] = CreateBoolean(Information)

    global t_start;
    global t_end;
    global t_end_rule_based
    
    Temp_table = Information;            % Create temp var 
    Temp_table.Day = dateshift(Information.Var1,'start','day'); % Create day variable
    Temp_table = timetable2table(Temp_table);   
    Temp_table = removevars(Temp_table,{'Spot','BaseLoad', 'PV', 'P_Grid_electricity_buy','P_PV_electricity_buy','P_PV_electricity_sell'}); % If you add more variables to the timetable remember remove them here also then.
    Clock_times = rowfun(@ClockTimes, Temp_table,'GroupingVariable','Day','OutputVariableNames',{'time1','time2'});
    Clock_times2 = rowfun(@ClockTimes_rule_based, Temp_table,'GroupingVariable','Day','OutputVariableNames',{'time1'});
    Information.boolean_flag = Information.PV*0;
    Information.boolean_flag2 = Information.PV*0;
    for i = 1:height(Clock_times)
        % Boolean flags for the optimization control
        Information.boolean_flag(timerange(Clock_times.time1(i), Clock_times.time1(i) + hours(24-t_start),'closed')) = 1;
        Information.boolean_flag(timerange(Clock_times.time2(i), Clock_times.time2(i) + hours(t_end-1),'closed')) = 1;
        Information.boolean_flag2(timerange(Clock_times2.time1(i), Clock_times2.time1(i) + hours(24-t_end_rule_based+1),'closed')) = 1; % Boolean flags for the reference control
    end
end

function [t1,t2] = ClockTimes(t,x)
    global t_start;
    global t_end;
    t1 = t(t_start + 1);
    t2 = t(t_end - (t_end-1));
end

function [t1] = ClockTimes_rule_based(t,x)
    global t_start_rule_based;
    t1 = t(t_start_rule_based+1);
  
end


function [] = plot_figures(Information)

    FileFormat = '.eps'; % '.eps' or '.tif'

    PrintCommand = '-depsc'; % EPS color

    PrintCommand = '-dtiff'; % TIFF 24-bit

    Resolution = '-r600'; % .tiff only


    FontSize = 10;
    FontName = 'Times New Roman';
    LineWidth = 1;
    MarkerSize = 2;
    Units = 'centimeters';

    Location = 'SouthWest';
    FigSizeWidth = 20;
    FigSizeHeight = 10;
    FigPosLeft = 0;
    FigPosBottom = 0;
    GridLineStyle = ':';
    %set(gcf,'Units',Units,'Position',[FigPosLeft FigPosBottom FigSizeWidth FigSizeHeight]);
        
    Information.HrOfDay = hour(Information.Var1);
    
    % Optimization control plots
    byHr_PV_battery = varfun(@sum,Information(:,{'E_pv_battery','HrOfDay'}),...
    'GroupingVariables','HrOfDay','OutputFormat','table');
    bar(byHr_PV_battery{:,{'HrOfDay'}},byHr_PV_battery{:,{'sum_E_pv_battery'}}, 'FaceColor',[1 0.85 0])
    title("Energy flow from PV to battery (Optimization control)")
        ylabel("kWh/year",'FontName',FontName,...
       'FontSize',FontSize);
    ax1 = gca;
    ax1.XTick = 0:1:23;
    ax1.XTickLabel = {'00','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23'};
    ax1.XTickLabelRotation = 45;
    set(gca,'FontName',FontName,'FontSize',FontSize,'GridLineStyle',GridLineStyle);

    figure
    byHr_grid_battery = varfun(@sum,Information(:,{'E_grid_battery','HrOfDay'}),...
    'GroupingVariables','HrOfDay','OutputFormat','table');
    bar(byHr_grid_battery{:,{'HrOfDay'}},byHr_grid_battery{:,{'sum_E_grid_battery'}}, 'FaceColor',[1 0.85 0])
    title("Energy flow from grid to battery (Optimization control)")
        ylabel("kWh/year",'FontName',FontName,...
       'FontSize',FontSize);
    ax1 = gca;
    ax1.XTick = 0:1:23;
    ax1.XTickLabel = {'00','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23'};
    ax1.XTickLabelRotation = 45;
    set(gca,'FontName',FontName,'FontSize',FontSize,'GridLineStyle',GridLineStyle);
    
    figure
    byHr_PV_grid = varfun(@sum,Information(:,{'E_pv_grid','HrOfDay'}),...
    'GroupingVariables','HrOfDay','OutputFormat','table');
    bar(byHr_PV_grid{:,{'HrOfDay'}},byHr_PV_grid{:,{'sum_E_pv_grid'}}, 'FaceColor',[1 0.85 0])
    title("Energy flow from PV to grid (Optimization control)")
    ylabel("kWh",'FontName',FontName,...
       'FontSize',FontSize);
    ax1 = gca;
    ax1.XTick = 0:1:23;
    ax1.XTickLabel = {'00','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23'};
    ax1.XTickLabelRotation = 45;
    set(gca,'FontName',FontName,'FontSize',FontSize,'GridLineStyle',GridLineStyle);
    
    figure
    byHr_PV_load = varfun(@sum,Information(:,{'E_pv_load','HrOfDay'}),...
    'GroupingVariables','HrOfDay','OutputFormat','table');
    bar(byHr_PV_load{:,{'HrOfDay'}},byHr_PV_load{:,{'sum_E_pv_load'}}, 'FaceColor',[1 0.85 0])
    title("Energy flow from PV to load (Optimization control)")
        ylabel("kWh/year",'FontName',FontName,...
       'FontSize',FontSize);
    ax1 = gca;
    ax1.XTick = 0:1:23;
    ax1.XTickLabel = {'00','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23'};
    ax1.XTickLabelRotation = 45;
    set(gca,'FontName',FontName,'FontSize',FontSize,'GridLineStyle',GridLineStyle);
    
    figure
    byHr_grid_load = varfun(@sum,Information(:,{'E_grid_load','HrOfDay'}),...
    'GroupingVariables','HrOfDay','OutputFormat','table');
    bar(byHr_grid_load{:,{'HrOfDay'}},byHr_grid_load{:,{'sum_E_grid_load'}})
    title("Energy flow from grid to load (Optimization control)")
        ylabel("kWh/year",'FontName',FontName,...
       'FontSize',FontSize);
    ax1 = gca;
    ax1.XTick = 0:1:23;
    ax1.XTickLabel = {'00','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23'};
    ax1.XTickLabelRotation = 45;
    set(gca,'FontName',FontName,'FontSize',FontSize,'GridLineStyle',GridLineStyle);
    
    % Reference control plots
    byHr_PV_battery_ref = varfun(@sum,Information(:,{'E_pv_battery_ref','HrOfDay'}),...
    'GroupingVariables','HrOfDay','OutputFormat','table');
    bar(byHr_PV_battery_ref{:,{'HrOfDay'}},byHr_PV_battery_ref{:,{'sum_E_pv_battery_ref'}}, 'FaceColor',[1 0.85 0])
    title("Energy flow from PV to battery (Reference control)")
        ylabel("kWh/year",'FontName',FontName,...
       'FontSize',FontSize);
    ax1 = gca;
    ax1.XTick = 0:1:23;
    ax1.XTickLabel = {'00','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23'};
    ax1.XTickLabelRotation = 45;
    set(gca,'FontName',FontName,'FontSize',FontSize,'GridLineStyle',GridLineStyle);
    
    figure
    byHr_grid_battery_ref = varfun(@sum,Information(:,{'E_grid_battery_ref','HrOfDay'}),...
    'GroupingVariables','HrOfDay','OutputFormat','table');
    bar(byHr_grid_battery_ref{:,{'HrOfDay'}},byHr_grid_battery_ref{:,{'sum_E_grid_battery_ref'}}, 'FaceColor',[1 0.85 0]);
    title("Energy flow from grid to battery (Reference control)")
        ylabel("kWh/year",'FontName',FontName,...
       'FontSize',FontSize);
    ax1 = gca;
    ax1.XTick = 0:1:23;
    ax1.XTickLabel = {'00','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23'};
    ax1.XTickLabelRotation = 45;
    set(gca,'FontName',FontName,'FontSize',FontSize,'GridLineStyle',GridLineStyle);
    
    figure
    byHr_PV_grid_ref = varfun(@sum,Information(:,{'E_pv_grid_ref','HrOfDay'}),...
    'GroupingVariables','HrOfDay','OutputFormat','table');
    bar(byHr_PV_grid_ref{:,{'HrOfDay'}},byHr_PV_grid_ref{:,{'sum_E_pv_grid_ref'}}, 'FaceColor',[1 0.85 0])
    title("Energy flow from PV to grid (Reference control))")
        ylabel("kWh/year",'FontName',FontName,...
       'FontSize',FontSize);
    ax1 = gca;
    ax1.XTick = 0:1:23;
    ax1.XTickLabel = {'00','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23'};
    ax1.XTickLabelRotation = 45;
    set(gca,'FontName',FontName,'FontSize',FontSize,'GridLineStyle',GridLineStyle);
    
    figure
    byHr_PV_load_ref = varfun(@sum,Information(:,{'E_pv_load_ref','HrOfDay'}),...
    'GroupingVariables','HrOfDay','OutputFormat','table');
    bar(byHr_PV_load_ref{:,{'HrOfDay'}},byHr_PV_load_ref{:,{'sum_E_pv_load_ref'}}, 'FaceColor',[1 0.85 0])
    title("Energy flow from PV to load (Reference control)")
        ylabel("kWh/year",'FontName',FontName,...
       'FontSize',FontSize);
    ax1 = gca;
    ax1.XTick = 0:1:23;
    ax1.XTickLabel = {'00','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23'};
    ax1.XTickLabelRotation = 45;
    set(gca,'FontName',FontName,'FontSize',FontSize,'GridLineStyle',GridLineStyle);
    
    figure
    byHr_grid_load_ref = varfun(@sum,Information(:,{'E_grid_load_ref','HrOfDay'}),...
    'GroupingVariables','HrOfDay','OutputFormat','table');
    bar(byHr_grid_load_ref{:,{'HrOfDay'}},byHr_grid_load_ref{:,{'sum_E_grid_load_ref'}}, 'FaceColor',[1 0.85 0])
    title("Energy flow from grid to load (Reference control)")
        ylabel("kWh/year",'FontName',FontName,...
       'FontSize',FontSize);
    ax1 = gca;
    ax1.XTick = 0:1:23;
    ax1.XTickLabel = {'00','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23'};
    ax1.XTickLabelRotation = 45;
    set(gca,'FontName',FontName,'FontSize',FontSize,'GridLineStyle',GridLineStyle);
end