%% Setting up the model
clear all 
close all

%external currents
tau = 0.1; % 100 ms time constant

% synaptic weights
W_11 = 0.2609 * 10^-9; 
W_22 = 0.2609 * 10^-9; 
W_12 = 0.0497 * 10^-9; 
W_21 = 0.0497 * 10^-9; 

dt=0.1/1000 ;
t=[0:dt:7000/1000] ; 

a = 270 * 10^9;
b = 108;
d = 0.15;

gam = 0.641; 

trialNum = 50;

S1=zeros(trialNum, length(t));
S2=zeros(trialNum, length(t));


I1 = 0.3310 * 10^-9;
I2 = 0.3304 * 10^-9;


%% Question 1

% The Model

for trial = 1:trialNum


    I_noise1 =  0.03 * 10 ^-9 * randn(size(t)); %normal noise on inputs to population 1.
    I_noise2 =  0.03 * 10 ^-9 * randn(size(t)); % normal noise on inputs to population 2.

    S1(trial, 1) = 0.1;
    S2(trial, 1) = 0.1;

    for l=2:length(t)

        %calculating synaptic currents
        Isyn1 = (W_11 * S1(trial, l-1)) - (W_12 * S2(trial, l-1)) + I1 + I_noise1(l);
        Isyn2 = (W_22 *S2 (trial, l-1)) - (W_21 * S1(trial, l-1)) + I2 + I_noise2(l);

        %calculating firing rates
        r1 = (a * Isyn1 - b) / (1-exp(-d * ( a * Isyn1 - b))); % firing rate
        r2 = (a * Isyn2 - b) / (1-exp(-d * ( a * Isyn2 - b))); % firing rate

        %calculating activity of the populations
        dS1dt = -S1(trial, l-1)/tau + ((1-S1(trial, l-1)) * gam * r1);
        dS2dt = -S2(trial, l-1)/tau + ((1-S2(trial, l-1)) * gam * r2);
        S1(trial, l) = S1(trial, l-1) + dt * dS1dt; 
        S2(trial, l) = S2(trial, l-1) + dt * dS2dt;


    end

end

%plot
figure(1)
hold on
for i = 1:trialNum
    plot(t*1000, S1(i, :), 'b')
    plot(t*1000,  S2(i, :), 'r')

end
hold off


xlabel("Time (ms)")
ylabel("Firing rate (Hz)")
title('Allison Woodward- Timecourse of firing rates r1 (blue) and r2 (red)')


%% Question 2

frGoal = 0.5; % goal firing rate

pop1WinCount = 0; %times that population 1 reaches 0.5 first
pop2WinCount = 0;

for trial = 1:trialNum

    reached = 0; % is set to 1 when a population reaches 0.5

    for i = 2:length(t)
        
        %if pop 1 reaches 0.5 add a count to population 1 wins
        if(S1(trial, i) >= frGoal)
            pop1WinCount = pop1WinCount+1 ;
            reached = 1;

        end

        %if pop 2 reaches 0.5 add a count to population 2 wins
        if(S2(trial, i) >= frGoal)
            reached = 1;
        end

        %loop broken if either population reaches firing rate of 0.5
        if reached == 1
            break
        end

    end

end

% Display the result
fprintf('Percentage of trials where population 1 reaches a firing rate of 0.5 first: %.2f%%\n', (pop1WinCount/50)*100);

%% Question 3

reachedTime = 0;
reachedCount = 0;

for trial = 1:trialNum

    reached1 = 0; 
    reached2 = 0;

    for i = 2:length(t)
        
        %if pop 1 reaches 0.5, the time it is reached is added to total
        if(reached1 == 0 && S1(trial, i) >= frGoal)
            reachedTime = reachedTime + i ;
            reachedCount = reachedCount + 1;
            reached1 = 1;

        end

        %if pop 2 reaches 0.5 add a count to population 2 wins
        if(reached2 == 0 && S2(trial, i) >= frGoal)
            reachedTime = reachedTime + i ;
            reachedCount = reachedCount + 1;
            reached2 = 1;
        end

        %loop broken if either population reaches firing rate of 0.5
        if(reached1 == 1 && reached2 == 1)
            break
        end

    end

end

avgTimeReached = reachedTime / reachedCount / 10;

% Display the result
fprintf(['Average time it takes for either of the populations to reach the firing rate of 0.5: %d ms\n'], round(avgTimeReached));

%% Questions 4-7 are included as written responses in the word document
% as well as commented as variables below in question 8






%% Question 8

%run a loop of the original model for times with each set of inputs

tau = 0.1; % 100 ms time constant

% synaptic weights
W_11 = 0.2609 * 10^-9; 
W_22 = 0.2609 * 10^-9; 
W_12 = 0.0497 * 10^-9; 
W_21 = 0.0497 * 10^-9; 

dt=0.1/1000 ;
t=[0:dt:7000/1000] ; 

a = 270 * 10^9;
b = 108;
d = 0.15;

gam = 0.641; 

trialNum = 50;

results = zeros(4, 2);

for x = 1:4

    S1=zeros(trialNum, length(t));
    S2=zeros(trialNum, length(t));

    if x == 1
        % Question 4 external currents
        I1 = 0.3304 * 10^-9;
        I2 = 0.3304 * 10^-9;
    end

    if x == 2
        % Question 5 external currents
        I1 = 0.3305 * 10^-9;
        I2 = 0.3304 * 10^-9;
    end

    if x == 3
        % Question 6 external currents
        I1 = 0.3308 * 10^-9;
        I2 = 0.3304 * 10^-9;
    end 

    if x == 4 
        % Question 7 external currents
        I1 = 0.3310 * 10^-9;
        I2 = 0.3300 * 10^-9;
    end


    for trial = 1:trialNum



        I_noise1 =  0.03 * 10 ^-9 * randn(size(t)); %normal noise on inputs to population 1.
        I_noise2 =  0.03 * 10 ^-9 * randn(size(t)); % normal noise on inputs to population 2.

        S1(trial, 1) = 0.1;
        S2(trial, 1) = 0.1;



        for l=2:length(t)

            %calculating synaptic currents
            Isyn1 = (W_11 * S1(trial, l-1)) - (W_12 * S2(trial, l-1)) + I1 + I_noise1(l);
            Isyn2 = (W_22 *S2 (trial, l-1)) - (W_21 * S1(trial, l-1)) + I2 + I_noise2(l);

            %calculating firing rates of populations
            r1 = (a * Isyn1 - b) / (1-exp(-d * ( a * Isyn1 - b))); % firing rate
            r2 = (a * Isyn2 - b) / (1-exp(-d * ( a * Isyn2 - b))); % firing rate

            %calculating activity of populations
            dS1dt = -S1(trial, l-1)/tau + ((1-S1(trial, l-1)) * gam * r1);
            dS2dt = -S2(trial, l-1)/tau + ((1-S2(trial, l-1)) * gam * r2);

            S1(trial, l) = S1(trial, l-1) + dt * dS1dt; 
            S2(trial, l) = S2(trial, l-1) + dt * dS2dt;


        end

    end



    % calculate percentage

    frGoal = 0.5; % goal firing rate

    pop1WinCount = 0; %times that population 1 reaches 0.5 first
    pop2WinCount = 0;

    for trial = 1:trialNum

        reached = 0; % is set to 1 when a population reaches 0.5

        for i = 2:length(t)
        
            %if pop 1 reaches 0.5 add a count to population 1 wins
            if(S1(trial, i) >= frGoal)
                pop1WinCount = pop1WinCount+1 ;
                reached = 1;

            end

            %if pop 2 reaches 0.5 add a count to population 2 wins
            if(S2(trial, i) >= frGoal)

                reached = 1;
            end

            %loop broken if either population reaches firing rate of 0.5
            if reached == 1
                break
            end

        end

    end

    %calculate percent
    percent = (pop1WinCount/50)*100;

    % storing results
    results(x, 1) = I1 - I2;
    results(x, 2) = percent;

end


figure(2)
hold on
for i = 1:4
    plot(results(i,1), results(i,2),'r.', 'LineWidth', 2, 'MarkerSize', 25)
end
hold off


xlabel("Difference Between external inputs I1 and I2")
ylabel("Percent of trials I1 reaches 0.5 first")
title('Allison Woodward- Population 1 winning trials as a difference of I1 and I2')