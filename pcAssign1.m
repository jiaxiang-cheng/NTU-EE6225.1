%% Process Control Assignment 1 - written by CHENG Jiaxiang

%*************************************************************************%
%  Before running the following codes, you need to run the simulink       %
%  firstly. Then the step response data of the process can be obtained.   %
%                                                                         %
%  Step 1: Clear all                                                      %
%  Step 2: Run originalProcess.slx                                        %
%  Step 3: Run pcAssign1.m                                                % 
%  Step 4: Run resultsTesting.slx                                         %
%  ( you can also test the Nyquist results by nyquistTest.m )             %  
%*************************************************************************%

%% First order plus time delay in TIME domain (using LSM)

% here to get the data from step response
data(:,1) = out.logsout{1}.Values.time % col 1 is time series
data(:,2) = out.logsout{1}.Values.data % col 2 is output

% here to get the data from relay feedback
data(:,3) = out.logsout{3}.Values.time % col 1 is time series
data(:,4) = out.logsout{3}.Values.data % col 2 is output

% here to get the data from sustained oscillation
data(:,5) = out.logsout{2}.Values.time % col 1 is time series
data(:,6) = out.logsout{2}.Values.data % col 2 is output

N = 100; % the num of sampling data for LSM calculation
A = 1;  % the amplitude of the step signal

% calculate the integration at every sampled time
integ = zeros(N,1);
for j = 1 : N
    for i = 1 : j
        h = (data(i+1,2) + data(i,2))/2;
        Ts = data(i+1,1) - data(i,1); % integration step size
        integ(j,1) = integ(j,1) + h*Ts;
    end
end      % using Midpoint Rule

% using the least square method to obtain the results
for k = 1 : N
    Y(k,1) = -integ(k,1);
    Y(k,2) = -A;
    Y(k,3) = data(k,1)*A;
end   % LSM matrix 1
T = data(1:N,2);   % LSM matrix 2
theta = inv(Y'*Y)*Y'*T; % obtain the final matrix Theta

% calculate the final parameters for first order plus delay
a = theta(1,1); b = theta(3,1);
L_1 = theta(2,1)/theta(3,1);  % time delay constant L
Kp_1 = b/a;  % gain of the transfer function
Ts_1 = 1/a;  % time constant for the first order
%*****************************************************%
%  Final result 1:                                    %
%                      1.052153                       %
%           G1(s) = ！！！-！！！！ e^(-1.593855s)    %
%                    7.021793s + 1                    %
%*****************************************************%

%% First order plus time delay in FREQUENCE domain (using LSM)

tf = N;
Yinf = 1; Ytf = data(tf,2);

Num2 = 5;
for i = 1 : Num2
    integR = 0;
    integI = 0;
    wi = 10^(i-Num2);
    for t = 1 : tf
        alpha = data(t+1,1)-data(t,1);
        deltaFtf = data(t,2)-Yinf;
        integR = integR + alpha * deltaFtf * sin(wi*data(t,1));
        integI = integI + alpha * deltaFtf * cos(wi*data(t,1));
    end
    Gw(i,3) = wi;
    Gw(i,4) = Yinf + wi*integR;
    Gw(i,5) = wi*integI;
    Gw(i,1) = sqrt(Gw(i,4)^2 + Gw(i,5)^2);
    Gw(i,2) = atan2(Gw(i,5),Gw(i,4));
end
for i = 1 : Num2
    wi = 10^(i - Num2);
    Yw(i,1) = - wi^2 * Gw(i,1)^2;
    Yw(i,2) = 1;
    Tw(i,1) = Gw(i,1)^2;
end

thetaw = inv(Yw'*Yw)*Yw'*Tw;
aw = sqrt(thetaw(1,1));
bw = sqrt(thetaw(2,1));
Lw = 0;
for i = 1 : Num2
    Lw = Lw + (-Gw(i,2)-atan(aw*Gw(i,3)))/Gw(i,3);
end
L_2 = Lw/Num2;
Kp_2 = bw;
Ts_2 = aw;
%*****************************************************%
%  Final result 2:                                    %
%                      1.003909                       %
%           G2(s) = ！！！-！！！！ e^(-4.301368s)    %
%                    2.200410s + 1                    %
%*****************************************************%

%% First order plus time delay using sustained oscillation

h_3 = 0.01; % Set up in the relay feedback

% Find the time delay L in oscillation:
i = 2; 
while (data(i,4) - data(i-1,4) == 0)
    i = i + 1;
end
L_3 = data(i,3);

% Find the period Pu in oscillation:
flag = 1; i = 0; num_max = 6;
while (flag == 1)
    for j = 2 : N
        if (data(j,6)>data(j+1,6))&&(data(j,6)>data(j-1,6))
            i = i + 1;
            peak (1,i) = data(j,5);
            peak (2,i) = data(j,6);
        else if (data(j,6)<data(j+1,6))&&(data(j,6)<data(j-1,6))
                i = i + 1;
                peak (1,i) = data(j,5);
                peak (2,i) = data(j,6);
            end
        end
        if (i == num_max + 1)
            flag = 0;
        end
    end
end
for n = 2 : num_max
    half_T(1,n-1) = peak(1,n+1) - peak(1,n);   
end
Pu = half_T(1,num_max-1)*2;
a_3 = peak(2,2); % Magnitude of oscillation

% Calculate the ultimate frequency and gain:
Ku_3 = 4 * h_3/ (a_3 * pi);
wu_3 = 2 * pi/ Pu;

% Calculate the parameters in first order plus delay:
Kp_3 = Kp_2; % Adopt the calculated gain in method 2
Ts_3 = sqrt((Kp_3*Ku_3)^2-1)/wu_3;
%*****************************************************%
%  Final result 3:                                    %
%                      1.003909                       %
%           G3(s) = ！！！-！！！！ e^(-3.0s)         %
%                    2.520361s + 1                    %
%*****************************************************%

%% Results Display

fprintf('TIME domain via LSM:\n')
fprintf('delay L = %f\ngain Kp = %f\ntime constant Ts = %f\n',L_1,Kp_1,Ts_1)
fprintf('\n')
fprintf('FREQUENCY domain via LSM:\n')
fprintf('delay L = %f\ngain Kp = %f\ntime constant Ts = %f\n',L_2,Kp_2,Ts_2)
fprintf('\n')
fprintf('Relay Feedback:\n')
fprintf('delay L = %f\ngain Kp = %f\ntime constant Ts = %f\n',L_3,Kp_3,Ts_3)


