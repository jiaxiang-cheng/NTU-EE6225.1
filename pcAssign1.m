%% Process Control Assignment 1 - written by CHENG Jiaxiang

%*************************************************************************%
%  Before running the following codes, you need to run the simulink       %
%  firstly. Then the step response data of the process can be obtained.   %
%  Run : pcSim.slx                                                        %
%*************************************************************************%

%% First order plus time delay in TIME domain (using LSM)

% here to get the data from step response
data(:,1) = out.logsout{1}.Values.time % col 1 is time series
data(:,2) = out.logsout{1}.Values.data % col 2 is output

N = 150; % the num of sampling data for LSM calculation
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
L = theta(2,1)/theta(3,1);  % time delay constant L
K = b/a;  % gain of the transfer function
ts = 1/a;  % time constant for the first order
%*****************************************************%
%  Final result:  ( Tf = 50; N = 80 )                 %
%                      1.034998                       %
%           G1(s) = ！！！-！！！！ e^(-1.681025s)    %
%                    6.358993s + 1                    %
%*****************************************************%
%  Final result:  ( Tf = 100; N = 150 )               %
%                      1.007976                       %
%           G1(s) = ！！！-！！！！ e^(-1.682253s)    %
%                    5.672527s + 1                    %
%*****************************************************%
% num = [0.08 1];
% den = [0.16 1.44 3.88 3.6 1];
% G = tf(num,den);
% num1 = [K]; den1 = [ts 1]; G1 = tf(num1,den1);
% figure(1) 
% subplot(1,2,1),nyquist(G1);
% subplot(1,2,2),nyquist(G);

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
Lw = Lw/Num2;

% s = tf('s');
fprintf('TIME domain via LSM:\n')
fprintf('delay L = %f\ngain K = %f\ntime cons ts = %f\n',L,K,ts)
%numt = [1];
%dent = [1 1];
%Gt = tf(numt,dent)
fprintf('FREQUENCY domain via LSM:\n')
fprintf('delay L = %f\ngain K = %f\ntime cons ts = %f\n',Lw,bw,aw)
% numf = [1];
% denf = [1 1];
% Gf = tf(numf,denf)
