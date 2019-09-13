num01 = [0.08 1];
den01 = [2 1];
G01 = tf(num01, den01);
num02 = [1];
den02 = [1 1];
G02 = tf(num02, den02);
num03 = [1];
den03 = [0.4 1];
G03 = tf(num03, den03);
num04 = [1];
den04 = [0.2 1];
G04 = tf(num04, den04);
G012 = series(G01,G02);
G0123 = series(G012,G03);
G0 = series(G0123,G04);
G0.ioDelay = 3;

num1 = [1.052153];
den1 = [7.021793 1];
G1 = tf(num1, den1);
G1.ioDelay = 1.593855;

num2 = [1.003909];
den2 = [2.200410 1];
G2 = tf(num2, den2);
G2.ioDelay = 4.301368;

num3 = [1.003909];
den3 = [2.520361 1];
G3 = tf(num3, den3);
G3.ioDelay = 3.000000;

figure(1)
subplot(2,2,1), nyquist(G0); title('Original Process');
subplot(2,2,2), nyquist(G1); title('LSM in Time Domain');
subplot(2,2,3), nyquist(G2); title('LSM in Frequence Domain');
subplot(2,2,4), nyquist(G3); title('Sustained Oscillation with Relay Feedback');