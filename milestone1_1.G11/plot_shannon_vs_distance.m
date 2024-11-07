clear all;
values = [4.4857 4.8638 4.2187 4.0455 3.8821 3.5242 2.8834 2.1811  1.6253 ]*10^3;
x = [1 2 3 4 5 6 7 8 9]*2;


figure; 
plot(x,values);
xlabel('distance (cm)');
ylabel('capacity (bit/sec)')
sgtitle('Signal capacity in function of distance from speaker')