clc, clear, close all
x1 = [0,    0,    0.865];  
x2 = [2.29, 0,    1.27];   
x3 = [2.29, 2.52, 0.865]; 
x4 = [0,    2.52, 1.27];

load('RF.mat'); 
%%
% -------------- Time and Coordinates ----------
% time:     IMU Measurement
% coord3:   Ground Truth of 3D Coordinates
% z,z_prev: 3D coordinates for Simulation


[coord3, radial] = get3Dcoord(x1, x2, x3, x4, rtime);

% Extract acceleration data for x, y, z axis, name by acc_data_(Axis)

%
%%
% Measurement Vector [acc_x, acc_y, acc_z, orientation, angular velocity]

% State vector [xdotdot, ydotdot, orientation, angular velocity]
%x = [0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01]'*ones(1,length(accT));
%o1 a1 o2 a2 o3 a3
x = zeros(15, length(rtime)); x(1,1) = 1.03;
% Intial P matrix
P = eye(length(x(:,1)));
% covariance for w
Q = diag([0.1, 0.1, 10, 0.1, 0.1, 10, 0.1, 0.1, 10, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,]);

%%
% covariance for v
RR = [300, 300, 300, 300, 300, 300, 300, 300]; 
factor = 1;

i = 1;
j = 1;
l = 1;
for m = 1:1:length(rtime)
    %[index, i, j, k] = getIndex(i, j, k, mag_time, gyro_time, acc_time, m, time);
 
    [y, l, index, len] = getyNPVA(r1, r2, r3, r4, rtime, l, m);

    %y
    if m == 1
        x(:,m) = x(:,m);
    else
        F = getF(rtime, m);
        %x(:,m-1)
        x(:,m) = F*x(:,m-1);
        %F
        %xx = x(:,m)
        xx = x(:,m);
        P = F*P*F' + Q;
        H = getJacoN(y,x(:,m),index,factor);

        R = getR(RR, index);
        % back tracking
        G = P*H'/(H*P*H' + R); %lambda
        % e
        x(:,m) = x(:,m) + G*(y - getH(x(:,m), index, factor));
        P = (eye(length(x(:,1))) - G*H)*P;
    end
    
    
    
    
end
%%
% figure
% subplot(311), plot(time, x(1,:), 'LineWidth', 3); hold on; plot(accT(2:end), posX, 'LineWidth', 2); hold on; plot(time, PP1, 'LineWidth', 2); legend('estimated','integration','ground truth'); title('Position Along x Axis [m]'); xlabel('t [s]'); ylabel('position [m]'); grid on; grid minor;%hold on; plot(time, XX(1,:))
% subplot(312), plot(time, x(7,:), 'LineWidth', 2); hold on; plot(accT(2:end), posY, 'LineWidth', 2); legend('estimated','integration'); title('Position Along y Axis [m]'); xlabel('t [s]'); ylabel('position [m]'); grid on; grid minor;%hold on; plot(time, XX(4,:))
% subplot(313), plot(time, x(4,:), 'LineWidth', 2); hold on; plot(accT(2:end), posZ, 'LineWidth', 2); legend('estimated','integration'); title('Position Along z Axis [m]'); xlabel('t [s]'); ylabel('position [m]'); grid on; grid minor;%hold on; plot(time, XX(7,:))
% csvwrite('data/x.csv',x);   csvwrite('data/XX.csv',XX);   csvwrite('data/phi_gt.csv', phi_gt); 
% csvwrite('data/T3.csv',x);   csvwrite('data/td3.csv',XX); 


%%
[pos1, vel1] = doubleInte(diff(rtime(1:end)), x(3,2:end));
[pos2, vel2] = doubleInte(diff(rtime(1:end)), x(6,2:end));
[pos3, vel3] = doubleInte(diff(rtime(1:end)), x(9,2:end));

N = 0;
%%
figure
subplot(9,1,1), plot(rtime(1:end-N), x(1,:), 'b', 'LineWidth', 2); hold on; plot(rtime, coord3(1,1:end), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated acceleration along x axis'); grid on;
subplot(9,1,2), plot(rtime(1:end-N), x(2,:), 'b', 'LineWidth', 2); hold on; plot(rtime, coord3(4,1:end), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated acceleration along y axis'); grid on;
subplot(9,1,3), plot(rtime(1:end-N), x(3,:), 'b', 'LineWidth', 2); hold on; plot(rtime, coord3(7,1:end), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated acceleration along z axis'); grid on;
subplot(9,1,4), plot(rtime(1:end-N), x(4,:), 'b', 'LineWidth', 2); hold on; plot(rtime, coord3(2,1:end), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated orientation along x axis'); grid on;
subplot(9,1,5), plot(rtime(1:end-N), x(5,:), 'b', 'LineWidth', 2); hold on; plot(rtime, coord3(5,1:end), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated angular velocity along x axis'); grid on;
subplot(9,1,6), plot(rtime(1:end-N), x(6,:), 'b', 'LineWidth', 2); hold on; plot(rtime, coord3(8,1:end), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated orientation along y axis'); grid on;
subplot(9,1,7), plot(rtime(1:end-N), x(7,:), 'b', 'LineWidth', 2); hold on; plot(rtime, coord3(3,1:end), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated angular velocity along y axis'); grid on;
subplot(9,1,8), plot(rtime(1:end-N), x(8,:), 'b', 'LineWidth', 2); hold on; plot(rtime, coord3(6,1:end), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated orientation along z axis'); grid on;
subplot(9,1,9), plot(rtime(1:end-N), x(9,:), 'b', 'LineWidth', 2); hold on; plot(rtime, coord3(9,1:end), 'r', 'LineWidth', 2); legend('estimated', 'gt'); title('Estimated angular velocity along z axis'); grid on;


%%
rmsErrorX = 1000*rms(x(1,:)-coord3(1,1:end-N))
rmsErrorY = 1000*rms(x(4,:)-coord3(2,1:end-N))
rmsErrorZ = 1000*rms(x(7,:)-coord3(3,1:end-N))