clc; clear; close all;

%% 说明***************************************************************
%{
    只考虑phi_4p的广义坐标，假设phi_2p的运动规律是已知的匀速转动
    文献中第一组参数的仿真
    文献中的图是从曲柄在水平位置为起始点
    程序是从共线位置为起始点
    不考虑打纬阻力的外力矩
    结果：
    杆4的角速度峰值对应相等    
%}





%% 基本参数输入***************************************************************
L_2p = 29*1e-3; % 曲柄
L_3 = 37*1e-3; % 连杆
L_4 = 188*1e-3; % 摇杆
d = 189*1e-3; % 机架
w_2p0 = 100; % 曲柄的角速度
a_2p0 = 0; % 曲柄的角加速度
T = 2*pi/w_2p0;

I_2p = 0.005; % 转动惯量
I_4 = 0.025;
I_4p = 0.35;

C_2 = 1e5; % 扭转刚度 C_2 = 1e5  C_4 = 5e4
C_4 = 5e4;

k_2 = 7; % 阻尼系数  k_2 = 7  k_4 = 29 
k_4 = 29;

%% 黄皮书运动学及摇杆与曲柄的关系***************************************************************
YD = zeros(361,18);
phi_2p0 = acos(((L_2p+L_3)^2+d^2-L_4^2)/(2*(L_2p+L_3)*d)); % 1.3799 rad = 79.0624 °
% phi_2p0 = 0;

for i =1:361
    phi1_s = (i-1)*(pi/180)+phi_2p0; % 曲柄的转角，从前止心位置算起
    YD(i,1) = phi1_s;
    C_R = CRANK_ROCKER(L_2p,L_3,L_4,d,0,0,phi1_s,w_2p0,0,+1);
    YD(i,2)=C_R(1); % 摇杆的角位移 增量17.7499°
    YD(i,3)=C_R(2); % 摇杆的角速度
    YD(i,4)=C_R(3); % 摇杆的角加速度
    YD(i,5)=C_R(4); % 连杆的角位移
    YD(i,6)=C_R(5); % 连杆的角速度
    YD(i,7)=C_R(6); % 连杆的角加速度
    % 下面是u24和v24的两种计算方法: 验证正确
    u24 = C_R(2)/w_2p0;
    YD(i,8) = u24;
    v24=(C_R(3)-u24*a_2p0)/w_2p0;
    YD(i,9)=v24;
    % u_24=(L_2p*sin(phi1_s-YD(i,5)))/(L_4*sin(YD(i,2)-YD(i,5)));
    % YD(i,10)=u_24;
    % v_24=(L_3*YD(i,6)^2+(w_2p0^2)*L_2p*cos(phi1_s-YD(i,5))-(YD(i,3)^2)*L_4*cos(YD(i,2)-YD(i,5))) / (w_2p0*L_4*sin(YD(i,2)-YD(i,5)));
    % YD(i,11)=v_24;
    YD(i,18) = (i-1)*T/360; % 时间
end 

% 模拟外力矩M4 (经纱片对筘座的力矩)**************************************************************************
LJ = zeros(361,2);
M0 = 0; % 取值0或500
alpha = 0.26; % 自己取得，打纬脉冲周期
for i = 1:361
    x = YD(i,1); % 曲柄转角
    LJ(i,1) = x;
    if x>=phi_2p0 && x<=phi_2p0+alpha/2
        M4 = M0*sin(pi*(x-phi_2p0+alpha/2)/alpha);
    elseif x>=phi_2p0+alpha/2 && x <=(phi_2p0+2*pi-alpha/2)
        M4 = 0;
    else
        M4 = M0*sin(pi*(x-phi_2p0-2*pi+alpha/2)/alpha);
    end
    LJ(i,2)=M4;
end
% figure;
% plot(YD(:,18), LJ(:,2));
% title('构建外力矩M_4p');


%% 微分方程求解phi_4p和phi_4p`****************************************************************
tspan = [0, T];
x0 = [YD(1,2); YD(1,3)];


opts = odeset('RelTol',1e-4,'AbsTol',1e-7);
% 'RelTol'-相对误差容限 1e-3（默认） % 'AbsTol'-绝对误差容限 1e-6（默认）


[t, x] = ode45(@(t, x) sub_odefun(t, x, I_4p, C_4, k_4, L_2p, L_3, L_4, d, phi_2p0, w_2p0, M0, alpha), tspan, x0, opts);
alpha_4pp = (-C_4.*(x(:,1)-interp1(YD(:,18), YD(:,2), t(:,1))) - interp1(YD(:,18), LJ(:,2), t(:,1)) - k_4.*(x(:,2)-interp1(YD(:,18), YD(:,3), t(:,1)))) / I_4p;


phi_4p = interp1(t(:,1), x(:,1), YD(:,18)); 
omega_4p = interp1(t(:,1), x(:,2), YD(:,18)); 
alpha_4p = (-C_4.*(phi_4p-YD(:,2)) - LJ(:,2) - k_4.*(omega_4p-YD(:,3))) / I_4p;





% figure(2);
% plot(YD(:,18), YD(:,2));
% hold on;
% plot(t(:,1), x(:,1));
% legend('杆4的角位移', '插值前杆4p的角位移', 'Location', 'northeast');
% title('角位移对比');
% ylim([2.7, 3.2])

figure(3);
plot(YD(:,18), YD(:,2));
hold on;
plot(YD(:,18), phi_4p(:,1));
legend('杆4的角位移', '插值后杆4p的角位移', 'Location', 'northeast');
title('角位移对比');
ylim([2.7, 3.2])

figure(4);
plot(YD(:,18), YD(:,3));
hold on;
plot(YD(:,18), omega_4p(:,1));
legend('杆4的角速度', '插值后杆4p的角速度', 'Location', 'northeast');
title('角速度对比');
ylim([-30, 30])
% figure(5);
% plot(YD(:,18), YD(:,3));
% hold on;
% plot(t(:,1), x(:,2));
% legend('杆4的角速度', '插值前杆4p的角速度', 'Location', 'northeast');
% title('角速度对比');
% ylim([-30, 30])



figure(6);
plot(YD(:,18), YD(:,4));
hold on;
plot(YD(:,18), alpha_4p(:,1));
legend('杆4的角加速度', '插值后杆4p的角加速度', 'Location', 'northeast');
title('角加速度对比');
ylim([-4000, 5000])
txt = { ['杆4最大值', num2str(max(YD(:,4))), ' 杆4最小值', num2str(min(YD(:,4)))] , ['杆4P最大值', num2str(max(alpha_4p(:,1))), ' 杆4P最小值', num2str(min(alpha_4p(:,1)))]};
text(0.015, 2500, txt);
% figure(7);
% plot(YD(:,18), YD(:,4));
% hold on;
% plot(t(:,1), alpha_4pp(:,1));
% legend('杆4的角加速度', '插值前杆4p的角加速度', 'Location', 'northeast');
% title('角加速度对比');
% ylim([-3000, 4000])


