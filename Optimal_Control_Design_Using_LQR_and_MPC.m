%% 
% 
%% *1.LQR-parameters*

% פרמטרים
alpha = 0.5;              % קיבול תרמי יחסי
emf = 0.1;                % חום מנותב [kg/s]
num_stones = 15;          % מספר האבנים
Texh = 200;               % טמפרטורה חיצונית [°C]
Tint = 190;               % טמפרטורה התחלתית [°C]
Tdesired = 200;           % טמפרטורה רצויה (200 מעלות)

x0 = Tint * ones(num_stones, 1) - Tdesired;  % מצב התחלתי (סטייה מהמטרה)

A = zeros(num_stones, num_stones); 
B = zeros(num_stones, 1);

% בניית מטריצת A
for i = 1:num_stones
    if i == 1
        A(i, i) = -alpha * emf;          % אבן ראשונה
    else
        A(i, i) = -alpha * emf;          % שאר האבנים
        A(i, i-1) = alpha * emf;         % חיבור לאבן קודמת
    end
end

% בניית מטריצה B
B(1) = alpha * emf;

% בניית מטריצה Q
Q = 0* eye(num_stones);
Q(15,15) =1;

R = 0.01;

%% *LQR-unconstrained : שימוש רק בתנאי התחלה* 


% פתרון בעיית CARE
P = care(A, B, Q, R);  
K = inv(R) * B' * P;    % מטריצת הבקרה האופטימלית

% טווח זמן לסימולציה
tspan = 0:0.01:1000;        

% הגדרת המערכת באמצעות מצבי מצב (state-space)
sys = ss(A - B*K, B, eye(num_stones), zeros(num_stones,1));

% חישוב התגובה על פי תנאי התחלה x0
[ y,t,x] = initial(sys, x0, tspan);

% הוספת Tdesired ל-x כדי להחזיר לטמפרטורה האבסולוטית
x = x + Tdesired;

% צמצום על פי המגבלות של u(t) ו-T
T_max = 201.5;  % טמפרטורה מקסימלית

u_t1 = zeros(length(t), num_stones);  
x_fixed1 = x;

% ההגבלות הם רק על אות היציאה y 
for i = 1:length(t)
    u_t1(i,:) = -K * (x(i,:) - Tdesired)'+190;  

end

% גרפים
figure(1)
subplot(2,1,1)
plot(t, x_fixed1(:,15), 'LineWidth', 2);  % טמפרטורת האבן ה-15
title('Temperature of last stone(°C)');
xlabel('Time (s)');
ylabel('Temperature (°C)');
ylim([190 205]);  % שינוי טווח ציר ה-y
yline(200, '--r',  'LineWidth', 1.5);  % קו מקווקו אדום ב-200
yline(201.5, '--b', 'LineWidth', 1.5);  % קו מקווקו כחול ב-201.5
grid on;

subplot(2,1,2)
plot(t, u_t1(:,15), 'LineWidth', 2);  % מאמץ הבקרה לאורך זמן עבור האבן ה-15
title('Tin=deltaT+Texh');
xlabel('Time (s)');
ylabel('Tin(°C)');
grid on;
%% 
% 
%% *LQR-constrained*



% פתרון CARE וחישוב מטריצת K
P = care(A, B, Q, R);
K = inv(R) * B' * P;

% טווח זמן סימולציה
tspan = 0:0.1:1000;

% Constraints on control input
U_min = 0;      % Minimum control input
U_max = 80;     % Maximum control input

% Initialize arrays
u_t2 = zeros(length(tspan), 1);  % Control input over time (scalar)
x_fixed2 = zeros(length(tspan), num_stones);  % States over time

% Initial condition
x_fixed2(1, :) = x0';


for i = 1:length(tspan)-1
% u without constrain
    u_unconstrained = -K * x_fixed2(i, :)';  

%Apply control input constraints (scalar)
    u_t2(i) = min(max(u_unconstrained, U_min), U_max);
    
 % Simulate next state using constrained control input
    x_dot = A * x_fixed2(i, :)' + B * u_t2(i);  
    x_fixed2(i+1, :) = x_fixed2(i, :) + x_dot' * (tspan(2) - tspan(1));
   
end

% Final control input constraint
u_t2(end) = min(max(-K * x_fixed2(end, :)', U_min), U_max);
u_t2=u_t2+190;


% גרפים
figure(2)
subplot(2,1,1)
plot(tspan, x_fixed2(:, 15) + Tdesired, 'LineWidth', 2);  % Temperature of Stone 15
title('Temperature of last stone(°C)');
xlabel('Time (s)');
ylabel('Temperature (°C)');
ylim([180 210]);  % Adjust y-axis range
yline(200, '--r', 'LineWidth', 1.5);  % Target temperature
yline(201.5, '--b', 'LineWidth', 1.5);  % Maximum temperature constraint
grid on;

subplot(2,1,2)
plot(tspan, u_t2, 'LineWidth', 2);  % Control input over time
title('Tin=deltaT+Texh');
xlabel('Time (s)');
ylabel('Tin(°C)');
ylim([180 280]);
yline(190, '--r', 'LineWidth', 1.5);  % Minimum control limit
yline(270, '--b', 'LineWidth', 1.5);  % Maximum control limit
grid on;

%% *unconstrained VS.constrained*

figure(2)
subplot(2,1,1)
plot(t, u_t1(:,15), 'LineWidth', 2, 'Color', 'g');  % מאמץ הבקרה לאורך זמן
hold on;
plot(tspan, u_t2,  'LineWidth', 2, 'Color', 'm');  % Tin
title('Tin=deltaT+Texh');
xlabel('Time (s)');
ylabel('Tin(°C)');
legend('T (Catalyst Temperature)', 'Tin (Control Input)');
grid on;
yline(190, '--r',  'LineWidth', 1.5);  % קו מקווקו אדום ב-200
yline(270, '--b', 'LineWidth', 1.5);  % קו מקווקו כחול ב-201.5
legend('LQR idial', 'LQR real');  % הוספת מקרא עם השמות המתאימים
hold off

subplot(2,1,2)
plot(t, x_fixed1(:,end), 'LineWidth', 2, 'Color', 'g'); % צבע כחול עבור x_fixed1
hold on;
plot(tspan, x_fixed2(:,end)+Tdesired, 'LineWidth', 2, 'Color', 'm'); % צבע אדום עבור x_fixed2
title('Temperature of last stone(°C)');
xlabel('Time (s)');
ylabel('Temperature (°C)');
ylim([180, 210]);
yline(200, '--r',  'LineWidth', 1.5);  % קו מקווקו אדום ב-200
yline(201.5, '--b', 'LineWidth', 1.5);  % קו מקווקו כחול ב-201.5
grid on;
legend('LQR idial', 'LQR real');  % הוספת מקרא עם השמות המתאימים

%% *MPC* 

% פרמטרים
alpha = 0.5;              % קיבול תרמי יחסי
emf = 0.1;                % חום מנותב [kg/s]
num_stones = 15;          % מספר האבנים
Texh = 200;               % טמפרטורה חיצונית [°C]
Tint = 190;               % טמפרטורה התחלתית [°C]
Tdesired = 200; 
num_stones=15;
% הגדרות מטריצות
C = eye(15); % Full state feedback
D = zeros(15,1);



A = zeros(num_stones, num_stones); 
B = zeros(num_stones, 1);

% בניית מטריצת A
for i = 1:num_stones
    if i == 1
        A(i, i) = -alpha * emf;          % אבן ראשונה
    else
        A(i, i) = -alpha * emf;          % שאר האבנים
        A(i, i-1) = alpha * emf;         % חיבור לאבן קודמת
    end
end

% בניית מטריצה B
B(1) = alpha * emf;
Ts = 5; 
sys_c = ss(A, B, C, D); % מערכת זמן רציף
sys_d = c2d(sys_c, Ts, 'zoh'); % המרה לדיסקרטית עם שימור אפס-החזקה (ZOH)

Ad = sys_d.A; % מטריצה A בדידה
Bd = sys_d.B; % מטריצה B בדידה


% הגבלות על אות בקרה 
umin = 0;
umax = 80;


x0 = Tint * ones(num_stones, 1) - Tdesired;  % מצב התחלתי (סטייה מהמטרה)
Tstop = 1000;      % Simulation time

% חישוב מטריצת LQR
[K,Qp,~] = dlqr(Ad,Bd,Q,R, zeros(15,1));  % Alternative:  [K,Qp] = lqry(Plant,Q,R);
L = chol(Qp);
Plant = ss(Ad, Bd, [C; L], [D; zeros(size(L, 1), size(Bd, 2))], Ts);


%% Defining signal groups
Plant.InputGroup.MV = 1;           % Manipulated vars
Plant.OutputGroup.MO = 1:15;         % Measured outputs
Plant.OutputGroup.UO = 16:30;      % Unmeasured outputs

% הגדרת MPC
mpcObj = mpc(Plant, Ts);
mpcObj.PredictionHorizon = 20;
mpcObj.ControlHorizon = 20;

% הגדרת משקלים עבור ה-MPC
yWeights = sqrt(diag(Q))';
uWeights = sqrt(diag(R))';
mpcObj.Weights.OV = [yWeights zeros(1,15)];
mpcObj.Weights.MV = uWeights;

% הגבלות על אות בקרה
mpcObj.MV.Min = umin;
mpcObj.MV.Max = umax;


%% Terminal conditions
Y = struct('Weight',[zeros(1,15) ones(1,15)]);
U = struct('Weight',uWeights);
setterminal(mpcObj,Y,U);


% הגדרת Model Disturbance (מודל להפרעות ביציאה)
setoutdist(mpcObj, 'model', ss(zeros(2*num_stones,1)));  % הגדרת מודל להפרעות עבור 15 יציאות
setEstimator(mpcObj, [], eye(num_stones));  % הערכת מצב על פי מערכת עם 15 משתנים



% יצירת וקטור זמן והכנה לעיבוד התוצאות
Nt = ceil(Tstop/Ts);
t = NaN(1, Nt);
U = NaN(1, Nt);
X = NaN(num_stones, Nt);

% אתחול
curt = 0;
curX = x0;

% לולאת סימולציה עיקרית
for ind = 1:Nt
    
    % חישוב אות הבקרה עבור MPC
    xmpc = mpcstate(mpcObj);
    xmpc.Plant = curX;
    ym_cur = curX;
    r_cur = zeros(30, 1);
    [u, ~] = mpcmove(mpcObj, xmpc, ym_cur, r_cur);
    curU = u(1);
    
    % רישום תוצאות
    t(ind) = curt;
    U(ind) = curU;
    X(:,ind) = curX;

    % עדכון הדינמיקה של המודל
    curt = curt + Ts;
    curX = Ad*curX + Bd*curU;
end

% הוספת ערך 200 לתוצאות X
X = X + 200;

% גרפים 
figure(3)
subplot(2,1,1)
h1 = plot(t, X(15,:), 'b-', 'LineWidth', 1.5);  % MPC controller
hold on
h2 = plot(tspan, x_fixed2(:, 15) + Tdesired, 'LineWidth', 2);  % Constrained LQR control input
yline(200, '--r', 'LineWidth', 1.5);  % קו מקווקו אדום ב-200
yline(201.5, '--b', 'LineWidth', 1.5);  % קו מקווקו כחול ב-201.5
title('MPC vs LQR with the same weights');
legend([h1, h2], {'MPC controller', 'Constrained LQR'}, 'Location', 'Best');
xlabel('Time (s)');
ylabel('Plant Outputs');
grid on
hold off


subplot(2,1,2)
U = U + 190;
h1 = plot(t, U, 'r-', 'LineWidth', 1.5); % MPC control input
hold on
h2 = plot(tspan, u_t2, 'LineWidth', 2);  % Constrained LQR control input
grid on
xlabel('Time (s)');
ylabel('Tin(°C)');
ylim([180 280]);
yline(190, '--r', 'LineWidth', 1.5);  % Minimum control limit
yline(270, '--b', 'LineWidth', 1.5);  % Maximum control limit
legend([h1, h2], {'MPC controller', 'Constrained LQR'}, 'Location', 'Best');