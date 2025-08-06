clear
%==================limit cycle=================
lags = [1 2];
tspan = linspace(0,80,100);
y0 = [600,1000];
sol = dde23(@pn,lags,@(t)yhist(t,y0),tspan);

t=sol.x;
y=sol.y;

figure('Color', 'k'); % Black background for contrast
hold on

% Limit cycle in white
p1 = plot(y(1,60:end),y(2,60:end), 'w', 'LineWidth', 2.5);

box on

%=================nullcline=====================
max_x = 1100;
max_y = 1200;
x_lim = [0,max_x];
y_lim = [0,max_y];
[H_null,S_null]=nullcline(x_lim,y_lim);

% H-nullcline: blue-cyan
p2 = plot(H_null(:,1),H_null(:,2), '-', 'LineWidth', 2.5, 'Color', [0 0.7 1]);
% S-nullcline: orange-red
p3 = plot(S_null(:,1),S_null(:,2), '-', 'LineWidth', 2.5, 'Color', [1 0.3 0.2]);

%===================quiver======================
xList=linspace(1,max_x,15);
yList=linspace(1,max_y,15);
[X, Y] = meshgrid(xList, yList);
ts = 1:2;
LineLength=0.04; % Slightly shorter arrows for elegance

for i=1:numel(X)
    y00=[X(i),Y(i)];
    sol = dde23(@pn,lags,@(t)yhist(t,y00),ts);
    y=sol.y;
    U=y(1,end)-y(1,1);
    V=y(2,end)-y(2,1);
    quiver(X(i),Y(i),U*LineLength,V*LineLength,0, ...
        'Color', [0.7 0.7 0.7], 'LineWidth', 1, 'MaxHeadSize', 2, 'AutoScale', 'off');
end

ylim([0,max_y]);
xlim([0,max_x]);
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'FontSize', 14, 'FontName', 'Arial')
xlabel('H', 'Color', 'w', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('S', 'Color', 'w', 'FontSize', 16, 'FontWeight', 'bold')
title('DDE Phase Portrait', 'Color', 'w', 'FontSize', 18, 'FontWeight', 'bold')
legend([p1 p2 p3], {'Limit Cycle', 'H Nullcline', 'S Nullcline'}, ...
    'TextColor', 'w', 'FontSize', 12, 'Location', 'northeast', 'Box', 'off')

% =================time traj=====================
axes('Position',[.68 .7 .23 .23])
tspan = linspace(0,80,100);
y0 = [600,1000];
sol = dde23(@pn,lags,@(t)yhist(t,y0),tspan);
y=sol.y;
plot(tspan, y(2,:), 'Color', [1 0.5 0.5], 'LineWidth', 2)
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'FontSize', 10, 'FontName', 'Arial')
title('S(t)', 'Color', 'w', 'FontSize', 12)
box on

%=================Functions=====================
function dydt = pn(t,y,Z)
    ylag1 = Z(:,1);
    ylag2 = Z(:,2);
    as = 30.5;
    ah = 183;
    ah0 = 0.1;
    as0 = 0.1;
    beta = 4.6;
    dm = 0.3;
    dh = 3.8;
    ds = 0.2;
    Kh = 326;
    Ks = 185;
    n1 = 3;
    n2 = 4.8;
    K = beta/dm;
    dydt = [K*ah0 + K*ah*Ks^n2/(Ks^n2+ylag2(2)^n2)-dh*y(1);
            K*as0 + K*as*ylag1(1)^n1/(Kh^n1+ylag1(1)^n1)-ds*y(2)];
end

function s=yhist(~,y0)
    s = y0;
end

function [H_null,S_null]=nullcline(xlim,ylim)
    steps = 2000;
    x_list = linspace(xlim(1),xlim(2),steps+1);
    y_list = linspace(ylim(1),ylim(2),steps+1);
    m=1; n=1;
    for i=2:length(x_list)
        xitem_1 = x_list(i-1);
        xitem_2 = x_list(i);
        for j=2:length(y_list)
            yitem_1 = y_list(j-1);
            [x1,y1] = quiver_(xitem_1,yitem_1);
            [x2,y2] = quiver_(xitem_2,yitem_1);
            if sign(x1)~=sign(x2)
                H_null(m,1) = xitem_1; 
                H_null(m,2) = yitem_1; 
                m=m+1;
            elseif sign(y1)~=sign(y2)
                S_null(n,1) = xitem_1; 
                S_null(n,2) = yitem_1;
                n=n+1;
            end
        end
    end
end

function [u,v]=quiver_(H,S)
    as = 30.5;
    ah = 183;
    ah0 = 0.1;
    as0 = 0.1;
    beta = 4.6;
    dm = 0.3;
    dh = 3.8;
    ds = 0.2;
    Kh = 326;
    Ks = 185;
    n1 = 3;
    n2 = 4.8;
    K = beta/dm;
    u = K*ah0 + K*ah*Ks.^n2./(Ks^n2+S.^n2)-dh*H;
    v =  K*as0+K*as*H.^n1./(Kh^n1+H.^n1)-ds*S;
end
