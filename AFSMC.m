clc; clear; close all;

r = 0.1; b = 0.8; d = 0.5;
m = 300; IA = 20; IB = 20;
k1 = 7; k2 = 7; k3 = 5; rho = 4; sigma = diag([2 2]);

X0 = [2.2; 1.7; atan2(2.2,1.7); 0; 0; 0; 0; 0; 0; 0];
tspan = 0:0.001:10;

fis = mamfis('Name','FuzzyKS');
fis = addInput(fis, [-2 2], 'Name','ev');
fis = addMF(fis, 'ev', 'trapmf', [-2 -2 -1.5 -1], 'Name','NB');
fis = addMF(fis, 'ev', 'trimf', [-1.5 -1 -0.5], 'Name','NM');
fis = addMF(fis, 'ev', 'trimf', [-1 -0.5 0], 'Name','NS');
fis = addMF(fis, 'ev', 'trimf', [-0.5 0 0.5], 'Name','ZO');
fis = addMF(fis, 'ev', 'trimf', [0 0.5 1], 'Name','PS');
fis = addMF(fis, 'ev', 'trimf', [0.5 1 1.5], 'Name','PM');
fis = addMF(fis, 'ev', 'trapmf', [1 1.5 2 2], 'Name','PB');
fis = addOutput(fis, [0 0.4], 'Name','ks');
fis = addMF(fis, 'ks', 'trimf', [0 0 0.1], 'Name','ZO');
fis = addMF(fis, 'ks', 'trimf', [0 0.1 0.2], 'Name','PS');
fis = addMF(fis, 'ks', 'trimf', [0.1 0.2 0.3], 'Name','PM');
fis = addMF(fis, 'ks', 'trimf', [0.2 0.3 0.4], 'Name','PB');
rules = [
    "ev==NB => ks=PB";
    "ev==NM => ks=PM";
    "ev==NS => ks=PS";
    "ev==ZO => ks=ZO";
    "ev==PS => ks=PS";
    "ev==PM => ks=PM";
    "ev==PB => ks=PB"
];
fis = addRule(fis, rules);

N = length(tspan);
tau1_arr = zeros(N,1);
tau2_arr = zeros(N,1);
vd_arr = zeros(N,1);
wd_arr = zeros(N,1);
e_c = zeros(N,1);
ex_arr = zeros(N, 1);
ey_arr = zeros(N, 1);
etheta_arr = zeros(N, 1);
vc_prev = [0; 0];
vcdot_prev = [0; 0];
s_arr  = zeros(N,2);
ev_arr = zeros(N,2);
int_ev = zeros(1,2);
ev_prev=[0,0];
K_mat = [1/r, -b/(2*r);
         1/r,  b/(2*r)];

X = X0.';
for i = 2:N
    t = tspan(i);
    if t < 4
        mc = 100; I = 200;
    else
        mc = 200; I = 250;
    end
    q = X(i-1,1:5).';
    dq = X(i-1,6:10).';
    x=q(1); y=q(2); theta = q(3); phi1=q(4); phi2=q(5);
    vx=dq(1); vy=dq(2); w=dq(3); dphi1=dq(4); dphi2=dq(5);

    ref_x = 2*cos(t);
    ref_y = 2*sin(t);
    ref_x_dot = -2*sin(t);
    ref_y_dot =  2*cos(t);
    ref_theta = atan2(ref_y_dot, ref_x_dot);
    % ref_theta=deg2rad(t);
    vr = 2; wr = 1;
    qd=[q(1); q(2); q(3)];
    eq = [ref_x - q(1); ref_y - q(2); wrapToPi(ref_theta - q(3))];

    R_theta = [cos(theta), sin(theta), 0; -sin(theta), cos(theta), 0; 0, 0, 1];
    ec = R_theta * eq;

    ex_arr(i) = ec(1);
    ey_arr(i) = ec(2);
    etheta_arr(i) = ec(3);

    vd = k1*ec(1) + vr*cos(ec(3));
    wd = wr + k2*vr*ec(2) + k3*vr*sin(ec(3));
    vc = [vd; wd];

    v=dq(4:5);

    v_left = dq(4);
    v_right = dq(5);
    vd_actual = (r/2) * (v_left + v_right);
    wd_actual = (r/b) * (v_right - v_left);
    VD= [vd_actual; wd_actual];

    ev = vc - VD;

    J = [r/2*cos(q(3)), r/2*cos(q(3));
         r/2*sin(q(3)), r/2*sin(q(3));
         -r/b, r/b;
         1, 0;
         0, 1];

    dJ = [ -r/2*sin(q(3))*dq(3), -r/2*sin(q(3))*dq(3);
            r/2*cos(q(3))*dq(3),  r/2*cos(q(3))*dq(3);
            0, 0;
            0, 0;
            0, 0 ];

    M = zeros(5,5);
    m_t = m + mc;
    M(1,1) = m_t;
    M(2,2) = m_t;
    M(3,3) = I+(m_t*d^2);
    M(4,4) = IA;
    M(5,5) = IB;
    M(1,3) = -m_t * d * sin(theta);
    M(3,1) = -m_t * d * sin(theta);
    M(2,3) = m_t * d * cos(theta);
    M(3,2) = m_t * d * cos(theta);

    Vm = zeros(5,5);
    Vm(1,3) = -m_t * d * dq(3)^2 * cos(theta);
    Vm(2,3) = -m_t * d * dq(3)^2 * sin(theta);
    Vm(3,1) = m_t * d * dq(3)^2 * cos(theta);
    Vm(3,2) = m_t * d * dq(3)^2 * sin(theta);

    E = [zeros(3,2); eye(2)];

    tau_d = [0; 0; 0; 4*sin(t); 4*sin(t)];

    if i > 1
        dt = tspan(i) - tspan(i-1);
        int_ev = int_ev + (ev' + ev_prev) * dt / 2;
    end
    ev_prev = ev';
    s = ev + (rho * int_ev');
    s_arr(i,:)  = s.';
    ev_arr(i,:) = ev.';
    ev_input = ev(1);
    ev_input = min(max(ev_input, -2), 2);
    ks = evalfis(fis, ev_input);
     % ks=1;
     
    H = J' * M * J * K_mat;
    N_mat = J' * ((M * dJ) + (Vm * J ))* K_mat;
    P_mat = J' * E;

    if i == 2
        vcdot = [0; 0];
    else
        vcdot = (vc - vc_prev) / (dt);
    end
    vc_prev = vc;

    tau_eq =  inv(P_mat)*(H * vcdot) + inv(P_mat)*(N_mat*VD) - inv(P_mat)*J'*tau_d + (rho *inv(P_mat)* H * ev);
    tau_sw = inv(P_mat) * H * sigma * sign(s);
    tau = tau_eq + (ks * tau_sw);

    tau1_arr(i) = tau(1);
    tau2_arr(i) = tau(2);
    vd_arr(i) = vd_actual;
    wd_arr(i) = wd_actual;
    V=K_mat*vc;
    v_dot = inv(J'*M*J)*J'*((-M*dJ*V-Vm*dq)+(E*tau+tau_d));

    dq_new = dq;
    dq_new(4:5) = dq(4:5) + v_dot * dt;

    x     = q(1); y = q(2); theta = q(3);
    phi1  = q(4); phi2 = q(5);

    theta_new = theta + r/b * (dq_new(5) - dq_new(4)) * dt;
    x_new     = x + r/2 * cos(theta_new) * (dq_new(4) + dq_new(5)) * dt;
    y_new     = y + r/2 * sin(theta_new) * (dq_new(4) + dq_new(5)) * dt;
    phi1_new  = phi1 + dq_new(4) * dt;
    phi2_new  = phi2 + dq_new(5) * dt;
    omega1 = dq_new(4);
    omega2 = dq_new(5);

    v_new = r/2 * (omega1 + omega2);
    w = r/b * (omega2 - omega1);

    dq_new(1) = v_new * cos(theta_new);
    dq_new(2) = v_new * sin(theta_new);
    dq_new(3) = w;
    q_new = [x_new; y_new; theta_new; phi1_new; phi2_new];
    X(i,:) = [q_new; dq_new]';
end

figure; plot(X(:,1), X(:,2), 'b', 'LineWidth',1.5); hold on
plot(2*cos(tspan), 2*sin(tspan),'r--','LineWidth',1.2)
legend('Robot Trajectory', 'Reference Trajectory')
xlabel('x (m)'); ylabel('y (m)'); axis equal
title('Trajectory Tracking');

figure; plot(tspan, vd_arr, 'b','LineWidth',1.5)
xlabel('Time (s)'); ylabel('v_d (m/s)'); title('Desired Linear Velocity v_d');

figure; plot(tspan, wd_arr, 'r','LineWidth',1.5)
xlabel('Time (s)'); ylabel('w_d (rad/s)'); title('Desired Angular Velocity w_d');

figure; plot(tspan, tau1_arr, 'b', tspan, tau2_arr, 'r','LineWidth',1.2)
legend('\tau(Left)','\tau(Right)'); xlabel('Time (s)'); ylabel('Torque (Nm)');
title('Wheel Torques');

figure;
plot(tspan, ex_arr, 'b', 'LineWidth', 1.1, 'DisplayName', 'e_x (m)');
hold on;
plot(tspan, ey_arr, 'r', 'LineWidth', 1.1, 'DisplayName', 'e_y (m)');
plot(tspan, etheta_arr, 'k', 'LineWidth', 1.1, 'DisplayName', 'e_\theta (rad)');
hold off;
xlabel('Time (s)');
ylabel('Tracking Errors');
title('Tracking Errors (e_x, e_y, e_\theta)');
legend('show');
grid on;

saveVideo = true;
videoFile = 'robot_animation.mp4';
fps = 60;
dt_mean = mean(diff(tspan));
subsample = max(1, round(1/(dt_mean*fps)));
idx = 1:subsample:length(tspan);
refx = 2*cos(tspan);
refy = 2*sin(tspan);
base = b;
bodyL = 1.2*base;
bodyW = 0.8*base;
wheelL = 0.35*base;
wheelW = 0.12*base;
headLen = 0.5*bodyL;
xmin = min([X(:,1); refx(:)]) - 1.5*bodyL;
xmax = max([X(:,1); refx(:)]) + 1.5*bodyL;
ymin = min([X(:,2); refy(:)]) - 1.5*bodyL;
ymax = max([X(:,2); refy(:)]) + 1.5*bodyL;

fig = figure('Color','w');
axes('Position',[0.06 0.08 0.88 0.88]); hold on; axis equal;
plot(refx, refy, 'r--', 'LineWidth', 1.2, 'DisplayName','Reference');
trail = animatedline('LineWidth', 2.0, 'DisplayName','Path');
xlim([xmin xmax]); ylim([ymin ymax]);
xlabel('x (m)'); ylabel('y (m)');
title('Robot Motion Animation');
legend('Location','best'); grid on;

x = X(1,1); y = X(1,2); th = X(1,3);
[bx, by] = rectCorners(x, y, th, bodyL, bodyW, 0, 0);
[wlx, wly] = rectCorners(x, y, th, wheelL, wheelW, +0.0, +base/2);
[wrx, wry] = rectCorners(x, y, th, wheelL, wheelW, +0.0, -base/2);
bodyPatch  = patch(bx,  by,  [0.85 0.85 0.88], 'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 1.2);
leftPatch  = patch(wlx, wly, [0.3 0.3 0.3],  'EdgeColor', 'k');
rightPatch = patch(wrx, wry, [0.3 0.3 0.3],  'EdgeColor', 'k');
headArrow  = quiver(x, y, headLen*cos(th), headLen*sin(th), 0, 'LineWidth', 1.6, 'MaxHeadSize', 0.8);

if saveVideo
    vw = VideoWriter(videoFile, 'MPEG-4');
    vw.FrameRate = fps;
    open(vw);
end

for k = idx
    x = X(k,1);  y = X(k,2);  th = X(k,3);
    [bx, by] = rectCorners(x, y, th, bodyL, bodyW, 0, 0);
    set(bodyPatch,  'XData', bx,  'YData', by);
    [wlx, wly] = rectCorners(x, y, th, wheelL, wheelW, +0.1*wheelL, +base/2);
    [wrx, wry] = rectCorners(x, y, th, wheelL, wheelW, +0.1*wheelL, -base/2);
    set(leftPatch,  'XData', wlx, 'YData', wly);
    set(rightPatch, 'XData', wrx, 'YData', wry);
    set(headArrow, 'XData', x, 'YData', y, 'UData', headLen*cos(th), 'VData', headLen*sin(th));
    addpoints(trail, x, y);
    drawnow limitrate;
    if saveVideo
        frame = getframe(fig);
        writeVideo(vw, frame);
    end
    pause(0.02)
end

if saveVideo
    close(vw);
    fprintf('Saved video -> %s\n', videoFile);
end

figure; hold on; plot(ex_arr,ey_arr,'b','LineWidth',1.2); plot(ex_arr(2),ey_arr(2),'go','MarkerFaceColor','g','MarkerSize',6); plot(ex_arr(end),ey_arr(end),'rx','MarkerSize',8,'LineWidth',1.2); xlabel('e_x'); ylabel('e_y'); title('Phase Portrait (Baseline)'); grid on;
figure; hold on; plot(s_arr(:,1),s_arr(:,2),'b','LineWidth',1.2); plot(s_arr(2,1),s_arr(2,2),'go','MarkerFaceColor','g','MarkerSize',6); plot(s_arr(end,1),s_arr(end,2),'rx','MarkerSize',8,'LineWidth',1.2); xlabel('s_v'); ylabel('s_\omega'); title('Sliding Surface Phase Plot (Baseline)'); grid on;
figure; subplot(2,1,1); plot(tspan,s_arr(:,1),'b','LineWidth',1.1); hold on; plot(tspan(2),s_arr(2,1),'go','MarkerFaceColor','g','MarkerSize',6); plot(tspan(end),s_arr(end,1),'rx','MarkerSize',8,'LineWidth',1.2); xlabel('t (s)'); ylabel('s_v'); title('Sliding Surface Components (Baseline)'); subplot(2,1,2); plot(tspan,s_arr(:,2),'b','LineWidth',1.1); hold on; plot(tspan(2),s_arr(2,2),'go','MarkerFaceColor','g','MarkerSize',6); plot(tspan(end),s_arr(end,2),'rx','MarkerSize',8,'LineWidth',1.2); xlabel('t (s)'); ylabel('s_\omega');

XA = zeros(N,10); XA(1,:) = X0.';
tau1_A = zeros(N,1); tau2_A = zeros(N,1);
vd_A = zeros(N,1); wd_A = zeros(N,1);
ex_A = zeros(N,1); ey_A = zeros(N,1); eth_A = zeros(N,1);
sA_arr = zeros(N,2); evA_arr = zeros(N,2);
c_hat = [1;1]; gamma = [5;5]; eta = [1;1]; cmin = [0.1;0.1]; cmax = [50;50]; phi = [0.05;0.05];
satfun = @(z,ph) max(-1,min(1,z./ph));
vc_prev_A = [0;0]; int_ev_A = zeros(1,2); ev_prev_A = [0,0];
c1_hist = zeros(N,1); c2_hist = zeros(N,1); c1_hist(1)=c_hat(1); c2_hist(1)=c_hat(2);

for i = 2:N
    t = tspan(i);
    if t < 4
        mc = 100; I = 200;
    else
        mc = 200; I = 250;
    end
    q = XA(i-1,1:5).'; dq = XA(i-1,6:10).';
    x = q(1); y = q(2); theta = q(3); v_left = dq(4); v_right = dq(5);

    ref_x = 2*cos(t); ref_y = 2*sin(t);
    ref_x_dot = -2*sin(t); ref_y_dot = 2*cos(t);
    ref_theta = atan2(ref_y_dot, ref_x_dot);
    vr = 2; wr = 1;
    eq = [ref_x - q(1); ref_y - q(2); wrapToPi(ref_theta - q(3))];
    R_theta = [cos(theta), sin(theta), 0; -sin(theta), cos(theta), 0; 0, 0, 1];
    ec = R_theta*eq;
    ex_A(i) = ec(1); ey_A(i) = ec(2); eth_A(i) = ec(3);

    vd = k1*ec(1) + vr*cos(ec(3));
    wd = wr + k2*vr*ec(2) + k3*vr*sin(ec(3));
    vc = [vd; wd];

    vd_actual = (r/2)*(v_left + v_right);
    wd_actual = (r/b)*(v_right - v_left);
    VD = [vd_actual; wd_actual];

    ev = vc - VD;
    if i > 1
        dt = tspan(i) - tspan(i-1);
        int_ev_A = int_ev_A + (ev' + ev_prev_A)*dt/2;
    end
    ev_prev_A = ev';
    s = ev + (rho*int_ev_A');
    sA_arr(i,:) = s.'; evA_arr(i,:) = ev.';

    J = [r/2*cos(q(3)), r/2*cos(q(3));
         r/2*sin(q(3)), r/2*sin(q(3));
         -r/b, r/b;
         1, 0;
         0, 1];
    dJ = [ -r/2*sin(q(3))*dq(3), -r/2*sin(q(3))*dq(3);
            r/2*cos(q(3))*dq(3),  r/2*cos(q(3))*dq(3);
            0, 0;
            0, 0;
            0, 0 ];

    M = zeros(5,5); m_t = m + mc;
    M(1,1)=m_t; M(2,2)=m_t; M(3,3)=I+(m_t*d^2); M(4,4)=IA; M(5,5)=IB;
    M(1,3)=-m_t*d*sin(theta); M(3,1)=M(1,3);
    M(2,3)= m_t*d*cos(theta); M(3,2)=M(2,3);

    Vm = zeros(5,5);
    Vm(1,3)=-m_t*d*dq(3)^2*cos(theta);
    Vm(2,3)=-m_t*d*dq(3)^2*sin(theta);
    Vm(3,1)= m_t*d*dq(3)^2*cos(theta);
    Vm(3,2)= m_t*d*dq(3)^2*sin(theta);

    E = [zeros(3,2); eye(2)];
    tau_d = [0;0;0;4*sin(t);4*sin(t)];

    if i == 2
        vcdot = [0;0];
    else
        vcdot = (vc - vc_prev_A)/dt;
    end
    vc_prev_A = vc;

    H = J'*M*J*K_mat;
    N_mat = J'*((M*dJ) + (Vm*J))*K_mat;
    P_mat = J'*E;

    tau_eq = inv(P_mat)*(H*vcdot) + inv(P_mat)*(N_mat*VD) - inv(P_mat)*J'*tau_d + (rho*inv(P_mat)*H*ev);
    c_hat = c_hat + dt.*(gamma.*abs(s) - eta.*c_hat);
    c_hat = min(cmax, max(cmin, c_hat));
    tau_sw_adapt = inv(P_mat)*H*(diag(c_hat)*sigma*satfun(s,phi));
    tauA = tau_eq + tau_sw_adapt;

    tau1_A(i) = tauA(1); tau2_A(i) = tauA(2);
    vd_A(i) = vd_actual; wd_A(i) = wd_actual;

    Vc = K_mat*vc;
    v_dotA = inv(J'*M*J)*J'*(-M*dJ*Vc - Vm*dq + (E*tauA + tau_d));

    dq_new = dq; dq_new(4:5) = dq(4:5) + v_dotA*dt;
    theta_new = theta + r/b*(dq_new(5)-dq_new(4))*dt;
    x_new = x + r/2*cos(theta_new)*(dq_new(4)+dq_new(5))*dt;
    y_new = y + r/2*sin(theta_new)*(dq_new(4)+dq_new(5))*dt;
    phi1_new = q(4) + dq_new(4)*dt; phi2_new = q(5) + dq_new(5)*dt;
    omega1 = dq_new(4); omega2 = dq_new(5);
    v_new = r/2*(omega1+omega2); w_new = r/b*(omega2-omega1);
    dq_new(1) = v_new*cos(theta_new);
    dq_new(2) = v_new*sin(theta_new);
    dq_new(3) = w_new;
    q_new = [x_new; y_new; theta_new; phi1_new; phi2_new];
    XA(i,:) = [q_new; dq_new]';
    c1_hist(i)=c_hat(1); c2_hist(i)=c_hat(2);
end

figure; hold on; plot(X(:,1),X(:,2),'b','LineWidth',1.5); plot(XA(:,1),XA(:,2),'m','LineWidth',1.5); plot(2*cos(tspan),2*sin(tspan),'r--','LineWidth',1.2); axis equal; legend('Baseline','Adaptive','Reference'); xlabel('x (m)'); ylabel('y (m)'); title('Trajectory: Baseline vs Adaptive');

figure; subplot(3,1,1); plot(tspan,ex_arr,'b',tspan,ex_A,'m','LineWidth',1.1); xlabel('t (s)'); ylabel('e_x'); legend('Base','Adapt'); title('Errors'); subplot(3,1,2); plot(tspan,ey_arr,'b',tspan,ey_A,'m','LineWidth',1.1); xlabel('t (s)'); ylabel('e_y'); legend('Base','Adapt'); subplot(3,1,3); plot(tspan,etheta_arr,'b',tspan,eth_A,'m','LineWidth',1.1); xlabel('t (s)'); ylabel('e_\theta'); legend('Base','Adapt');

figure; plot(tspan,tau1_arr,'b',tspan,tau2_arr,'r','LineWidth',1.0); hold on; plot(tspan,tau1_A,'m--',tspan,tau2_A,'k--','LineWidth',1.0); legend('\tau_L Base','\tau_R Base','\tau_L Adapt','\tau_R Adapt'); xlabel('t (s)'); ylabel('Torque (Nm)'); title('Wheel Torques: Baseline vs Adaptive');

figure; hold on; plot(ex_arr,ey_arr,'b','LineWidth',1.2); plot(ex_A,ey_A,'m','LineWidth',1.2); plot(ex_arr(2),ey_arr(2),'go','MarkerFaceColor','g','MarkerSize',6); plot(ex_arr(end),ey_arr(end),'rx','MarkerSize',8,'LineWidth',1.2); plot(ex_A(2),ey_A(2),'mo','MarkerFaceColor','m','MarkerSize',6); plot(ex_A(end),ey_A(end),'kx','MarkerSize',8,'LineWidth',1.2); xlabel('e_x'); ylabel('e_y'); legend('Base','Adapt','Base start','Base end','Adapt start','Adapt end'); title('Phase Portrait (e_x vs e_y)'); grid on;

figure; hold on; plot(s_arr(:,1),s_arr(:,2),'b','LineWidth',1.2); plot(sA_arr(:,1),sA_arr(:,2),'m','LineWidth',1.2); plot(s_arr(2,1),s_arr(2,2),'go','MarkerFaceColor','g','MarkerSize',6); plot(s_arr(end,1),s_arr(end,2),'rx','MarkerSize',8,'LineWidth',1.2); plot(sA_arr(2,1),sA_arr(2,2),'mo','MarkerFaceColor','m','MarkerSize',6); plot(sA_arr(end,1),sA_arr(end,2),'kx','MarkerSize',8,'LineWidth',1.2); xlabel('s_v'); ylabel('s_\omega'); legend('Base','Adapt','Base start','Base end','Adapt start','Adapt end'); title('Sliding Surface Phase Plot'); grid on;

figure; subplot(2,1,1); plot(tspan,s_arr(:,1),'b','LineWidth',1.1); hold on; plot(tspan(2),s_arr(2,1),'go','MarkerFaceColor','g','MarkerSize',6); plot(tspan(end),s_arr(end,1),'rx','MarkerSize',8,'LineWidth',1.2); plot(tspan,sA_arr(:,1),'m','LineWidth',1.1); plot(tspan(2),sA_arr(2,1),'mo','MarkerFaceColor','m','MarkerSize',6); plot(tspan(end),sA_arr(end,1),'kx','MarkerSize',8,'LineWidth',1.2); xlabel('t (s)'); ylabel('s_v'); legend('Base','Base start','Base end','Adapt','Adapt start','Adapt end'); title('Sliding Surface Components'); subplot(2,1,2); plot(tspan,s_arr(:,2),'b','LineWidth',1.1); hold on; plot(tspan(2),s_arr(2,2),'go','MarkerFaceColor','g','MarkerSize',6); plot(tspan(end),s_arr(end,2),'rx','MarkerSize',8,'LineWidth',1.2); plot(tspan,sA_arr(:,2),'m','LineWidth',1.1); plot(tspan(2),sA_arr(2,2),'mo','MarkerFaceColor','m','MarkerSize',6); plot(tspan(end),sA_arr(end,2),'kx','MarkerSize',8,'LineWidth',1.2); xlabel('t (s)'); ylabel('s_\omega');

figure; plot(tspan, vecnorm(sA_arr,2,2),'m','LineWidth',1.3); xlabel('t (s)'); ylabel('||s||_2'); title('Adaptive: Sliding Surface Norm');

figure; plot(tspan,c1_hist,'m','LineWidth',1.2); hold on; plot(tspan,c2_hist,'k','LineWidth',1.2); xlabel('t (s)'); ylabel('c_i'); legend('c_1','c_2'); title('Adaptive Gains');

function [Xc, Yc] = rectCorners(cx, cy, theta, L, W, ox, oy)
    xh = L/2; yh = W/2;
    rectB = [ -xh, -yh;
               xh, -yh;
               xh,  yh;
              -xh,  yh ];
    rectB = rectB + [ox, oy];
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    rectW = (R * rectB.').';
    rectW(:,1) = rectW(:,1) + cx;
    rectW(:,2) = rectW(:,2) + cy;
    Xc = [rectW(:,1); rectW(1,1)];
    Yc = [rectW(:,2); rectW(1,2)];
end
