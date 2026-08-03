clc; clear; close all;

r = 0.1; b = 0.8; d = 0.5;
m = 300; IA = 20; IB = 20;
k1 = 7; k2 = 7; k3 = 5; 

Kp_BS = diag([10 11]);   
Ki_BS = diag([0.5 0.3]);
Kd_BS = diag([0 0]);
u_sat = 7000;

X0 = [2.2; 1.7; atan2(2.2,1.7); 0; 0; 0; 0; 0; 0; 0];
tspan = 0:0.001:10;

N = length(tspan);
tau1_arr = zeros(N,1);
tau2_arr = zeros(N,1);
vd_arr = zeros(N,1);
wd_arr = zeros(N,1);
ex_arr = zeros(N,1); 
ey_arr = zeros(N,1); 
etheta_arr = zeros(N,1);
vc_prev = [0; 0];
int_ev = zeros(1,2);
ev_prev = [0,0];
K_mat = [1/r, -b/(2*r);
         1/r,  b/(2*r)];
VD_prev = [0;0];          
iw_lim = 20;               

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
    theta = q(3);

    ref_x = 2*cos(t);
    ref_y = 2*sin(t);
    ref_x_dot = -2*sin(t);
    ref_y_dot =  2*cos(t);
    ref_theta = atan2(ref_y_dot, ref_x_dot);

    vr = 2; wr = 1;
    eq = [ref_x - q(1); ref_y - q(2); wrapToPi(ref_theta - q(3))];
    R_theta = [cos(theta), sin(theta), 0; -sin(theta), cos(theta), 0; 0, 0, 1];
    ec = R_theta * eq;

    ex_arr(i) = ec(1); 
    ey_arr(i) = ec(2); 
    etheta_arr(i) = ec(3);

    vd = k1*ec(1) + vr*cos(ec(3));
    wd = wr + k2*vr*ec(2) + k3*vr*sin(ec(3));
    vc = [vd; wd];
 
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

    dt = tspan(i) - tspan(i-1);
    if i > 2
        int_ev = int_ev + 0.5*((ev' + ev_prev))*dt;         
        int_ev = max(-iw_lim*ones(1,2), min(iw_lim*ones(1,2), int_ev));
    end
    ev_prev = ev';                                            

    if i == 2
        vcdot = [0;0];   VDdot = [0;0];
    else
        vcdot = (vc - vc_prev)/dt;
        VDdot = (VD - VD_prev)/dt;
    end
    vc_prev = vc;
    VD_prev = VD;

    a_des_body = vcdot + Kp_BS*ev + Ki_BS*(int_ev') + Kd_BS*(-VDdot);  
    a_des_phi  = K_mat * a_des_body;                                 

    P_mat = J' * E;
    rhs   = J'*( M*J*a_des_phi + M*dJ*[dq(4);dq(5)] + Vm*dq  ); 
    tau   = pinv(P_mat) * rhs;              
    tau   = max(-u_sat, min(u_sat, tau));              

    tau1_arr(i) = tau(1);
    tau2_arr(i) = tau(2);
    vd_arr(i) = vd_actual;
    wd_arr(i) = wd_actual;

    V     = [dq(4); dq(5)];  
    v_dot = (J'*M*J)\( J'*( -M*dJ*V - Vm*dq + (E*tau + tau_d) ) );

    dq_new = dq;
    dq_new(4:5) = dq(4:5) + v_dot * dt;

    theta_new = q(3) + r/b * (dq_new(5) - dq_new(4)) * dt;
    x_new     = q(1) + r/2 * cos(theta_new) * (dq_new(4) + dq_new(5)) * dt;
    y_new     = q(2) + r/2 * sin(theta_new) * (dq_new(4) + dq_new(5)) * dt;
    phi1_new  = q(4) + dq_new(4) * dt;
    phi2_new  = q(5) + dq_new(5) * dt;
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
legend('\tau_(Left)','\tau_(Right)'); xlabel('Time (s)'); ylabel('Torque (Nm)');
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

    set(headArrow, 'XData', x, 'YData', y, ...
                   'UData', headLen*cos(th), 'VData', headLen*sin(th));

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
