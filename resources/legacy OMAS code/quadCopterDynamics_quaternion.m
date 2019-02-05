
% UPDATE THE QUADROTOR STATE VECTOR
function [dXdt] = quadCopterDynamics_quaternion(obj,in,pn_desired,psi_desired)

% State variables
pn=in(1:3);
vn=in(4:6);
qn=in(7:10);
omega=in(11:13);

% MODEL DYNAMIC PROPERTIES
g = obj.DYNAMICS.g;
e3 = obj.DYNAMICS.e3;
m = obj.DYNAMICS.m;
M = obj.DYNAMICS.M;
K = obj.DYNAMICS.K;

% control law
%----------case 1: equilibrium input + noise
% f=m*g+4*randn/2;
% % tau=zeros(3,1)+[0.1 0 0]'; % will rotate and crash!
% tau=zeros(3,1)+0.5*randn(3,1)/2; % will rotate and crash!
%----------case 2: linear LQR based on the linearized model

x_desireConst = [pn_desired;zeros(3,1);[0;0;psi_desired];zeros(3,1)];

R = OMAS_axisTools.quaternionToRotationMatrix(qn); % By default, q is mapping G to B
R = R';

% THE ROTATION ANGLES OF THE BODY
[eta] = obj.fcn_EulerFromRotation(R);

delta_x = [pn;vn;eta;omega] - x_desireConst;
delta_u = -K*delta_x; % delta u=K*delta x

u   = [m*g;0;0;0] + delta_u;
f   = u(1);  %+randn/10;
tau = u(2:4);%+randn/10; % < VALID UNTIL HERE

% nonlinear dynamics based on rotation dynamics
pn_dot = vn;
vn_dot = f/m*R*e3-g*e3;
q_dot = OMAS_axisTools.quaternionDifferential(qn,omega);
omega_dot = inv(M)*(tau-obj.fcn_skewSymmetric(omega)*M*omega);

% INTEGRATE THE STATES
dXdt = [pn_dot;vn_dot;q_dot;omega_dot];
end