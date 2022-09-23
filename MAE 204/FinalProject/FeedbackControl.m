function [Ve,Xerr] = FeedbackControl(Tse,Tsed,Tsednext,Kp,Ki,dt)
%% Function Description:
% The third component you will write is a function, FeedbackControl that calculates the task-space feedforward
% plus feedback control law discussed in class, and also appearing as Equation (11.16) and (13.37) in the textbook.
%% Inputs:
% The current actual end-effector configuration X (aka Tse)
% The current reference end-effector configuration Xd (aka Tse,d)
% The reference end-effector configuration at the next timestep, Xd,next (aka Tse,d,next)
% The PI gain matrices Kp and Ki
% The timestep ∆t between reference trajectory configurations
%% Outputs:
% The commanded end-effector twist V expressed in the end-effector frame {e}.

%% Compute Error
Xerr = se3ToVec(MatrixLog6(Tse^-1*Tsed));

%% Integral of Error
global intXerr;
if isempty(intXerr)
    intXerr = Xerr*dt;
else
    intXerr = intXerr + Xerr*dt;
end

%% Compute End-Effector Twist
Ve = Adjoint(Tse^-1*Tsed)*1/dt*se3ToVec(MatrixLog6(Tsed^-1*Tsednext))+Kp*Xerr+Ki*intXerr;
%Ve = Kp*se3ToVec(Xerr)+Ki*se3ToVec(intXerr);

end