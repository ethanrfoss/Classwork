function [T,gripState] = TrajectoryGenerator(Tseinitial,Tscinitial,Tscfinal,Tcegrasp,Tcestandoff,k)
%% Function Descirption:
% For this component, you will write a function called TrajectoryGenerator to create the reference (desired) tra-
% jectory for the end-effector frame {e}. This trajectory should consist of eight concatenated trajectory segments, de-
% scribed below. Each trajectory segment begins and ends at rest. This function is likely to use either ScrewTrajectory
% or CartesianTrajectory from the Modern Robotics code library.
%% Inputs:
% Tseinitial - The initial configuration of the end-effector
% Tscinitial - The initial configuration of the cube
% Tscfinal - The desired final configuration of the cube
% Tcegrasp - The configuration of the end-effector relative to the cube while grasping
% Tcestandoff - The standoff configuration of the end-effector above the cube, before and after grasping, relative to the cube
% k - The number of trajectory reference configurations per 0.01 seconds
%% Outputs:
% T - Cell array of Transformations at each time step
% gripState - State of gripper. 1 for closed, 0 for open
%% Trajectory Steps:
% 1. Move the gripper from its initial configuration to a ”standoff” configuration a few cm above the block.
% 2. Move the gripper down to the grasp position.
% 3. Close the gripper.
% 4. Move the gripper back up to the ”standoff” configuration.
% 5. Move the gripper to a ”standoff” configuration above the final configuration.
% 6. Move the gripper to the final configuration of the object.
% 7. Open the gripper.
% 8. Move the gripper back to the ”standoff” configuration.

%% Trajectory 1: Move the gripper from its initial configuration to a ”standoff” configuration a few cm above the block. 0 to 10 sec
Tf = 10; N = Tf*k/.01;
T = ScrewTrajectory(Tseinitial, Tscinitial*Tcestandoff, Tf, N, 5);
gripState = zeros(1,N);
%% Trajectory 2: Move the gripper down to the grasp position. 10 to 12 sec
Tf = 2; N = Tf*k/.01;
Tnew = ScrewTrajectory(Tscinitial*Tcestandoff, Tscinitial*Tcegrasp, Tf, N, 5);
T = [T Tnew];
gripState = [gripState zeros(1,N)];
%% Trajectory 3: Close the gripper. 12 to 13 sec
Tf = 1; N = Tf*k/.01; % Takes about .65 seconds to close the gripper
Tnew = cell(1,N);
for i = 1:N
    Tnew{i} = T{end};
end
T = [T Tnew];
gripState = [gripState ones(1,N)];
%% Trajectory 4: Move the gripper back up to the ”standoff” configuration. 13 to 15 sec
Tf = 2; N = Tf*k/.01;
Tnew = ScrewTrajectory(Tscinitial*Tcegrasp,Tscinitial*Tcestandoff, Tf, N, 5);
T = [T Tnew];
gripState = [gripState ones(1,N)];
%% Trajectory 5: Move the gripper to a ”standoff” configuration above the final configuration.15 to 25 sec
Tf = 10; N = Tf*k/.01;
Tnew = ScrewTrajectory(Tscinitial*Tcestandoff,Tscfinal*Tcestandoff, Tf, N, 5);
T = [T Tnew];
gripState = [gripState ones(1,N)];
%% Trajectory 6: Move the gripper to the final configuration of the object. 25 to 27 sec
Tf = 2; N = Tf*k/.01;
Tnew = ScrewTrajectory(Tscfinal*Tcestandoff,Tscfinal*Tcegrasp, Tf, N, 5);
T = [T Tnew];
gripState = [gripState ones(1,N)];
%% Trajectory 7: Open the gripper. 27 to 28 sec
Tf = 1; N = Tf*k/.01; % Takes about .65 seconds to close the gripper
Tnew = cell(1,N);
for i = 1:N
    Tnew{i} = T{end};
end
T = [T Tnew];
gripState = [gripState zeros(1,N)];
%% Trajectory 8: Move the gripper back to the ”standoff” configuration. 28 to 30 sec
Tf = 2; N = Tf*k/.01;
Tnew = ScrewTrajectory(Tscfinal*Tcegrasp,Tscfinal*Tcestandoff, Tf, N, 5);
T = [T Tnew];
gripState = [gripState zeros(1,N)];

%% Save Transformation States to Matrix and Save to File
for i = 1:length(T)
    Tse(i,:) = [T{i}(1,1) T{i}(1,2) T{i}(1,3) T{i}(2,1) T{i}(2,2) T{i}(2,3) T{i}(3,1) T{i}(3,2) T{i}(3,3) T{i}(1,4) T{i}(2,4) T{i}(3,4) gripState(i)];
end
writematrix(Tse,'C:\Users\ethan\OneDrive\Documents\MATLAB\MAE 204\FinalProject\Tse.csv');


%% Copellia Viedo Link
% https://drive.google.com/drive/folders/1rGz7sa2HcHKVPKZYq9HZyhNNsDw576cz?usp=sharing
end
