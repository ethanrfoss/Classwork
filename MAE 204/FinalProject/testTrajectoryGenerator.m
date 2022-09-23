function testTrajectoryGenerator
%% Trajectroy Generator Tester

%% Transformation Matrices:
% Tseinitial:
Tseinitial = [0 0 1 0;0 1 0 0;-1 0 0 .5;0 0 0 1];
% Tscinitial:
Tscinitial = [1 0 0 1;0 1 0 0;0 0 1 .025;0 0 0 1];
% Tscfinal:
Tscfinal = [0 1 0 0;-1 0 0 -1;0 0 1 .025;0 0 0 1];
% Tcestandoff:
Tcestandoff = [0 0 1 -.035;0 1 0 0;-1 0 0 .25;0 0 0 1];
% Tcegrasp:
Tcegrasp = [-.5 0 sqrt(3)/2 0;0 1 0 0;-sqrt(3)/2 0 -.5 0;0 0 0 1];

%% Run Trajectroy Generator:
k = 1;
[T,gripState] = TrajectoryGenerator(Tseinitial,Tscinitial,Tscfinal,Tcegrasp,Tcestandoff,k);

%% Plotter:
h = figure; hold on; axis equal; view(45,45);
xlabel('X[dm]');ylabel('Y[dm]');zlabel('Z[dm]');
% Space Frame:
qs = quiver3(0,0,0,1,0,0,1,'k');
quiver3(0,0,0,0,1,0,1,'k');
quiver3(0,0,0,0,0,1,1,'k');
% Tseinital:
qsei = quiver3(Tseinitial(1,4)*10,Tseinitial(2,4)*10,Tseinitial(3,4)*10,Tseinitial(1,1),Tseinitial(2,1),Tseinitial(3,1),1,'g');
quiver3(Tseinitial(1,4)*10,Tseinitial(2,4)*10,Tseinitial(3,4)*10,Tseinitial(1,2),Tseinitial(2,2),Tseinitial(3,2),1,'g');
quiver3(Tseinitial(1,4)*10,Tseinitial(2,4)*10,Tseinitial(3,4)*10,Tseinitial(1,3),Tseinitial(2,3),Tseinitial(3,3),1,'g');
% Tscinitial:
qsci = quiver3(Tscinitial(1,4)*10,Tscinitial(2,4)*10,Tscinitial(3,4)*10,Tscinitial(1,1),Tscinitial(2,1),Tscinitial(3,1),1,'y');
quiver3(Tscinitial(1,4)*10,Tscinitial(2,4)*10,Tscinitial(3,4)*10,Tscinitial(1,2),Tscinitial(2,2),Tscinitial(3,2),1,'y');
quiver3(Tscinitial(1,4)*10,Tscinitial(2,4)*10,Tscinitial(3,4)*10,Tscinitial(1,3),Tscinitial(2,3),Tscinitial(3,3),1,'y');
% Tscfinal:
qscf = quiver3(Tscfinal(1,4)*10,Tscfinal(2,4)*10,Tscfinal(3,4)*10,Tscfinal(1,1),Tscfinal(2,1),Tscfinal(3,1),1,'r');
quiver3(Tscfinal(1,4)*10,Tscfinal(2,4)*10,Tscfinal(3,4)*10,Tscfinal(1,2),Tscfinal(2,2),Tscfinal(3,2),1,'r');
quiver3(Tscfinal(1,4)*10,Tscfinal(2,4)*10,Tscfinal(3,4)*10,Tscfinal(1,3),Tscfinal(2,3),Tscfinal(3,3),1,'r');
% % Tcestandoff:
% qsei = quiver3(Tcestandoff(1,4)*10,Tcestandoff(2,4)*10,Tcestandoff(3,4)*10,Tcestandoff(1,1),Tcestandoff(2,1),Tcestandoff(3,1),1,'g');
% quiver3(Tcestandoff(1,4)*10,Tcestandoff(2,4)*10,Tcestandoff(3,4)*10,Tcestandoff(1,2),Tcestandoff(2,2),Tcestandoff(3,2),1,'g');
% quiver3(Tcestandoff(1,4)*10,Tcestandoff(2,4)*10,Tcestandoff(3,4)*10,Tcestandoff(1,3),Tcestandoff(2,3),Tcestandoff(3,3),1,'g');
% % Tcegrasp:
% qsei = quiver3(Tcegrasp(1,4)*10,Tcegrasp(2,4)*10,Tcegrasp(3,4)*10,Tcegrasp(1,1),Tcegrasp(2,1),Tcegrasp(3,1),1,'g');
% quiver3(Tcegrasp(1,4)*10,Tcegrasp(2,4)*10,Tcegrasp(3,4)*10,Tcegrasp(1,2),Tcegrasp(2,2),Tcegrasp(3,2),1,'g');
% quiver3(Tcegrasp(1,4)*10,Tcegrasp(2,4)*10,Tcegrasp(3,4)*10,Tcegrasp(1,3),Tcegrasp(2,3),Tcegrasp(3,3),1,'g');

%% Cubes:
cubeX = [.25 .25 .25 .25 -.25 -.25 -.25 -.25]; cubeY = [.25 .25 -.25 -.25 .25 .25 -.25 -.25]; cubeZ = [.25 -.25 .25 -.25 .25 -.25 .25 -.25];
cubehull = convhull(cubeX,cubeY,cubeZ);
tC = trisurf(cubehull,cubeX+Tscinitial(1,4)*10,cubeY+Tscinitial(2,4)*10,cubeZ+Tscinitial(3,4)*10,'FaceColor','g','EdgeColor','g');
%% Iterative plot:
i = 1;
qT1 = quiver3(T{i}(1,4)*10,T{i}(2,4)*10,T{i}(3,4)*10,T{i}(1,1),T{i}(2,1),T{i}(3,1),1,'b');
qT2 = quiver3(T{i}(1,4)*10,T{i}(2,4)*10,T{i}(3,4)*10,T{i}(1,2),T{i}(2,2),T{i}(3,2),1,'b');
qT3 = quiver3(T{i}(1,4)*10,T{i}(2,4)*10,T{i}(3,4)*10,T{i}(1,3),T{i}(2,3),T{i}(3,3),1,'b');
frame = getframe(h);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
imwrite(imind,cm,[cd '\SavedGifs\' 'Robot' '.gif'],'gif','DelayTime',0, 'Loopcount',inf); 

for i = 2:length(T)
    set(qT1,'XData',T{i}(1,4)*10,'YData',T{i}(2,4)*10,'ZData',T{i}(3,4)*10,'UData',T{i}(1,1),'VData',T{i}(2,1),'WData',T{i}(3,1));
    set(qT2,'XData',T{i}(1,4)*10,'YData',T{i}(2,4)*10,'ZData',T{i}(3,4)*10,'UData',T{i}(1,2),'VData',T{i}(2,2),'WData',T{i}(3,2));
    set(qT3,'XData',T{i}(1,4)*10,'YData',T{i}(2,4)*10,'ZData',T{i}(3,4)*10,'UData',T{i}(1,3),'VData',T{i}(2,3),'WData',T{i}(3,3));
    drawnow;
    if gripState(i) == 1
        Tsc = T{i}*(Tcegrasp)^-1;
        Tsc(1:3,4) = Tsc(1:3,4)*10;
        C = Tsc*[cubeX;cubeY;cubeZ;ones(1,length(cubeX))];
        set(tC,'XData',C(1,:),'YData',C(2,:),'ZData',C(3,:));
    end
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    imwrite(imind,cm,[cd '\SavedGifs\' 'Robot' '.gif'],'gif','DelayTime',0,'WriteMode','append'); 
end