function plotPend(xk,s,saveName)
%% Shapes for plotting:
clear body pend1 pend2;
body(1) = polyshape([3/4*s.L1 3/4*s.L1 -3/4*s.L1 -3/4*s.L1],[.25*s.L2 -.25*s.L2 -.25*s.L2 .25*s.L2]);
thickness = s.L1/32;
pend1(1) = polyshape([thickness thickness -thickness -thickness],[0 s.L1 s.L1 0]);
pend2(1) = polyshape([thickness thickness -thickness -thickness],[0 s.L2 s.L2 0]);
for i = 2:length(xk(1,:))
    if isnan(xk(1,i))
        break;
    end
    body(i).Vertices = body(1).Vertices +[xk(1,i) 0];
    pend1(i).Vertices = rotate(pend1(1),xk(2,i)*180/pi).Vertices + [xk(1,i)+s.L1/4 0];
    pend2(i).Vertices = rotate(pend2(1),xk(3,i)*180/pi).Vertices + [xk(1,i)-s.L1/4 0];
end

h = figure; hold on; axis equal; tit = title({'Full System - LTV Feedback and Estimation on Finite Horizon', '+ LQR and Ricatti on Infinite Horizon',sprintf('t = 0')});
xlim([min(xk(1,1))-s.L1-.2 max(xk(1,1))+s.L1+.2]); ylim([-s.L1-.2, s.L1+.2]);
plot([min(xk(1,:))-s.L1-.2 max(xk(1,:))+s.L1+.2],[0 0]);
pb = plot(body(1));
pp1 = plot(pend1(1));
pp2 = plot(pend2(1));
for i = 1:length(body)
    xlim([min(xk(1,i))-s.L1-.2 max(xk(1,i))+s.L1+.2]);
    set(pb,'Shape',body(i));
    set(pp1,'Shape',pend1(i));
    set(pp2,'Shape',pend2(i));
    drawnow;
    
    set(tit,'String',{'Full System - LTV Feedback and Estimation on Finite Horizon', '+ LQR and Ricatti on Infinite Horizon',sprintf('t = %.2f',(i-1)*.01)});
    if nargin>2
        frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if i == 1
            imwrite(imind,cm,['C:\Users\ethan\OneDrive\Documents\MATLAB\MAE 200\Final Project\' saveName '.gif'],'gif','DelayTime',0, 'Loopcount',inf); 
        else
            imwrite(imind,cm,['C:\Users\ethan\OneDrive\Documents\MATLAB\MAE 200\Final Project\' saveName '.gif'],'gif','DelayTime',0,'WriteMode','append'); 
        end
    end
end
end