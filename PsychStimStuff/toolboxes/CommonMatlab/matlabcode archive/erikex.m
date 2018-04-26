                                                                     
                                                                     
                                                                     
                                             
%clean up whatever's going on
clf
clear all
clc

%define a unit circle
pts = 0:.05:2*pi;
xPts=sin(pts);
yPts=cos(pts);
subplot(1,2,1)
plot(xPts,yPts,'.')
axis(3*[-1 1 -1 1])
axis equal
grid on

%define a random matrix
A=1.5*randn(2);

%set to 1 to enforce symmetric matrix
sym=0
if sym
    A=(A+A')/2;
end
A(2,1)=0;
%print out the matrix
A
title(sprintf('A=[%g %g; %g %g]',A(1),A(2),A(3),A(4)))

%do the transform
newPts=A*[xPts;yPts];
% newPts=A*newPts;

hold on
plot(newPts(1,:),newPts(2,:),'r.')

%set to 1 to draw the connectors between original points and their images
connect=1;
if connect
    for i=1:length(pts)
        plot([xPts(i) newPts(1,i)],[yPts(i) newPts(2,i)],'k')
    end
end

%set to 1 to draw the eigenbasis (only works for symmetric A)
drawEig=1;
% if sym & drawEig
if drawEig
    
    %do spectral decomposition and sort by eigenvalue magnitude
    [vects vals]=eig(A);
    vals=diag(vals);
    [sortVals order]=sort(abs(vals));
    order=flipud(order);
    sortVals=vals(order);
    vects=vects(:,order);

    %show eigenvectors
    plot([0 vects(1,1)],[0 vects(2,1)],'b','LineWidth',7)
    plot([0 vects(1,2)],[0 vects(2,2)],'b','LineWidth',7)    

    %set to 1 to scale eigenvectors by their eigenvalues
    drawVal=1;
    if drawVal
        plot([0 sortVals(1)*vects(1,1)],[0 sortVals(1)*vects(2,1)],'r','LineWidth',3)
        plot([0 sortVals(2)*vects(1,2)],[0 sortVals(2)*vects(2,2)],'r','LineWidth',3)     
    end

    %set to 1 to remove correlations by switching to eigenbasis
    whiten=0;
    if whiten
        whitened=vects'*newPts;
        plot(whitened(1,:),whitened(2,:),'g.')

        %set to 1 to undo scaling action of A
        normalize=0;
        if normalize
            normalized=diag(1./sortVals)*whitened; %the inverse of a diagonal matrix is the same matrix with the reciprocal of each entry on the diagonal
            plot(normalized(1,:),normalized(2,:),'k.','MarkerSize',20)
        end
        
        %set to 1 to undo rotation action of A -- this completes the
        %inversion of A
        rerotate=0;
        if rerotate
            rerotated=vects*normalized;
            plot(rerotated(1,:),rerotated(2,:),'y.','MarkerSize',20)
        end
    end
end

% subplot(1,2,2)
% 
% pauseTime=.01;
% numRuns=200;
% numSteps=1000;
% stepSize=.01;
% plotPct=.1;
% 
% inds=1:1/plotPct:numSteps;
% 
% currSize=3;
% for j=1:numRuns
%     x=randn(2,1);
%     y=zeros(2,numSteps);
%     for i=1:numSteps
%         y(:,i)=x;
%         x = stepSize*A*x + x;
%     end
%     %colors go blue->red
%     scatter(y(1,inds),y(2,inds),5,1:length(inds));
%     newSize=mean(mean(abs(y)))/100;
%     if newSize>currSize
%         %currSize=newSize;
%     end
%     axis(currSize*[-1 1 -1 1])
%     axis square
%     hold on
%     pause(pauseTime)
% end
    
% discussion points

%note 4 free parameters in 2x2 operator
% - can pick a direction
% - stretch/flip along it
% - stretch/flip along the orthogonal
% - "twist" of points (which maps to which)

%only 3 free parameters in 2x2 symmetric matrix
% - so which parts of above parameters is lost?
% - equivalent question -- what can general operators do that symmetric
%   ones can't?
% - how is this reflected in real vs. complex eigenvalues/vectors?

%eigenvectors of a symmetric matrix form an orthonormal basis
%an orthonormal matrix is a rotation matrix.  why?
%what do eigenvalues near 0 mean?
%what do negative eigenvalues mean?
%why is A'+A symmetric for any A?


