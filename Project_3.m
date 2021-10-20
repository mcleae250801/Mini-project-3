close all; clear all;

zz=zeros(3,1); ex = [1;0;0]; ey = [0;1;0]; ez = [0;0;1];

load S_sphere_path

%% plot sphere
% plot the spherical S
figure(1);plot3(p_S(1,:),p_S(2,:),p_S(3,:),'rx','linewidth',3);
%{
p_s(1,:) = x 
p_s(2,:) = y
p_s(3,:) = z
%}

xlabel('x');ylabel('y');zlabel('z');
hold on;
% add a 3d sphere
surf(X,Y,Z)
% make it transparent
alpha .5
axis(r*[-1 1 -1 1 0 2]);axis('square');
view(120,10);

% convert to equal path length grid
diffS=vecnorm(diff(p_S')');
%path length
ls=[0 cumsum(diffS)];
%final path length
lf=sum(diffS);
N=100;
l=(0:lf/N:lf);

pS=interp1(ls,p_S',l,'spline')';
% plot it out again with equal path length
figure(2);plot3(pS(1,:),pS(2,:),pS(3,:),'rx','linewidth',3);
xlabel('x');ylabel('y');zlabel('z');
hold on;
% 3d sphere
surf(X,Y,Z)
% make it transparent
alpha 0.4
axis(r*[-1 1 -1 1 0 2]);axis('square');
view(120,10);

% check the path length is indeed equal
dlvec=vecnorm(diff(pS')');
figure(3);plot(dlvec,'x')
dl=mean(dlvec);
disp(max(abs(dlvec-dl)));

% save it as the path file
save S_sphere_path_uniform l pS



%% Part 1
%END EFFECTOR FRAME
% find the end effector frame
pc=r*ez; %center of the ball
N=length(pS);
xT=zeros(3,N);zT=zeros(3,N);yT=zeros(3,N);
quat=zeros(4,N); % allocate space for unit quaternion representation of R_0T}

B_1 = [];
B_2 = [];
BC =[]

for i=1:N
    xT(:,i)=(pS(:,i)-pc)/norm(pS(:,i)-pc);
    if i<N
        yT(:,i)=(pS(:,i+1)-pS(:,i));
    else
        yT(:,i)=yT(:,i-1);
    end
    yT(:,i)=yT(:,i)-yT(:,i)'*xT(:,i)*xT(:,i);
    yT(:,i)=yT(:,i)/norm(yT(:,i));
    zT(:,i)=cross(xT(:,i),yT(:,i));
    R=[xT(:,i) yT(:,i) zT(:,i)];
    quat(:,i)=R2q(R);
     
    [B1,B2] = subprob2(-ez, ey, R*ex, ex);
    B3_1 = subprob1(ex, ez, (rot(ey,-B2(1))*rot(ez,-B1(1))*R*ez));
    B3_2 = subprob1(ex, ez, (rot(ey,-B2(2))*rot(ez,-B1(2))*R*ez));
   
    
    
    BS_1 = [B1(1); B2(1); B3_1];
    BS_2 = [B1(2); B2(2); B3_2];
    
    B_1 = [B_1 BS_1];
    B_2 = [B_2 BS_2];
    
    BC = [BC vee(logm(R))];
    
end

figure(5)
plot(ls, quat(1,:))
hold on
title("Quaternion vector")
plot(ls, quat(2,:))
plot(ls, quat(3,:))
plot(ls, quat(4,:))
ylabel('')
xlabel('lambda')
legend('q1', 'q2', 'q3', 'q4')



figure(6)
plot(ls,B_1(1,:))
hold on
title("Euler Solution 1")
plot(ls,B_1(2,:))
plot(ls,B_1(3,:))
xlabel('lambda')
ylabel('rads')
legend('Β1', 'Β2', 'Β3')


figure(7)
plot(ls,B_2(1,:))
title("Euler Solution 2")
hold on
plot(ls,B_2(2,:))
plot(ls,B_2(3,:))
xlabel('lambda')
ylabel('rads')
legend('Β1', 'Β2', 'Β3')


figure(8)
plot(ls,BC(1,:))
title("Euler Solution using B=theta*k")
hold on
plot(ls,BC(2,:))
plot(ls,BC(3,:))
xlabel('lambda')
ylabel('rads')
legend('Β1', 'Β2', 'Β3')


% plot out the end effector frame
m=5;
%matlab plottransform command plots a frame at a given location
figure(1);h=plotTransforms(pS(:,1:m:end)',quat(:,1:m:end)')

set(h,'LineWidth',.5);

% ABB IRB 1200 parameters

L1=399.1;
L2=348;
L3=42;
L4=451;
L5=82;

% P
p01=0*ex+L1*ez;
p12=zz;
p23=L2*ez;
p34=L3*ez+L4*ex;
p45=zz;
p56=zz;
p6T=L5*ex;


% H
h1=ez;
h2=ey;
h3=ey;
h4=ex;
h5=ey;
h6=ex;


% Final transformation

%% part 2
% define abb 1200 robot using POE convention
irb1200.P=[p01 p12 p23 p34 p45 p56 p6T]/1000;
irb1200.H=[h1 h2 h3 h4 h5 h6];
irb1200.joint_type=[0 0 0 0 0 0];
%irb1200.R6T=R6T;

% define collision body for abb 1200
radius=.01;
[irb1200_rbt,colLink]=defineRobot(irb1200,radius);
 
% 
% choose the inverse kinematics solution
%

for i=1:N
    % specify end effector SE(3) frame
    %Td{i}=[[xT(:,i) yT(:,i) zT(:,i)]*R6T' pS(:,i);[0 0 0 1]];
    Td{i}=[[xT(:,i) yT(:,i) zT(:,i)] pS(:,i);[0 0 0 1]];
    irb1200.T=Td{i};
    %
    irb1200=invkinelbow(irb1200); % << you need to supply this!!!
    %
    for k=1:8
        q(:,i,k)=irb1200.q(:,k);
    end
        % check forward kinematics to make sure the IK solution is correct
    for k=1:8
        irb1200.q=q(:,i,k);
        irb1200=fwddiffkiniter(irb1200);
        T{i,k}=irb1200.T;
    end

end

% choose the pose to visualize
ksol=1

for i=1:N
    % show robot pose (ever 5 frames)
    if mod(i,5)==0
        disp(norm(T{i,ksol}-Td{i}));
        figure(2);show(irb1200_rbt,q(:,i,ksol),'collision','on');
        view(150,10);
    end
end
% compute end effector linear and angular velcoity 

lsdot=.01;
for i=1:N-1
    dt(i) = (ls(i+1)-ls(i))/lsdot;
    for k=1:8
        qdot(:,i,k)=(q(:,i+1,k)-q(:,i,k))/dt(i);
        Ri1=T{i+1,k}(1:3,1:3);
        Ri=T{i,k}(1:3,1:3);
        w(:,i,k)=vee(Ri1*Ri'-eye(3,3))/dt(i);
        pi1=T{i+1,k}(1:3,4);
        pi=T{i,k}(1:3,4);
        v(:,i,k)=(pi1-pi)/dt(i);
    end
end

for k=1:8
   for i=1:6
       maxqdot(i,k)=max(qdot(i,:,k));
   end
   fprintf('maximum qdot for pose %d \n', k);
   disp(maxqdot(:,k)');   
end


%% Part 3
maxV = zeros(1,8);
sumq3 = zeros(1,8);
for i=1:8
    for j=1:N
        if (max(abs(q(:,j,i))) > maxV(i))
            maxV(i) = max(abs(q(:,j,i)));
        end
        sumq3(i) = sumq3(i) + q(3,j,i);
    end
end

min_q3 = find(sumq3==min(sumq3));
disp("Pose(s) with the smallest variation of q3:")
fprintf("Pose %i \n", min_q3);

min_max = find(maxV==min(maxV));
disp("Pose(s) with the fastest path speed:")
fprintf("Pose %i\n",min_max);


%% fuctions

function q=R2q(R)
 q=zeros(4,1);
  q(1)=.5*sqrt(trace(R)+1);
  if abs(q(1))<1e-5
    [k,theta]=R2kth(R);
    q(2:4)=k;
  else
    q(2:4)=vee(R-R')/4/q(1);
  end
end



function robot = invkinelbow(robot)
    ex = [1;0;0]; ey = [0;1;0]; ez = [0;0;1];
    T=robot.T;
    R0T = T(1:3,1:3);
    p0T = T(1:3,4);
    
    h1=ez;
    h2=ey;
    h3=ey;
    h4=ex;
    h5=ey;
    h6=ex;

    
    p01 = robot.P(:,1);
    p12 = robot.P(:,2);
    p23 = robot.P(:,3);
    p34 = robot.P(:,4);
    p45 = robot.P(:,5);
    p56 = robot.P(:,6);
    p6T = robot.P(:,7);
    
    
   p16_0 = p0T-p01-R0T*p6T;
   
   %solves for q3_1 and q3_2
   q3 = subprob3(h3,-p34,p23,norm(p16_0));
   %q3_1 = q3(1);
   %q3_2 = q3(2);
   
   %         sol 1      sol 2       sol 1      sol 2
   q3 = [q3(1) q3(1) q3(2) q3(2) q3(1) q3(1) q3(2) q3(2)];
   
   %solves Rotataion y matrix for q3_1
   %Ry_q3_1 = rot(ey,q3_1);
   %solves Rotation y matrix for q3_2
   %Ry_q3_2 = rot(ey,q3_2);
   %p2 for subproblem
   %p2_1 = (p23+Ry_q3_1*p34);
   %p2_2 = (p23+Ry_q3_2*p34);
   
   %solves for solution 1 of q1 and q2
   [q1_1,q2_1] = subprob2(-h1,h2,p16_0,(p23+(rot(h3,q3(1))*p34)));
   %solves for solution 2 of q1 and q2
   [q1_2,q2_2] = subprob2(-h1,h2,p16_0,(p23+(rot(h3,q3(3))*p34)));
   
   %      sol 1   sol 2  sol 3   sol 4        
   q1 = [q1_1(1) q1_1(2) q1_2(1) q1_2(2) q1_1(1) q1_1(2) q1_2(1)  q1_2(2)];
   q2 = [q2_1(1) q2_1(2) q2_2(1) q2_2(2) q2_1(1) q2_1(2) q2_2(1)  q2_2(2)];
   
   R03_1 = rot(h1,q1(1))* rot(h2,q2(1))*rot(h3, q3(1));
   R03_2 = rot(h1,q1(2))* rot(h2,q2(2))*rot(h3, q3(2));
   R03_3 = rot(h1,q1(3))* rot(h2,q2(3))*rot(h3, q3(3));
   R03_4 = rot(h1,q1(4))* rot(h2,q2(4))*rot(h3, q3(4));
   
   
   %q4 and q5 using solution 1 
   [q4_1,q5_1] = subprob2(-h4,h5,(R03_1')* R0T * h6,h6);
   [q4_2,q5_2] = subprob2(-h4,h5,(R03_2')* R0T * h6,h6);
   [q4_3,q5_3] = subprob2(-h4,h5,(R03_3')* R0T * h6,h6);
   [q4_4,q5_4] = subprob2(-h4,h5,(R03_4')* R0T * h6,h6);
   
   q4 = [q4_1(1) q4_2(1) q4_3(1) q4_4(1) q4_1(2) q4_2(2) q4_3(2) q4_4(2)];
   q5 = [q5_1(1) q5_2(1) q5_3(1) q5_4(1) q5_1(2) q5_2(2) q5_3(2) q5_4(2)];
   
   
   P2_1 = R03_1*rot(h4,q4(1))*rot(h5,(q5(1)));
   P2_2 = R03_2*rot(h4,q4(2))*rot(h5,(q5(2))); 
   P2_3 = R03_3*rot(h4,q4(3))*rot(h5,(q5(3)));
   P2_4 = R03_4*rot(h4,q4(4))*rot(h5,(q5(4)));
   P2_5 = R03_1*rot(h4,q4(5))*rot(h5,(q5(5)));
   P2_6 = R03_2*rot(h4,q4(6))*rot(h5,(q5(6)));
   P2_7 = R03_3*rot(h4,q4(7))*rot(h5,(q5(7)));
   P2_8 = R03_4*rot(h4,q4(8))*rot(h5,(q5(8)));
   
   q6_1 = subprob1(h6,h5,P2_1'*R0T*h5);
   q6_2 = subprob1(h6,h5,P2_2'*R0T*h5);
   q6_3 = subprob1(h6,h5,P2_3'*R0T*h5);
   q6_4 = subprob1(h6,h5,P2_4'*R0T*h5);
   q6_5 = subprob1(h6,h5,P2_5'*R0T*h5);
   q6_6 = subprob1(h6,h5,P2_6'*R0T*h5);
   q6_7 = subprob1(h6,h5,P2_7'*R0T*h5);
   q6_8 = subprob1(h6,h5,P2_8'*R0T*h5);
   
   q6 = [q6_2 q6_1 q6_4 q6_3 q6_6 q6_5 q6_8 q6_7];
   
   q_sol = [q1; q2; q3; q4; q5; q6];
    
   robot.q = q_sol
        
end

function robot = fwddiffkiniter(robot)
    q=robot.q;
    n=length(robot.q);

    T=eye(4,4);
    joint_type=robot.joint_type;
    if ~exist('joint_type');robot.joint_type=zeros(1,n);end

    for i=1:n
        h=robot.H(1:3,i);
        if robot.joint_type(i)==0
            R=expm(hat(h)*q(i));p=robot.P(1:3,i);
        else
            R=eye(3,3);p=robot.P(1:3,i)+q(i)*h;
        end
        T=T*[R p;zeros(1,3) 1];
    end    
    
    robot.T=T*[eye(3,3) robot.P(1:3,n+1);zeros(1,3) 1];

end

function R = rot(k,theta)
    k=k/norm(k);
    R=eye(3,3)+sin(theta)*hat(k)+(1-cos(theta))*hat(k)*hat(k);
end 



%
% solve for q subtended between p1 and p2
%    k determines the sign of q
%
% input: k,p1,p2 as R^3 column vectors
% output: q (scalar)
%

function q=subprob0(k,p1,p2)

if ((k'*p1)>sqrt(eps)|(k'*p2)>sqrt(eps))
  error('k must be perpendicular to p and q');
end

p1=p1/norm(p1);
p2=p2/norm(p2);

q=2*atan2(norm(p1-p2),norm(p1+p2));

if k'*(cross(p1,p2))<0
  q=-q;
end 
end

%
% q=subprob1(k,p1,p2)
%
% solve for q from
%
% exp(k x q) p1 = p2
%
% input: k,p1,p2 as R^3 column vectors
% output: q (scalar)
%

function q=subprob1(k,p1,p2)

p2=p2/norm(p2)*norm(p1);

if norm(p1-p2)<sqrt(eps);q=0;return;end
  
k=k/norm(k);
pp1=p1-(p1'*k)*k;
pp2=p2-(p2'*k)*k;

epp1=pp1/norm(pp1);
epp2=pp2/norm(pp2);

q=subprob0(k,epp1,epp2);
%q=atan2(k'*(cross(epp1,epp2)),epp1'*epp2);
end

%[q1,q2]=subprob2(k1,k2,p1,p2)
%
% solve for theta1 and theta2 from
%
% exp(k1 x q1) p1 = exp(k2 x q2) p2 
%  
% input: k1,k2,p1,p2 as R^3 column vectors
%
% output: q1 and q2 as 2x1 columns corresponding to the two solutions
%

function [q1,q2]=subprob2(k1,k2,p1,p2)

p2=p2/norm(p2)*norm(p1);
k12=k1'*k2;
pk1=p1'*k1;
pk2=p2'*k2;

% check if solution exists

if abs(k12^2-1)<eps;theta1=[];theta2=[];
    q1=[NaN;NaN];q2=[NaN;NaN];
    disp('no solution (k1 and k2 are collinear)');
    return;
end

a=[1 -k12; -k12 1]*[pk1;pk2]/(1-k12^2);

% 
% check if solution exists
%
cond=(norm(p1)^2-norm(a)^2-2*a(1)*a(2)*k12);

% special case: 1 solution
if abs(cond)<eps;
  v=[k1 k2]*a;
  q1a=subprob1(k1,p1,v);
  q2a=subprob1(k2,p2,v);
  q1=[q1a;q1a];
  q2=[q2a;q2a];
end

% special case: no solution
if cond<0
    q1=[NaN NaN];q2=[NaN NaN];
    disp('no solution (two cones do not intersect)');
    return;
end

gamma=sqrt(cond)/norm(cross(k1,k2));

% general case: 2 solutions

q1=zeros(2,1);
q2=zeros(2,1);

v1=[k1 k2 cross(k1,k2)]*[a;gamma];
v2=[k1 k2 cross(k1,k2)]*[a;-gamma];
q1(1)=subprob1(k1,p1,v1);
q1(2)=subprob1(k1,p1,v2);

q2(1)=subprob1(k2,p2,v1);
q2(2)=subprob1(k2,p2,v2);

end

%
% q=subprob3(k,p1,p2,d)
%
% solve for theta from
%
% norm(p2-exp(k x q) p1) = d
%
% input: k,p1,p2 as R^3 column vectors, delta: scalar
% output: q (2x1 vector, 2 solutions)
%

function q=subprob3(k,p1,p2,d)

pp1=p1-k'*p1*k;
pp2=p2-k'*p2*k;
dpsq=d^2-(k'*(p1-p2))^2;

if dpsq<0;theta=[NaN;NaN];return;end

if dpsq==0;theta=subprob1(k,pp1/norm(pp1),pp2/norm(pp2));return;end
  
bb=(norm(pp1)^2+norm(pp2)^2-dpsq)/(2*norm(pp1)*norm(pp2));
if abs(bb)>1; theta=[NaN;NaN];return;end

phi=acos(bb);

q0=subprob1(k,pp1/norm(pp1),pp2/norm(pp2));
q=zeros(2,1);

q(1)=q0+phi;
q(2)=q0-phi;

end

function khat = hat(k)
  
  khat=[0 -k(3) k(2); k(3) 0 -k(1); -k(2) k(1) 0];

end

function k = vee(K)

k=[-K(2,3);K(1,3);-K(1,2)];

end
