% MEC313 - Project 2
% ------------------------------
% Code written by Andrei Roibu
% Created: 15/04/2019
% Last Modified: 12/05/2019
% ------------------------------
% This script perfoms the relevant symbolic, numerical and statistical
% calculations for the MEC313 Project 2
% ------------------------------

% We start by clearing the workspace

clear all;

% In order to do a sanity check of the hand calculations for the stuffiness
% matrix of a rectangular element with one Gauss point, we perform a
% symbolic calculation and print the result. 

syms w h E niu t
D = [1,niu,0;niu,1,0;0,0,(1-niu)/2];
D = E/(1-niu^2)*D;
J = [w/2,0;0,h/2];
detJ = det(J);
B = [-1/(2*w),0,1/(2*w),0,1/(2*w),0,-1/(2*w),0;
    0,-1/(2*h),0,-1/(2*h),0,1/(2*h),0,1/(2*h);
    -1/(2*h),-1/(2*w),-1/(2*h),1/(2*w),1/(2*h),1/(2*w),1/(2*h),-1/(2*w)];
H = 4;
K = t * H * det(J) * transpose(B) * D * B

k1 = K;
v1 = [1,2,3,4,5,6,7,8]; % The position vectors
k2 = K;
v2 = [7,8,5,6,11,12,9,10];
k3 = sym(zeros(12));
for i = v1
    for j= v1
        k3(i,j) = k3(i,j) + k1(find(v1 == i),find(v1 ==j));
    end
end

for i = v2
    for j= v2
        k3(i,j) = k3(i,j) + k2(find(v2 == i),find(v2 ==j));
    end
end
k3

% Once the symbolic calculation is done, we calculate the solution of the
% displacement and strain distributions using two rectangular elements. 

w = 3.477e-2;
h = 8.455e-3;
E = 71e9;
niu = 0.3;
D = [1,niu,0;niu,1,0;0,0,(1-niu)/2];
D = E/(1-niu^2)*D;
J = [w/2,0;0,h/2];
detJ = det(J);
B = [-1/(2*w),0,1/(2*w),0,1/(2*w),0,-1/(2*w),0;
    0,-1/(2*h),0,-1/(2*h),0,1/(2*h),0,1/(2*h);
    -1/(2*h),-1/(2*w),-1/(2*h),1/(2*w),1/(2*h),1/(2*w),1/(2*h),-1/(2*w)];
H = 4;
t = 0.02;
K = t * H * det(J) * transpose(B) * D * B;

% After calculating the stiffness matrix for one element, we obtain the
% global stiffness matrix. 

k1 = K;
v1 = [1,2,3,4,5,6,7,8]; % The position vectors
k2 = K;
v2 = [7,8,5,6,11,12,9,10];
k3 = zeros(12);
for i = v1
    for j= v1
        k3(i,j) = k3(i,j) + k1(find(v1 == i),find(v1 ==j));
    end
end

for i = v2
    for j= v2
        k3(i,j) = k3(i,j) + k2(find(v2 == i),find(v2 ==j));
    end
end
K_global = k3;

u1= -0.048e-3; v1 = 0.191e-3; u2 = 0.072e-3 ; v2 = 0.183e-3; u5= 0.077e-3; v5 = 0.191e-3; u6= -0.04e-3; v6= 0.183e-3;
ks = k3;

% Having obtained the global stiffness matrix, we define the equations
% characterising the unknows, which are themselves defined as symbols.
% Then, using the 'solve' function, we find calculate their values. 

syms u3 v3 u4 v4

var = u1*ks(5,1)+v1*ks(5,2)+u2*ks(5,3)+v2*ks(5,4)+u3*ks(5,5)+v3*ks(5,6)+u4*ks(5,7)+v4*ks(5,8)+u5*ks(5,9)+v5*ks(5,10)+u6*ks(5,11)+v6*ks(5,12);
var2 = u1*ks(6,1)+v1*ks(6,2)+u2*ks(6,3)+v2*ks(6,4)+u3*ks(6,5)+v3*ks(6,6)+u4*ks(6,7)+v4*ks(6,8)+u5*ks(6,9)+v5*ks(6,10)+u6*ks(6,11)+v6*ks(6,12);
var3 = u1*ks(7,1)+v1*ks(7,2)+u2*ks(7,3)+v2*ks(7,4)+u3*ks(7,5)+v3*ks(7,6)+u4*ks(7,7)+v4*ks(7,8)+u5*ks(7,9)+v5*ks(7,10)+u6*ks(7,11)+v6*ks(7,12);
var4 = u1*ks(8,1)+v1*ks(8,2)+u2*ks(8,3)+v2*ks(8,4)+u3*ks(8,5)+v3*ks(8,6)+u4*ks(8,7)+v4*ks(8,8)+u5*ks(8,9)+v5*ks(8,10)+u6*ks(8,11)+v6*ks(8,12);

sol = solve([var == 0,var2 ==0,var3 == 0, var4==0],[u3,v3,u4,v4]);

u3 = double(sol.u3);
v3 = double(sol.v3);
u4 = double(sol.u4);
v4 = double(sol.v4);

U1 = [u1;v1;u2;v2;u3;v3;u4;v4];
U2 = [u4;v4;u3;v3;u6;v6;u5;v5];

% Having calculated the displacement values, we calculate the strain
% values, together with their interpolation functions for both the x and y
% directions. 

strains1 = B*U1;
strains2 = B*U2;

e1x = strains2(1);
e2x = strains1(1);
e1y = strains2(2);
e2y = strains1(2);
y1 = -38.5775 / 1000;
y2 = -47.0325 / 1000;

a2x = ((y1-y2)/(e1x-e2x))^-1;
a1x = e1x - a2x*y1;

a2y = ((y1-y2)/(e1y-e2y))^-1;
a1y = e1y - a2y*y1;

% Using the interpolation functions, the hand-calculated strain values are
% evaluated at the same points as the experimental strain values. The
% values are printed for easier comparison.

Y = [-34.35,-37.98,-41.60,-45.22,-48.85,-51.26];
Y = Y/1000;
Ex = a1x + a2x*Y;
Ey = a1y + a2y*Y;

ex = [-0.002, -0.0014, -0.0005, 0.0008, 0.0016, 0.002];
ey = [0.0008,0.00058,0.00021,-0.00025,-0.00053,-0.00068];

Y = [-34.35,-37.98,-41.60,-45.22,-48.85,-51.26];
Y = abs((Y+34.35))/1000;

Exx1234 = (1/w)*(1-Y(4:6)/h)*(u2-u1) + (Y(4:6)/(w*h))*(u3-u4);
Exx4365 = (Y(1:3)/(w*h))*(u3-u4) + (1/w)*(Y(1:3)/h - 1)*(u6-u5);
Exx = [rot90(rot90(Exx1234)),rot90(rot90(Exx4365))];

X = 0:0.0000000001:w;
Eyy = mean((1/h)*(1-X/w)*(v4-v1) + X/(w*h)*(v3-v2)); %1234
Eyy2 = mean((1/h)*(1-X/w)*(v5-v4) + X/(w*h)*(v6-v3)); %4365

Y = [-34.35,-37.98,-41.60,-45.22,-48.85,-51.26]/1000;

fig1=figure();
plot(Ex,Y,'ks',Exx,Y,'md',ex,Y,'r*');
hold on;
ylabel('Y coordinate (m)');
xlabel('\epsilon_x_x strain');
legend('Hand Calculation GP','Hand Calculation SF','DIC','location','northeast');
grid on;
grid minor;
fig1.PaperUnits='centimeters';
fig1.PaperPosition=[0,0,14.8,7.4];
saveas(fig1,'Image1','png');

fig2=figure();
plot(Ey,Y,'ks',ey,Y,'r*');
hold on;
ylabel('Y coordinate (m)');
xlabel('\epsilon_y_y strain');
legend('Hand Calculation','DIC','location','southeast');
grid on;
grid minor;
fig2.PaperUnits='centimeters';
fig2.PaperPosition=[0,0,14.8,7.4];
saveas(fig2,'Image2','png');
% 
% After obtaining the strain distributions using two rectangular elements
% in ANSYS, the values are imputed and compared with the expermental data
% using the graphical method.

Ex_ansys_linear = flip([3.4512e-003, 1.9596e-003, 4.721e-004, -9.7761e-004, -2.4124e-003, -3.365e-003]);
Ey_ansys_linear = flip([-1.0324e-003, -5.8487e-004, -1.3862e-004, 2.9629e-004, 7.2673e-004, 1.0125e-003]);

Ex_ansys_quad = flip([1.7429e-003,1.0015e-003,2.622e-004,-4.6414e-004, -1.186e-003, -1.6653e-003]);
Ey_ansys_quad = flip([-5.1715e-004, -2.9475e-004, -7.2959e-005, 1.4494e-004, 3.6151e-004, 5.0528e-004]);

fig3=figure();
plot(Ex_ansys_linear,Y,'bo');
hold on;
plot(Ex_ansys_quad,Y,'gd',ex,Y,'r*');
ylabel('Y coordinate (m)');
xlabel('\epsilon_x_x strain');
legend('Linear Interp.','Quad. Interp.','DIC','location','northeast');
grid on;
grid minor;
fig3.PaperUnits='centimeters';
fig3.PaperPosition=[0,0,14.8,7.4];
saveas(fig3,'Image3','png');

fig4=figure();
plot(Ey_ansys_linear,Y,'bo')
hold on;
plot(Ey_ansys_quad,Y,'gd',ey,Y,'r*');
ylabel('Y coordinate (m)');
xlabel('\epsilon_y_y strain');
legend('Linear Interp.','Quad. Interp.','DIC','location','southeast');
grid on;
grid minor;
fig4.PaperUnits='centimeters';
fig4.PaperPosition=[0,0,14.8,7.4];
saveas(fig4,'Image4','png');

% Finally, all obtained datasets are compared using the graphical, L2 and
% MSE methods. 

fig5=figure();
plot(Ex,Y,'ks');
hold on;
plot(Exx,Y,'md',Ex_ansys_linear,Y,'bo');
plot(Ex_ansys_quad,Y,'gd',ex,Y,'r*');
ylabel('Y coordinate (m)');
xlabel('\epsilon_x_x strain');
legend('Hand Calculation GP','Hand Calculation SF','Linear Interp.','Quad. Interp.','DIC','location','northeast');
grid on;
grid minor;
fig5.PaperUnits='centimeters';
fig5.PaperPosition=[0,0,14.8,7.4];
saveas(fig5,'Image5','png');

fig6=figure();
plot(Ey,Y,'ks');
hold on;
plot(Ey_ansys_linear,Y,'bo')
plot(Ey_ansys_quad,Y,'gd',ey,Y,'r*');
ylabel('Y coordinate (m)');
xlabel('\epsilon_y_y strain');
legend('Hand Calculation','Linear Interp.','Quad. Interp.','DIC','location','southeast');
grid on;
grid minor;
fig6.PaperUnits='centimeters';
fig6.PaperPosition=[0,0,14.8,7.4];
saveas(fig6,'Image6','png');


L21mx = immse(ex,Ex);
L22mx = immse(ex,Ex_ansys_linear);
L23mx = immse(ex,Ex_ansys_quad);
L24mx = immse(ex,Exx);
L21my = immse(ey,Ey);
L22my = immse(ey,Ey_ansys_linear);
L23my = immse(ey,Ey_ansys_quad);

L2mx = [L21mx;L22mx;L23mx;L24mx]
L2my = [L21my;L22my;L23my]

L21sx = sum((ex - Ex).^2);
L22sx = sum((ex - Ex_ansys_linear).^2);
L23sx = sum((ex - Ex_ansys_quad).^2);
L24sx = sum((ex - Exx).^2);
L21sy = sum((ey - Ey).^2);
L22sy = sum((ey - Ey_ansys_linear).^2);
L23sy = sum((ey - Ey_ansys_quad).^2);

L2sx = [L21sx;L22sx;L23sx;L24sx]
L2sy = [L21sy;L22sy;L23sy]