% Looks at the lines made by cylinder centers and constraints and compares
% them to the lines that are created by the segments in the monte carlo
% simulation. If all is well, then the points of the segments should land
% on the lines of the cylinder and constraints.
close all;
figure;
plot3(x1,y1,z1); %cylinders and constraints
hold;
plot3(VarName1, VarName2, VarName3, 'r'); %segments from monte carlo
axis([-600 600 11.58 11.59 -.1 .1]);