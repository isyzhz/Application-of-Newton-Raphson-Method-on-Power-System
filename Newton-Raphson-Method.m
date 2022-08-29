% Description and derivation of matrices and values in Youtube video
% Defining the Y matrix elements
Y11 = 14;       Y12 = 10;       Y13 = 4;
Y21 = 10;       Y22 = 15;       Y23 = 5;
Y31 = 4;        Y32 = 5;        Y33 = 9;

t11 = -pi/2;    t12 = pi/2;     t13 = pi/2;
t21 = pi/2;     t22 = -pi/2;    t23 = pi/2;
t31 = pi/2;     t32 = pi/2;     t33 = -pi/2;

% Defining the known states
d1 = -0.9; v1 = 1.0;  v3 = 1.01;

% Defining the power injections
p2 = -0.9; q2 = -0.5; p3 = 1.3-0.7;
%q3 not needed for calculation (unknown)

% Defining initial values of unknowns
d2 = 0;  d3 = 0;  v2 = 1.0;  % Default initial values
X = [d2; d3; v2];
fprintf('\nDefault initial values:\n');
fprintf('Delta BB2: %f\nDelta BB3: %f\nVoltage BB2: %f\n', X(1), X(2), X(3));
m = 0;  % Iterator
dd2 = 1;    dd3 = 1;    dv2 = 1;

while (dd2 > 0.001) || (dd3 > 0.001) || (dv2 > 0.001)
    
   fp2 = Y21*v1*v2*cos(t21+d1-d2) +Y22*v2^2*cos(t22) +Y23*v3*v2*cos(t23+d3-d2) -p2;
   fp3 = Y31*v1*v3*cos(t31+d1-d2) +Y32*v2*v3*cos(t32+d2-d3) +Y33*v3^2*cos(t33) -p3;
   fq2 = -Y21*v1*v2*sin(t21+d1-d2) -Y22*v2^2*sin(t22) -Y23*v3*v2*sin(t23+d3-d2) -q2;
   fx = [fp2; fp3; fq2];
   
   J11 = Y21*v1*v2*sin(t21+d1-d2) +Y23*v3*v2*sin(t23+d3-d2);
   J12 = -Y23*v3*v2*sin(t23+d3-d2);
   J13 = Y21*v1*cos(t21+d1-d2) +2*Y22*v2*cos(t22) +Y23*v3*cos(t23+d3-d2);
   
   J21 = -Y31*v2*v3*sin(t32+d2-d3);
   J22 = Y31*v1*v3*sin(t31+d1-d3) +Y32*v2*v3*sin(t32+d2-d3);
   J23 = Y32*v3*cos(t32+d2-d3);
   
   J31 = Y21*v2*v2*cos(t21+d1-d2) +Y23*v3*v2*cos(t23+d3-d2);
   J32 = -Y23*v3*v2*cos(t23+d3-d2);
   J33 = -Y21*v1*sin(t21+d1-d2) -2*Y22*v2*sin(t22) -Y23*v3*sin(23+d3-d2);
   
   J = [J11 J12 J13; J21 J22 J23; J31 J32 J33];
   
   %NR iteration - calculation of corrections
   X = X -J\fx; 
   dd2 = d2 - X(1);     dd3 = d3 - X(2);    dv2 = v2 - X(3);
   % Assign new values to unknown variables
   d2 = X(1);   d3 = X(2);      v2 = X(3);
   m = m+1;
   
   fprintf('\nIteration no.: %d \n', m);
   fprintf('Delta BB2: %f\nDelta BB3: %f\nVoltage BB2: %f\n', X(1), X(2), X(3));
end

% Source: https://www.youtube.com/watch?v=yTSvxgkyWmw