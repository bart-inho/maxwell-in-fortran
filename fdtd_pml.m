%% Finite Difference Time Domain modelling on a 2-D Staggered grid with PMLs
% Numerical Modelling, ETH Zurich.
clear all;close all;clc
%% Setting up the spatial & temporal model
% --- The true dimensions (in meters)
lx = 800;
ly = 800;

% --- The true total run-time (in seconds)
lt = 1;

% --- The true medium parameters (m/s and kg/m3)
vp  = 2000;
rho = 1500;

% --- The finite-difference parameters
dx = 2;
dy = dx;
dt = 0.5e-3;

% --- The associated discrete grid in space and time
nx = floor(lx/dx);
ny = floor(ly/dy);
nt = floor(lt/dt);
x = [0:nx-1]*dx;
y = [0:ny-1]*dy;
t = [0:nt-1]*dt;

% --- Source-time wavelet (supply central frequency fc, in Hertz)
ricker = @(fm,t) (1-2*pi^2*fm^2*t.^2) .* exp(-pi^2*fm^2.*t.^2);
fc = 20;
fs = ricker(fc, t-0.1 ); % Time-delayed Ricker wavelet (plot, if desired!)

% --- Source location (positioned at 200, 300)
[~,sidx] = min(abs(x - 200));
[~,sidy] = min(abs(y - 300));

% --- Receiver location (positioned at 450, 450)
[~,ridx] = min(abs(x - 450));
[~,ridy] = min(abs(y - 450));

% --- PML parameters
npml = 60;
pmlfac = 0.85;
pmlexp = 0.85;
qx = zeros(nx,ny);
qy = zeros(nx,ny);
for a = 1:npml
    qx(a     ,:     ) = pmlfac*(npml-a)^pmlexp; % left
    qx(nx-a+1,:     ) = pmlfac*(npml-a)^pmlexp; % right
    qy(:     ,a     ) = pmlfac*(npml-a)^pmlexp; % top
    qy(:     ,ny-a+1) = pmlfac*(npml-a)^pmlexp; % bottom
end
qx = [qx(:,1) qx];
qx = [qx(1,:);qx];
qy = [qy(:,1) qy];
qy = [qy(1,:);qy];

% --- PML corners
xpml = [npml, nx-npml, nx-npml, npml, npml];
ypml = [npml, npml, ny-npml, ny-npml, npml];

%% The computational loop
% Initialize fields
px = zeros(nx+1,ny+1);
py = zeros(nx+1,ny+1);
vx = zeros(nx+1,ny+1);
vy = zeros(nx+1,ny+1);
record = zeros(nt,1);

% Stencil update locations
X = 1:nx;
Y = 1:ny;

figure(1);
for a = 1:nt
    % (1) Inject source funtion
    px(sidx,sidy) = px(sidx,sidy) + 0.5*dt*fs(a)/dx/dy*vp^2;
    py(sidx,sidy) = py(sidx,sidy) + 0.5*dt*fs(a)/dx/dy*vp^2;
    
    % (2a) Update px
    px(X+1,Y+1) = (1 - dt*qx(X+1,Y+1)) .* px(X+1,Y+1) - (dt*rho*vp^2)*((vx(X+1,Y) - vx(X,Y))/dx);
    
    % (2b) Update py
    py(X+1,Y+1) = (1 - dt*qy(X+1,Y+1)) .* py(X+1,Y+1) - (dt*rho*vp^2)*((vy(X,Y+1) - vy(X,Y))/dy); 
    
    % (2c) Combine px+py;
    p = px + py;
    
    % (3a) Update ux
    vx(X,Y) = (1 - dt*qx(X,Y)) .* vx(X,Y) - (dt/rho)*(p(X+1,Y+1) - p(X,Y+1))/dx;
    
    % (3b) Update uy
    vy(X,Y) = (1 - dt*qy(X,Y)) .* vy(X,Y) - (dt/rho)*(p(X+1,Y+1) - p(X+1,Y))/dy;
    
    % Display every 5th computation
    if (mod(a,10)==0)
        imagesc(x,y,p');
        hold on;
        plot(x(ridx), y(ridy), 'kv', 'MarkerFaceColor','k');
        plot(x(xpml), y(ypml), 'k-.');
        hold off
        axis equal;
        title(['Timestep = ', num2str(a), 'dt']);
        xlabel('x [m]');
        ylabel('y [m]');
        drawnow
    end;
    record(a) = p(ridx,ridy);
end;

% --- Compare the recorded traces
figure(2)
source_receiver_offset = norm( [x(sidx)-x(ridx), y(sidy)-y(ridy) ] );
analytical_solution = analyt_2d(fs,dt,source_receiver_offset,vp)
plot( t, record, t, analytical_solution )

% --- The root-mean-square error:
rms ( record - analytical_solution )