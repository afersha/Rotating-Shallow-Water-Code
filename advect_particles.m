function [Xp] = advect_particles(Xp,u,v,dx,dy,dt)
    
% particle position vector for particle m: Xp(t)=[Xp(t).x_m,Xp(t).y_m]
% m = 1:np
% Inputs:  Xp:  structure with Xp.x and Xp.y
%          u, v:  velocity field
    
% advect with RK4

xp1 = dt*interpolate(Xp.x,Xp.y,u,dx,dy);
yp1 = dt*interpolate(Xp.x,Xp.y,v,dx,dy);

xp2 = dt*interpolate(Xp.x+xp1/2,Xp.y+yp1/2,u,dx,dy);
yp2 = dt*interpolate(Xp.x+xp1/2,Xp.y+yp1/2,v,dx,dy);

xp3 = dt*interpolate(Xp.x+xp2/2,Xp.y+yp2/2,u,dx,dy);
yp3 = dt*interpolate(Xp.x+xp2/2,Xp.y+yp2/2,v,dx,dy);

xp4 = dt*interpolate(Xp.x+xp3,Xp.y+yp3,u,dx,dy);
yp4 = dt*interpolate(Xp.x+xp3,Xp.y+yp3,v,dx,dy);

Xp.x = Xp.x + (xp1 + 2*xp2 + 2*xp3 + xp4)/6;
Xp.y = Xp.y + (yp1 + 2*yp2 + 2*yp3 + yp4)/6; 

%*********************************************************************

function FI = interpolate(x, y, F, dx, dy)

% Interpolate function F (velocity component u or v, on grid) to
% particle positions (x,y). Uses Lagrangian interpolation of order
% Iord (Iord = 1 ==> cubic interpolation).  See
% Duran Ch. 6, for example.  ax,ay below are x and y components
% of what he calls alpha, the fractional grid position.

% nx = 2*(kmax+1)
% dx = 2*pi/nx 

Iord = 1;
bump = 10^(-10); % Prevent NaNs

[nx,ny] = size(F);
FI = zeros(size(x));  % size of output field = # of particles

for m=1:length(x)
        
    % get x,y as fractions of nx (remove periodic wrap-around)
    xl = mod(x(m)/dx, nx);
    yl = mod(y(m)/dy, ny);
        
    % get indeces of left/bottom grid point of cell containing
    % particle, min(i0,j0) = 1,1
    i0 = 1 + floor(xl);
    j0 = 1 + floor(yl);
        
    % get fractional position within cell, 0 <= ax,ay < 1
    ax = 1 + xl - i0;
    ay = 1 + yl - j0;
    
    wx = ones(2*(Iord+1),1); wy = wx;
    for i=-Iord:Iord+1
        for j=-Iord:Iord+1
            if (i~=j) 
                wx(i+Iord+1) = wx(i+Iord+1)*(ax - j + bump)/(j - i);
                wy(i+Iord+1) = wy(i+Iord+1)*(ay - j + bump)/(j - i);
            end 
        end 
    end 
    
    for i=-Iord:Iord+1
        for j=-Iord:Iord+1
            ig = 1 + mod( i0 + i - 1, nx );
            jg = 1 + mod( j0 + j - 1, nx );
            FI(m) = FI(m) + wx(i+Iord+1)*wy(j+Iord+1)*F(ig,jg);
        end 
    end     
end
      

