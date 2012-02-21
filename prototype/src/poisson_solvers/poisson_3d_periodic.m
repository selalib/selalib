clear all
display('================================================================================================================================')
% Data

nx = 17;
ny = 16;
nz = 16;

dx = 2*pi/nx;
dy = 2*pi/ny;
dz = 2*pi/nz;

for k=1:nz
    z = (k-1)*dz;
    for j=1:ny
        y=(j-1)*dy;
        for i=1:nx
            x = (i-1)*dx;
            phi_an(i,j,k) = cos(x)*sin(y)*cos(z);
        end
    end
end
rho = 3*phi_an;

% Poisson solver

% FFTs in z-direction
for j=1:ny
    for i=1:nx
        hat_rho(i,j,:) = fft( rho(i,j,:) );
    end
end
for k=1:2:nz-1
    hat_rho(:,:,k) = -hat_rho(:,:,k);
end

% FFTs in y-direction
for k=1:nz
    for i=1:nx
        hat_rho(i,:,k) = fft( rho(i,:,k) );
    end
end

for j=1:2:ny-1
    hat_rho(:,j,:) = -hat_rho(:,j,:);
end

% FFTs in x-direction
for k=1:nz
    for j=1:ny
        hat_rho(:,j,k) = fft( rho(:,j,k) );
    end
end
for i=1:2:nx-1
    hat_rho(i,:,:) = -hat_rho(i,:,:);
end

hat_rho = hat_rho/(nx*ny*nz);

% Compute hat_phi, phi = inv_fft(hat_phi)
for k=-nz/2:nz/2-1
    for j=-ny/2:ny/2-1
        for i=-nx/2:nx/2-1
            if ( k==0 )
                hat_phi(i+nx/2+1,j+ny/2+1,k+nz/2+1) = 0;
            else
                hat_phi(i+nx/2+1,j+ny/2+1,k+nz/2+1) = hat_rho(i+nx/2+1,j+ny/2+1,k+nz/2+1) / ( 4*pi^2*( (i/nx)^2 + (j/ny)^2 + (k/nz)^2 ) );
            end
        end
    end
end

% Inverse FFTs in z-direction
for j=1:ny
    for i=1:nx
        phi(i,j,:) = ifft( hat_phi(i,j,:) );
    end
end
for k=1:2:nz-1
    phi(:,:,k) = -phi(:,:,k);
end

% Inverse FFTs in y-direction
for k=1:nz
    for i=1:nx
        phi(i,:,k) = ifft( hat_phi(i,:,k) );
    end
end
for j=1:2:ny-1
    phi(:,j,:) = -phi(:,j,:);
end

% Inverse FFTs in x-direction
for k=1:nz
    for j=1:ny
        phi(:,j,k) = ifft( hat_phi(:,j,k) );
    end
end
for i=1:2:nx-1
    phi(i,:,:) = -phi(i,:,:);
end

display(sum(sum(sum(abs(phi_an-phi))))/(nx*ny*nz))

