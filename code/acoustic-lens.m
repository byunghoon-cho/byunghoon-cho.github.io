% Requires k-wave library

freq = 5e3;                     % Source frequency [Hz]
source_to_slit = 0.2;           % Distance from source to lens [m]

Nx = 256;                       % Number of grid points in x (longitudinal) direction
Ny = 160;                       % Number of grid points in y direction
Nz = 160;                       % Number of grid points in z direction
dx = source_to_slit * 4 / Nx;   % Grid point spacing in the x direction [m]
dy = dx;                        % Grid point spacing in the y direction [m]
dz = dx;                        % Grid point spacing in the z direction [m]

kgrid = makeGrid(Nx, dx, Ny, dy, Nz, dz);

% Properties of medium (air)
c0 = 343;                       % Speed of sound [m/s]
rho0 = 1.225;                   % Density [kg/m^3]
medium.sound_speed = c0 * ones(Nx, Ny, Nz); % [m/s]
medium.density = rho0 * ones(Nx, Ny, Nz);   % [kg/m^3]
medium.alpha_power = 1.5;
medium.alpha_coeff = (3e-3) / (freq/1e6)^medium.alpha_power; % [dB/(MHz^y cm)]

t_end = dx * Nx / c0;           % [s]
CFL = 0.02;
[kgrid.t_array, dt] = makeTime(kgrid, medium.sound_speed, CFL, t_end);
Nt = length(kgrid.t_array);

% Point Source
source.p_mask = zeros(Nx, Ny, Nz);
source.p_mask(end - Nx / 4, Ny / 2, Nz / 2) = 1;

% Time varying sinusoidal wave
source_freq = freq;             % [Hz]
source_mag = 10;                % [Pa]
source.p = source_mags * sin(2 * pi * source_freq * kgrid.t_array);
source.p = filterTimeSeries(kgrid, medium, source.p);
sensor.mask = ones(Nx, Ny, Nz);

% Barrier/Lens
barrier_thickness = 6;          % Lens thickness [mm]
grid_thickness = round(barrier_thickness / 1e3 / dx); % [m]
radii_mm = [41.63, 72.87, 95.05, 113.6, 130.09, 145.22, 159.38, 172.8, 185.65]; % [mm]
radii_grids = round(radii_mm / 1e3 / dy);
barrier_mask = zeros(Nx, Ny, Nz);
disc1 = reshape(makeDisc(Ny, Nz, Ny / 2, Nz / 2, radii_grids(1)), 1, Ny, Nz);
disc2 = reshape(makeDisc(Ny, Nz, Ny / 2, Nz / 2, radii_grids(3)) & ~makeDisc(Ny, Nz, Ny / 2, Nz / 2, radii_grids(2)), 1, Ny, Nz);
disc3 = reshape(makeDisc(Ny, Nz, Ny / 2, Nz / 2, radii_grids(5)) & ~makeDisc(Ny, Nz, Ny / 2, Nz / 2, radii_grids(4)), 1, Ny, Nz);
disc4 = reshape(makeDisc(Ny, Nz, Ny / 2, Nz / 2, radii_grids(7)) & ~makeDisc(Ny, Nz, Ny / 2, Nz / 2, radii_grids(6)), 1, Ny, Nz);
disc5 = reshape(makeDisc(Ny, Nz, Ny / 2, Nz / 2, radii_grids(9)) & ~makeDisc(Ny, Nz, Ny / 2, Nz / 2, radii_grids(8)), 1, Ny, Nz);
template = disc1 | disc2 | disc3 | disc4 | disc5;
barrier_mask(Nx / 2:(Nx / 2 + grid_thickness - 1), :, :) = ~repmat(template, grid_thickness, 1, 1);

% Control (no barrier)
if ~control
    barrier_sound_speed = 3500;
    barrier_density = 700;
    medium.sound_speed(barrier_mask == 1) = barrier_sound_speed;
    medium.density(barrier_mask == 1) = barrier_density;
end

% Acoustic parameters to record
sensor.record = {'p_final', 'p_max_all', 'p_rms'};
sensor.record_start_index = floor(0.8 * Nt);

input_args = {'DisplayMask', barrier_mask | source.p_mask, 'DataCast', 'single'};
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

if ~control
    save('3D_request_05_withbarrier.mat', '-v7.3');
else
    save('3D_request_05_nobarrier.mat', '-v7.3');
end
