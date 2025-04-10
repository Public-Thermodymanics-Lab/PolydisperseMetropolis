% Initialize the HardSphere system
particle_number = 2000;
polydispersity = 0.0;
initial_pfrac = 0.2;
expand_rate = 0.0001;
step_size = 0.1;

% Create an instance of the HardSphere class
hs = HardSphere('Example Run',particle_number, polydispersity, initial_pfrac, expand_rate, step_size);
hs.hold(10000, 1000); %Equilibriate
hs.compress(0.5,1000,1000); % Compress to a dense packing fraction that still allows for dispersion
hs.hold(100000, 1000); %Hold again
hs.compress(0.55,1000,1000); % Compress to a denser state.
hs.hold(100000, 1000);
hs.compress(0.6,1000,1000); % Crystallization should have been reached by now
