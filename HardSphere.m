classdef HardSphere < handle
    properties
        file_name
        run_folder
        system_size_file
        particles_file
        radii_file
        msd_file
        packing_fraction_file
        box_size
        pos
        rad
        pos_no_jump
        orig_pos
        tick
        compressionticks
        expand_rate
        step_size
        msd_offset
        pfrac
    end
    
    methods
        function obj = HardSphere(name, particle_number, polydispersity, pfrac, expand_rate, step_size)
            % Initialize the system
            obj.file_name = name;
            obj.run_folder = sprintf('%s', obj.file_name);
            if ~exist(obj.run_folder, 'dir')
                mkdir(obj.run_folder);
            end
            obj.system_size_file = fullfile(obj.file_name, 'system_size.csv');
            obj.particles_file = fullfile(obj.file_name, 'particles.csv');
            obj.radii_file = fullfile(obj.file_name, 'radii.csv');
            obj.msd_file = fullfile(obj.file_name, 'msd.csv');
            obj.packing_fraction_file = fullfile(obj.file_name, 'packing_fraction.csv');
            obj.expand_rate = expand_rate;
            obj.step_size = step_size;
            obj.pfrac = pfrac;
            obj.msd_offset = 0;
            obj.compressionticks = 0;

            if exist(obj.system_size_file, 'file') && exist(obj.particles_file, 'file') && exist(obj.radii_file, 'file')
                % Load existing data
                obj.box_size = csvread(obj.system_size_file);
                pos_data = readmatrix(obj.particles_file, 'NumHeaderLines', 0);
                pos_data = pos_data(~any(isnan(pos_data), 2), :);
                latest_tick = max(pos_data(:,1));
                obj.tick = latest_tick;
                obj.pos = pos_data(pos_data(:, 1) == latest_tick, 2:end).';
                rad_data = readmatrix(obj.radii_file, 'NumHeaderLines', 0);
                rad_data = rad_data(~any(isnan(rad_data), 2), :);
                obj.rad = rad_data(rad_data(:, 1) == latest_tick, 2:end);
                obj.pos_no_jump = pos_data(pos_data(:, 1) == latest_tick, 2:end).';
                obj.orig_pos = pos_data(pos_data(:, 1) == latest_tick, 2:end).';
                msd_data = readmatrix(obj.msd_file, 'NumHeaderLines', 0);
                obj.msd_offset = msd_data(msd_data(:, 1) == latest_tick, 2);
                disp(obj.msd_offset);
                disp('Restarting simulation from existing files...');
            else
                % Generate new particles
                [obj.pos, obj.rad, obj.box_size] = obj.generate_particles(particle_number, polydispersity, pfrac);
                obj.orig_pos = obj.pos;
                obj.pos_no_jump = obj.pos;
                obj.tick = 0;
                csvwrite(obj.system_size_file, obj.box_size);
            end
        end

        function hold(obj, tick_length, data_interval)
            % Run carlo_move and record data at intervals
            for t = 1:tick_length
                [accept] = obj.carlo_move();
                if mod(t, data_interval) == 0
                    obj.record_data();
                    disp(accept);
                end
                obj.tick = obj.tick + 1;
            end
        end

        function compress(obj, target_pfrac, shrink_interval, data_interval)
            % Compress the system
            while obj.calculate_packing_fraction() < target_pfrac
                q = abs(obj.pos - permute(obj.pos, [3 2 1]));
                o = squeeze(vecnorm(min(q, obj.box_size - q), 2, 2)) < (obj.rad + obj.rad');
                su = sum(triu(o, 1), 'all');
                if su == 0 && obj.compressionticks > shrink_interval
                    [new_shrink_rate] = obj.rescale(obj.expand_rate);
                    obj.compressionticks = 0;
                    if new_shrink_rate == obj.expand_rate
                        disp([obj.tick, "shrunk"]);
                    else
                        obj.expand_rate = new_shrink_rate;
                        disp([obj.tick, "Failed to shrink, new shrink rate:"]);
                        disp(obj.expand_rate);
                    end
                else
                    obj.compressionticks = obj.compressionticks + 1;
                    if obj.compressionticks > 1e5 || obj.tick > 1e6
                        disp([obj.tick, 'Simulation ended']);
                        return;
                    end
                end
                [accept] = obj.carlo_move();
                if mod(obj.tick, data_interval) == 0
                    obj.record_data();
                    disp(accept);
                end
                obj.tick = obj.tick + 1;
            end
        end

        function decompress(obj, target_pfrac, data_interval)
            % Decompress the system by scaling down the radii
            while obj.calculate_packing_fraction() > target_pfrac
                obj.rad = obj.rad * (1-obj.expand_rate);
                [accept] = obj.carlo_move();
                if mod(obj.tick, data_interval) == 0
                    obj.record_data();
                    disp(accept);
                end
                obj.tick = obj.tick + 1;
            end
            disp('System decompressed.');
        end

        function record_data(obj)
            % Record data to CSV files
            disp(obj.tick);
            msd = obj.calculate_msd() + obj.msd_offset;
            pfrac = obj.calculate_packing_fraction();
            writematrix([obj.tick, msd], obj.msd_file, 'WriteMode', 'append');
            writematrix([obj.tick, pfrac], obj.packing_fraction_file, 'WriteMode', 'append');
            particle_data = [[obj.tick,obj.tick,obj.tick].', obj.pos(:,:,:).'];
            writematrix(particle_data, obj.particles_file, 'WriteMode', 'append');
            writematrix([obj.tick, obj.rad], obj.radii_file, 'WriteMode', 'append');
        end

        function [poses, rad, box_size] = generate_particles(~, particle_number, polydispersity, pfrac)
            poses = [];
            rad = [];
            m = 1;
            v = (polydispersity*m)^2;
            rads = lognrnd(log((m^2)/sqrt(v+m^2)), sqrt(log(v/(m^2)+1)),particle_number,1);
            disp(rads);
            s = 0;
            for i = 1:particle_number
                s = s + 4/3*rads(i)^3*pi;
            end
            disp(s);
            height = nthroot(s/pfrac,3);
            box_size = [height,height,height];
            for i = 1:particle_number
                pos = [];
                for k = 1:1000000
                    candidate = [rand() * box_size(1) , rand() * box_size(2) , rand() * box_size(3)];
                    if ~isempty(poses)
                        if ~check_overlap_pbc(candidate,rads(i),poses,rad,box_size)
                            pos = candidate;
                            break;
                        end
                    else
                        pos = candidate;
                        break;
                    end
                end
                if isempty(pos)
                    %disp('Failed to generate a particle');
                else
                    disp(["Generated a particle at", pos])
                    poses(end+1, :) = pos;
                    rad(end + 1) = rads(i);
                end
                %scatter3(poses(:,1),poses(:,2),poses(:,3),rad.^2*pi*250);
                drawnow();
            end
            disp(["Created",length(poses)," particles"]);
        end

        function [accept] = carlo_move(obj)
            s = 0;
            for i = 1:length(obj.pos)
                p = randi(length(obj.pos));
                v = (rand(1, 3)*2 - 1) * min(obj.rad) * 0.2;
                new_pos = mod(obj.pos(p,:) + v, obj.box_size);  % Apply PBC
                if  ~check_overlap_pbc(new_pos, obj.rad(p),obj.pos([1:p-1 , p+1:end],:),obj.rad([1:p-1, p+1:end]),obj.box_size)
                    obj.pos(p,:) = new_pos;  % Update position
                    obj.pos_no_jump(p,:) = obj.pos_no_jump(p,:) + v;
                    s = s + 1;
                end
            end
            accept = s/length(obj.pos);
        end

        function [scale] = rescale(obj, scale)
            new_rad = obj.rad*(1+scale);
            %disp(count/length(pos));
            q = abs(obj.pos-permute(obj.pos,[3 2 1]));
            o = squeeze(vecnorm(min(q,obj.box_size-q),2,2)) < (new_rad + new_rad.');
            if sum(triu(o,1),'all')/length(obj.rad) < 0.025
                obj.rad = new_rad;
                return
            else
                scale = scale*.9;
                return 
            end
        end

        function msd = calculate_msd(obj)
            msd = mean(vecnorm((obj.pos_no_jump-obj.orig_pos).^2,2,2));
        end

        function packing_fraction = calculate_packing_fraction(obj)
            total_particle_area = sum(4/3*obj.rad.^3*pi);  % Calculate total area
            box_vol = obj.box_size(1)*obj.box_size(2)*obj.box_size(3);
            packing_fraction = total_particle_area / box_vol;
        end
    end
end

function overlap = check_overlap_pbc(pos1, rad1, poses, rads, box_size)
    delta = abs(poses-pos1);
    overlap = any(vecnorm(min(delta,box_size - delta),2,2) < (rad1 + rads)');
end

