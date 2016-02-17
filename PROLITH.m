% ----------------------------------------------------------% 
% File name: PROLITH.m
% 
% Description:
% Object oriented wrapper for PROLITH lithography simulator
% ----------------------------------------------------------% 

classdef PROLITH
   % Used to establish and control connections to PROLITH via ActiveX server
   
   properties
       status;              % Status of ActiveX server
       server='';           % IP address of remote server to open ActiveX sever
       
       application;         % PROLITH application object
       document;            % PROLITH document object
       database;       

       sim_engine;          % PROLITH simulation object
       grid;
       waferProcessStack;
       mask;                % PROLITH mask object
       img_system;          % PROLITH imagining system object
       aberration;          % PROLITH aberration object

       final_step;
       
       focus = PPI_input(29102, @get_focus);
       wavelength = PPI_input(29139, @get_wavelength);
       NA = PPI_input(29105, @get_NA);
       CRAO = PPI_input(31318, @get_CRAO);
       azimuthal_angle = PPI_input(31319, @get_azimuthal);
       reduction = PPI_input(29141, @get_reduction);
       threshold = PPI_input(29116, @() NaN);
       feature_width = PPI_input(29996, @get_feature_width);
       pitch = PPI_input(30001, @get_pitch);
       source_radius = PPI_input(29106, @() NaN);
       target_cd;
       pupil_filter;
       mask_rotation;
       metrology_planes;
       sim_region;
   end
   
   properties (Constant)
       aerial_image = struct('contrast',14, ...
                             'cd', 17, ...
                             'placement_error', 18, ...
                             'intensity', 43, ...
                             'nils', struct('right', 140, 'left', 141, 'avg', 142), ...
                             'threshold', 29116);           
        diffraction = struct('amplitude', struct('xx', 49, 'xy', 50, 'yx', 51, 'yy', 52, 'kirch', 57), ...
                             'phase', struct('xx', 53, 'xy', 54, 'yx', 55, 'yy', 56, 'kirch', 58)); 
   end
   
   methods
    function self = PROLITH(varargin)
        % class constructor
        
        % parse variable arguments:
        % server - server to run PROLITH
        if ~isempty(varargin)
            self.server = varargin{1};
        end

        self = self.open();
    end
    
    function varargout = open(self)
        % opens ActiveX PROLITH link
        
        % set object status
        self.status = 'OPEN';
        
        % connect to PROLITH
        if ~isempty(self.server)
            self.application = actxserver('Prolith.Application', self.server);
        else
            self.application = actxserver('Prolith.Application');
        end
        
        self = self.getLowLevelObjects();

        if nargout
            varargout{1} = self;
        elseif ~isempty(inputname(1))
            assignin('caller', inputname(1), self)
        end       
    end
    
    function summary = get_summary(self)
        % get_summary()
        % Returns all PROLITH settings as a string
        summary = self.sim_engine.GetParameterSummaryString(); 
    end
    
    function varargout = getLowLevelObjects(self)
        global img_system
        global sim_engine
        global mask

        self.document = self.application.ActiveDocument;   
        self.database = self.application.Database;
        self.sim_engine = self.document.SimulationEngine;
        
        % instantiate image system object
        self.img_system = self.document.GetImagingSystem(1);
        
        % Load Zernike COMA_05 aberration and instantiate aberration object
        self.img_system.LoadDatabaseAberration('COMA_05', 1);
        self.aberration = self.img_system.GetAberration();   
        self.aberration.get('Coefficient', 7).set('Value', 0);
        
        self.mask = self.document.GetMask(1);
        self.grid = self.document.GetGridSize();
        
        % Remove simulation set inputs and outputs
        self.clear_inputs();
        self.clear_outputs();

        img_system = self.img_system;
        mask = self.mask;
        sim_engine = self.sim_engine;

        if nargout
            varargout{1} = self;
        elseif ~isempty(inputname(1))
            assignin('caller', inputname(1), self)
        end
    end
    
    function varargout = new_file(self, path)
        % new_file(path)
        % Create a new PROLITH document at the specified path.
        self.document = self.document.New();
        self = self.getLowLevelObjects();

        if nargout
            varargout{1} = self;
        elseif ~isempty(inputname(1))
            assignin('caller', inputname(1), self)
        end
    end
    
    function varargout = open_file(self, path)
        % open_file(path)
        % Open the PROLITH document at the specified path.
        if ~strcmp(self.document.path, path)
            self.document = self.document.get('Open', path);            
            self = self.getLowLevelObjects();
        end

        if nargout
            varargout{1} = self;
        elseif ~isempty(inputname(1))
            assignin('caller', inputname(1), self)
        end
    end
    
    function varargout = close_file(self)
        % close_file()
        % Closes the current PROLITH document.
        self.document.Close();
        self = self.getLowLevelObjects();
        
        if nargout
            varargout{1} = self;
        elseif ~isempty(inputname(1))
            assignin('caller', inputname(1), self)
        end
    end
    
    function varargout = import(self,varargin)
        % import(path, overwrite)
        % Imports a file at path into the PROLITH database.
        % 
        % overwrite - bool (optional: defaults to false)
        % Specifies whether or not existing database entries should be overwritten
        
        if numel(varargin) < 2
            self.database.Import(varargin{1},0);
        elseif ~varargin{2}
            self.database.Import(varargin{1},0);
        else
            self.database.Import(varargin{1},1);
        end
        self = self.getLowLevelObjects();
        
        if nargout
            varargout{1} = self;
        elseif ~isempty(inputname(1))
            assignin('caller', inputname(1), self)
        end
    end
    
    function rotate_mask(self)      
        % rotate_mask()
        % Rotate the current mask from the current state
        
        if self.mask.IsRotated
            self.mask.IsRotated = false; 
        else
            self.mask.IsRotated = true;
        end
    end

    function value = best_focus(self)
        % best_focus()
        % Get best focus from simulation set analysis tab
        % Note: must run appropriate simulation set for this
        % to be valid.
        value = self.sim_engine.get('BestFocus');
    end

    function value = depth_of_focus(self)
        % depth_of_focus()
        % Get depth of focus from simulation set analysis tab
        % Note: must run appropriate simulation set for this
        % to be valid.
        value = self.sim_engine.get('DepthOfFocus');
    end
   
    function varargout = set_source(self, varargin)
        % Set illuminator source shape
        % 
        % Shape Name - String
        % Accepted values:
        %   -conventional
        %   -dipole
        %   -annular
        %   -quadrupole
        %   -Database source
        % 
        % Dipole requires 'x' or 'y' orientation to be specified
        % Quadrupole requires 45 or 90 degree orientation to be specified
        
        if strcmp(lower(varargin{1}),'conventional')
            self.img_system.LoadParametricSource(10);
        elseif strcmp(lower(varargin{1}),'dipole')
            if strcmp(lower(varargin{2}), 'y')
                self.img_system.LoadParametricSource(71);
            else
                self.img_system.LoadParametricSource(70);
            end
        elseif strcmp(lower(varargin{1}),'annular')
            self.img_system.LoadParametricSource(30);
        elseif strcmp(lower(varargin{1}),'quadrupole')
            if varargin{2}==90
                self.img_system.LoadParametricSource(40);
            elseif varargin{2}==45
                self.img_system.LoadParametricSource(41);
            end
        else
           try
                self.img_system.LoadDatabaseSource(varargin{1});
                self = self.getLowLevelObjects();
           catch exception
                disp(exception.message);
           end 
        end
        
        if nargout
            varargout{1} = self;
        elseif ~isempty(inputname(1))
            assignin('caller', inputname(1), self)
        end
    end
    
    function varargout = set_source_coherence(self, varargin)
        % Set the coherence of the current source.
        % 
        % If source is a conventional partially coherent source, you must specify the degree of partial coherence (sigma)
        % If source is a dipole, you must specify the center and radius
        % If source is an annulus, you must supply the outer radius and inner radius        
        % If source is a quadrupole, you must supply the center and radius
        
        source = self.img_system.GetSource();
        source_type = lower(source.Name);
        
        if strcmp(source_type,'conventional--partially coherent')
            source.get('Radius').set('Value', varargin{1});
        elseif strcmp(source_type,'dipole')            
            source.get('Center').set('Value', varargin{1});
            source.get('Radius').set('Value', varargin{2});
        elseif strcmp(source_type,'annular')
            % We have to make sure that the inner radius is always less than the outer
            % so the values have to be set in the correct order
            
            inner_radius = source.get('InnerRadius');
            outer_radius = source.get('OuterRadius');

            outer_radius.set('Value', 0.99999);
            inner_radius.set('Value', 0.0);
            
            % This might be broken:
            %if inner_radius.get('Value') > varargin{1}
            %    outer_radius.set('Value', varargin{1});
            %    inner_radius.set('Value', varargin{2});         
            %else
            %    inner_radius.set('Value', varargin{2});                         
            %    outer_radius.set('Value', varargin{1});
            %end

            inner_radius.set('Value', varargin{1});
            outer_radius.set('Value', varargin{2});
        elseif strcmp(source_type,'quadrupole')
            source.get('Center').set('Value', varargin{1});
            source.get('Radius').set('Value', varargin{2});
        end
        
        if nargout
            varargout{1} = self;
        elseif ~isempty(inputname(1))
            assignin('caller', inputname(1), self)
        end
    end

    function load_mask(self, mask_name, varargin)
        % load_mask(self, mask_name, varargin)
        % Load mask from PROLITH database
        %
        % Optional arguments
        % pass_num - integer 1 through 5
        % dims - mask dimensions (integer, 2 or 3)
        % load_sim_region - boolean
        % load_cse - boolean

        p = inputParser;
        addOptional(p,'pass_num',1,@(v) isinteger(v) && v>=1 && v<=5);
        addOptional(p,'dims',2,@(v) v==2 || v==3);
        addOptional(p,'load_sim_region',true, @islogical);
        addOptional(p,'load_cse',true, @islogical);

        parse(p, varargin{:});
        args = p.Results;

        self.sim_engine.LoadDatabaseMask(args.pass_num, args.dims-2, mask_name, args.load_sim_region, args.load_cse);
    end
    
    function set_target_cd(self, metro_plane, target)
        % set_target_cd(metro_plane, target)
        % Set the target CD for metrology plane, metro_plane, to target.
        
        self.sim_engine.SetMetrologyPlaneLithoTargetCD(metro_plane, target);
    end

    function add_aberration(self,index, value)
        % add_aberration(index, value)
        % Add Zernike aberration to imaging system  
        self.aberration.get('Coefficient', index).set('Value', value);
    end
    
    function add_set_aberration(self,index, values)
        % add_set_aberration(index, values)
        % Add Zernike aberration to simulation set
        % 
        % values - vector
        % Specifies the desired start, stop, and step values for the simulation set
        self.aberration.get('Coefficient', index).set('Value', 0);
        ID = self.aberration.get('Coefficient', index).get('ID');
        self.add_set_input(ID, values);
    end

    function add_input(self, sim_input, values)
        % add_set_input(sim_input, values)
        % Add an input to a simulation set.
        % 
        % sim_input - integer
        % PROLITH input for desired input
        % 
        % values - vector
        % Specifies the desired start, stop, and step values for the simulation set
        

        if strcmp(class(sim_input), 'PPI_input')
            self.sim_engine.AddInput(sim_input.id, values(1), values(2), values(3), 0);     
        else
            self.sim_engine.AddInput(sim_input, values(1), values(2), values(3), 0);
        end      
    end

    function add_set_input(self, sim_input, values)
        % This function is depricated! Please use add_input(sim_input, values) instead.
        %
        % add_set_input(sim_input, values)
        % Add an input to a simulation set.
        % 
        % sim_input - integer
        % PROLITH input for desired input
        % 
        % values - vector
        % Specifies the desired start, stop, and step values for the simulation set
        
        if strcmp(class(sim_input), 'PPI_input')
            self.add_input(sim_input.id, values)
        else
            self.add_input(sim_input, values)
        end
    end

    function id = target_cd_id(self, metro_plane)
        % target_cd_id(metro_plane)
        % Get input ID from target CD from specified metrology plane

        metro_planes = self.document.GetMetrologyPlanes();
        plane = metro_planes.GetMetrologyPlane(metro_plane);
        id = plane.LithoTargetCD.ID;
    end

    function couple_inputs(self, input_id, coupled_id)
        % couple_inputs(input_id, coupled_id)
        % Couple two simulation set inputs

        if strcmp(class(input_id), 'PPI_input'), input_id = input_id.id; end
        if strcmp(class(coupled_id), 'PPI_input'), coupled_id = coupled_id.id; end
        self.sim_engine.SetInputCoupledID(input_id, coupled_id);
    end
    
    function add_output(self, output)
        % add_output(ID)
        % Add an output to a simulation set
        % 
        % output - integer
        % PROLITH output ID for desired output
        
        self.sim_engine.AddOutput(output);
    end

    function varargout = clear_inputs(self)
        % clear_inputs()
        % Remove all inputs from simulation set
        
        try
            self.sim_engine.RemoveAllInputs(); 
        catch err
            if strcmp(err.identifier,'MATLAB:COM:E0'), rethrow(err); end
            self = self.close();
            self = self.open();
            self = clear_inputs(self);
        end
        
        if nargout
            varargout{1} = self;
        elseif ~isempty(inputname(1))
            assignin('caller', inputname(1), self)
        end   
    end
    
    function varargout = clear_outputs(self)  
        % clear_outputs()
        % Remove all outputs from simulation set
        
        try     
            self.sim_engine.RemoveAllOutputs();   
        catch err
            if strcmp(err.identifier,'MATLAB:COM:E0'), rethrow(err); end
            self = self.close();
            self = self.open();
            self = clear_outputs(self);
        end
        
        if nargout
            varargout{1} = self;
        elseif ~isempty(inputname(1))
            assignin('caller', inputname(1), self)
        end
    end
    
    function varargout = run_async(self)
        % run_async()
        % Run simulation asynchronously
        
        try
            self.sim_engine.SimulationRun();
        catch Exception                        
            if isempty(regexp(Exception.message, 'SimulationRun: No simulation set inputs have been specified.  AddInput must be called first.','once'))
                rethrow(Exception)
            end
            self.sim_engine.SingleRunAsynchronous();
        end
        
        if nargout
            varargout{1} = self;
        elseif ~isempty(inputname(1))
            assignin('caller', inputname(1), self)
        end
    end

    function varargout = run(self)
        % run()
        % Run simulation syncronously
        
        try
            self.sim_engine.RunSimSet();
        catch Exception                 
            if isempty(regexp(Exception.message, 'SimulationRun: No simulation set inputs have been specified.  AddInput must be called first.','once'))
                rethrow(Exception)
            end
            self.sim_engine.SingleRun();
        end
        
        if nargout
            varargout{1} = self;
        elseif ~isempty(inputname(1))
            assignin('caller', inputname(1), self)
        end
    end
    
    function data = get_data(self, metro_plane, output)
        % data = get_data(metro_plane, ID)
        % Returns data from the last run simulation
        % 
        % metro_plane - string
        % Name of metrology plane
        % 
        % output - integer
        % PROLITH output output for desired output 
        
        num_sim = self.sim_engine.NumResultsRecords;

        if num_sim > 0
            data = zeros(num_sim,1);
            
            try
                for i=1:num_sim
                    data(i) = self.sim_engine.GetMetrologyPlaneSimSetResult(metro_plane, output, i-1);
                end
            catch Exception
                if ~strcmp(Exception.identifier,'MATLAB:UnableToConvert')
                    rethrow(Exception)
                end

                data = {};
                for i=1:num_sim
                    data(i) = self.sim_engine.GetMetrologyPlaneSimSetResult(metro_plane, output, i-1);
                end
            end
        else
            data = self.sim_engine.GetMetrologyPlaneSingleRunResult(metro_plane, output);
        end
    end

    function data = get_diffraction(self, ID, varargin)
        % data = get_diffraction(ID, pass*)
        % Returns the diffraction pattern data
        % 
        % ID - integer
        % PROLITH output ID for desired output 
        %
        % pass - integer (optional)
        % Exposure pass

        if ~isempty(varargin)
            pass = varargin{1};
        else
            pass = 1;
        end

        num_sim = self.sim_engine.NumResultsRecords;

        if num_sim > 0
            data = cell(num_sim,1);
            
            for i=1:num_sim
                val = invoke(self.sim_engine,'GetDiffractionPatternSimSetResult', ID, i-1);
                data{i} = val;
            end
        else
            data = self.sim_engine.GetDiffractionPatternSingleRunResult(pass, ID);
        end

    end
   
    function varargout = close(self)
       % closes PROLITH link
       
       % set object status
       self.status='CLOSED';
       
       % release PROLITH objects
       delete(self.application);
        
       if nargout
            varargout{1} = self;
        elseif ~isempty(inputname(1))
            assignin('caller', inputname(1), self)
        end
    end

    function generate_filter(self, x, amplitude, phase)
        %  Write PROLITH .FIL file from arrays amplitude and phase.
        %  x should be a pupil coordinate vector.

        %  Create pupil coordinate grid and convert
        %  phase function to degrees

        persistent A
        persistent P

        if ~isempty(A) && ~isempty(P) && isequal(A, amplitude) & isequal(P, phase)
            return
        end

        [X,Y] = meshgrid(x,x);
        A = amplitude;
        P = rad2deg(phase);

        %  Add the header
        rows = {};
        rows{1} = '[Version]';
        rows{2} = '15.0.1.15';
        rows{3} = '';
        rows{4} = '[Parameters]';
        rows{5} = 'Wavefront ;Pupil Filter Name';
        rows{6} = '1 ;0 = Radial Data, 1 = x,y Grid Data';
        rows{7} = sprintf('% 0.9f ;Step Size (x-y grid only)', x(2)-x(1));
        rows{8} = '[Data]';

        %  Transform gridded arrays into vectors
        X = reshape(X,[],1);
        Y = reshape(Y,[],1);
        A = reshape(A,[],1);
        P = reshape(P,[],1);

        %  Remove any points where both the amplitude function
        %  and the phase function are equal to zero
        X = X(A~=0 | P~=0);
        Y = Y(A~=0 | P~=0);
        temp = A(A~=0 | P~=0);
        P = P(A~=0 | P~=0);
        A = temp; 

        %  Add each coordinate to the file
        r=9;
        for i=1:numel(X)
            if ~isnan(P(i)) & ~isnan(A(i))
                rows{r} = sprintf('% 0.3f\t% 0.3f\t% 0.3f\t% 0.3f', X(i), Y(i), P(i), A(i));
                r = r+1;
            end
        end

        %  Write out the file
        filePath = fullfile(pwd, 'wavefront.FIL');
        file  = fopen(filePath, 'w');
        for r=1:numel(rows)
            fprintf(file,'% s\n',rows{r});
        end
        fclose(file);
    end        
    
    % -----------% 
    % | Setters |% 
    % -----------% 
    function self = set.pupil_filter(self, value)
        if ischar(value) && ~isempty(value)
            self.img_system.LoadDatabasePupilFilter(value);
            self.pupil_filter = value;
        else
            self.img_system.UnloadPupilFilter();
            self.pupil_filter = '';
        end
    end

    function self = set.focus(self, value)
        self.img_system.focus.value = value;
    end
        
    function self = set.wavelength(self, value)
        self.img_system.wavelength.value = value;
    end
    
    function self = set.CRAO(self, value)
        self.img_system.ChiefRayIncidentAngle.value = value;
    end

    function self = set.azimuthal_angle(self, value)
        self.img_system.ChiefRayAzimuthalAngle.value = value;
    end
    
    function self = set.NA(self, value)
        self.img_system.NumericalAperture.value = value;
    end
    
    function self = set.reduction(self, value)
        self.img_system.ReductionRatio.value = value;
    end
    
    %function self = set.threshold(self, value)        
    %    self.threshold.value = value;
    %    invoke(self.sim_engine, 'SetInput', self.threshold.id, 0, value);
    %end   
        
    function self = set.mask_rotation(self,value)
        if islogical(value)
            self.mask_rotation = value;
            self.mask.IsRotated = value;
        end
    end
    
    function self = set.final_step(self, value)
        if strcmp(lower(value),'aerial image')
            step_num = 3;
        elseif strcmp(lower(value),'reflectivity')
            step_num = 4;
        elseif strcmp(lower(value),'latent image')
            step_num = 5;
        elseif strcmp(lower(value),'PEB')
            step_num = 6;
        elseif strcmp(lower(value),'develop contours')
            step_num = 7;
        elseif strcmp(lower(value),'resist profile')
            step_num = 8;
        elseif strcmp(lower(value),'diffraction pattern')
            step_num = 12;
        elseif strcmp(lower(value),'image in resist')
            step_num = 13;
        elseif strcmp(lower(value),'mask')
            step_num = 14;
        elseif strcmp(lower(value),'etch')
            step_num = 16;
        elseif strcmp(lower(value),'topography')
            step_num = 18;
        end

        self.sim_engine.SetFinalState(step_num);
    end

    function self = set.feature_width(self, value)
        self.mask.FeatureWidth = value;
    end

    function self = set.pitch(self, value)
        self.mask.Pitch = value;
    end

    function self = set.target_cd(self, value)        
        metro_planes = self.document.GetMetrologyPlanes();
        plane = metro_planes.GetMetrologyPlane(value{1});
        plane.LithoTargetCD.Value = value{2};
    end

    function self = set.sim_region(self, value)
        invoke(self.sim_engine, 'SimulationRegion', value);
    end

    % -----------% 
    % | Getters |% 
    % -----------% 
    function value = get.final_step(self)
        value = self.final_step;
    end
    
    function value = get.mask_rotation(self)
       value = self.mask.IsRotated; 
    end
        
    function value = get.reduction(self)
       value = self.img_system.ReductionRatio.value; 
    end   
    
    function value = get.metrology_planes(self)
        value = self.sim_engine.GetAllMetrologyPlanes();
    end

    function value = get.sim_region(self)    
        value = self.sim_engine.SimulationRegion;
    end
   end
end

function value = get_focus()
    global img_system
    value = img_system.focus.value;
end

function value = get_wavelength()
    global img_system
    value = img_system.wavelength.value;
end

function value = get_NA()
    global img_system
    value = img_system.NumericalAperture.value;
end

function value = get_CRAO()
    global img_system
    value = img_system.ChiefRayIncidentAngle.value;
end

function value = get_azimuthal()
    global img_system
    value = img_system.ChiefRayAzimuthalAngle.value;
end

function value = get_reduction()
    global img_system
    value = img_system.ReductionRatio.value;
end

function value = get_feature_width()
    global mask
    value = mask.FeatureWidth.value;
end

function value = get_pitch()
    global mask
    value = mask.Pitch.value;
end