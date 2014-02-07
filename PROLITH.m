%----------------------------------------------------------%
%File name: PROLITH.m
%
%Description:
%Object oriented wrapper for PROLITH lithography simulator
%----------------------------------------------------------%

classdef PROLITH
   %Used to establish and control connections to PROLITH via ActiveX server
   
   properties
       status;              %Status of ActiveX server
       server='';           %IP address of remote server to open ActiveX sever
       
       application;         %PROLITH application object
       document;            %PROLITH document object
       database;       

       sim_engine;          %PROLITH simulation object
       grid;
       waferProcessStack;
       mask;                %PROLITH mask object
       img_system;          %PROLITH imagining system object
       aberration;          %PROLITH aberration object

       final_step;
       
       focus = struct('id', 29102, 'value', 0.0);
       wavelength = struct('id', 29139, 'value', 0.0);
       threshold = struct('id', 29116, 'value', 0.0);
       mask_rotation;
   end
   
   properties (Constant)
       aerial_image = struct('contrast',14, ...
                             'cd', 17, ...
                             'placement_error', 18, ...
                             'intensity', 43);
   end
   
   methods
    function self = PROLITH(varargin)
        %class constructor
        
        %parse variable arguments:
        %server - server to run PROLITH
        if (numel(varargin)>0)
            ind = findStr('server', varargin{1});
            self.server = varargin{2}(ind);            
            self.server = char(self.server);
        end
    end
    
    function self = open(self)
        %opens ActiveX PROLITH link
        
        %set object status
        self.status = 'OPEN';
        
        %connect to PROLITH
        if ~isempty(self.server)
            self.application = actxserver('Prolith.Application', self.server.address);
        else
            self.application = actxserver('Prolith.Application');
        end
        
        self = self.getLowLevelObjects();
    end
    
    function summary = get_summary(self)
		%get_summary()
		%Returns all PROLITH settings as a string
        summary = self.sim_engine.GetParameterSummaryString(); 
    end
    
    function self = getLowLevelObjects(self)
        self.document = self.application.ActiveDocument;   
        self.database = self.application.Database;
        self.sim_engine = self.document.SimulationEngine;
        
        %instantiate image system object
        self.img_system = self.document.GetImagingSystem(1);
        
        %Load Zernike COMA_05 aberration and instantiate aberration object
        self.img_system.LoadDatabaseAberration('COMA_05', 1);
        self.aberration = self.img_system.GetAberration();   
        self.aberration.get('Coefficient', 7).set('Value', 0);
        
        self.mask = self.document.GetMask(1);
        self.grid = self.document.GetGridSize();
        
        %Remove simulation set inputs and outputs
        self.clear_inputs();
        self.clear_outputs();
    end
    
    function self = new_file(self, path)
		%new_file(path)
		%Create a new PROLITH document at the specified path.
        self.document = self.document.New();
        self = self.getLowLevelObjects();
    end
    
    function self = open_file(self, path)
		%open_file(path)
		%Open the PROLITH document at the specified path.
        self.document = self.document.get('Open', path);            
        self = self.getLowLevelObjects();
    end
    
    function self = close_file(self)
		%close_file()
		%Closes the current PROLITH document.
        self.document.Close();
        self = self.getLowLevelObjects();
    end
    
    function self = import(self,varargin)
		%import(path, overwrite)
		%Imports a file at path into the PROLITH database.
		%
		%overwrite - bool (optional: defaults to false)
		%Specifies whether or not existing database entries should be overwritten
		
        if numel(varargin) < 2
            self.database.Import(varargin{1},0);
		elseif ~varargin{2}
			self.database.Import(varargin{1},0);
        else
            self.database.Import(varargin{1},1);
        end
    end
    
    function rotate_mask(self)		
		%rotate_mask()
		%Rotate the current mask from the current state
		
        if self.mask.IsRotated
            self.mask.IsRotated = false; 
        else
            self.mask.IsRotated = true;
        end
    end
   
    function self = set_source(self, varargin)
		%Set illuminator source shape
		%
		%Shape Name - String
		%Accepted values:
		%	-conventional
		%	-dipole
		%	-annular
		%	-quadrupole
		%	-Database source
		%
		%Dipole requires 'x' or 'y' orientation to be specified
		%Quadrupole requires 45 or 90 degree orientation to be specified
		
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
    end
    
    function self = set_source_coherence(self, varargin)
		%Set the coherence of the current source.
		%
		%If source is a conventional partially coherent source, you must specify the degree of partial coherence (sigma)
		%If source is a dipole, you must specify the center and radius
		%If source is an annulus, you must supply the outer radius and inner radius		
		%If source is a quadrupole, you must supply the center and radius
		
        source = self.img_system.GetSource();
        source_type = lower(source.Name);
        
        if strcmp(source_type,'conventional--partially coherent')
            source.get('Radius').set('Value', varargin{1});
        elseif strcmp(source_type,'dipole')            
            source.get('Center').set('Value', varargin{1});
            source.get('Radius').set('Value', varargin{2});
        elseif strcmp(source_type,'annular')
			%We have to make sure that the inner radius is always less than the outer
			%so the values have to be set in the correct order
			
			inner_radius = source.get('InnerRadius');
			outer_radius = source.get('OuterRadius');
			
			if inner_radius.get('Value') > varargin{1}
				outer_radius.set('Value', varargin{1});
				inner_radius.set('Value', varargin{2});			
			else
				inner_radius.set('Value', varargin{2});							
				outer_radius.set('Value', varargin{1});
			end
        elseif strcmp(source_type,'quadrupole')
            source.get('Center').set('Value', varargin{1});
            source.get('Radius').set('Value', varargin{2});
        end
    end
    
    function set_target_cd(self, metro_plane, target)
		%set_target_cd(metro_plane, target)
		%Set the target CD for metrology plane, metro_plane, to target.
		
        self.sim_engine.SetMetrologyPlaneLithoTargetCD(metro_plane, target);
    end

    function add_aberration(self,index, value)
		%add_aberration(index, value)
        %Add Zernike aberration to imaging system  
        self.aberration.get('Coefficient', index).set('Value', value);
    end
    
    function add_set_aberration(self,index, values)
		%add_set_aberration(index, values)
        %Add Zernike aberration to simulation set
		%
		%values - vector
		%Specifies the desired start, stop, and step values for the simulation set
        ID = self.aberration.get('Coefficient', index).get('ID');
        self.add_set_input(ID, values);
    end

    function add_set_input(self, ID, values)
		%add_set_input(ID, values)
		%Add an input to a simulation set.
		%
		%ID - integer
		%PROLITH input ID for desired input
		%
		%values - vector
		%Specifies the desired start, stop, and step values for the simulation set
		
       self.sim_engine.AddInput(ID, values(1), values(2), values(3), 0);
    end
    
    function add_output(self, ID)
		%add_output(ID)
		%Add an output to a simulation set
		%
		%ID - integer
		%PROLITH output ID for desired output
		
        self.sim_engine.AddOutput(ID);
    end

    function self = clear_inputs(self)
		%clear_inputs()
		%Remove all inputs from simulation set
		
        try
            self.sim_engine.RemoveAllInputs(); 
        catch err
            if strcmp(err.identifier,'MATLAB:COM:E0'), rethrow(err); end
            self = self.close();
            self = self.open();
            self = clear_inputs(self);
        end            
    end
    
    function self = clear_outputs(self)  
		%clear_outputs()
		%Remove all outputs from simulation set
		
        try     
            self.sim_engine.RemoveAllOutputs();   
        catch err
            if strcmp(err.identifier,'MATLAB:COM:E0'), rethrow(err); end
            self = self.close();
            self = self.open();
            self = clear_outputs(self);
        end
    end
    
    function self = run_async(self)
		%run_async()
        %Run simulation asynchronously
		
        try
            self.sim_engine.SimulationRun();
        catch Exception                        
            if isempty(regexp(Exception.message, 'SimulationRun: No simulation set inputs have been specified.  AddInput must be called first.','once'))
                rethrow(Exception)
            end
            self.sim_engine.SingleRunAsynchronous();
        end
    end

    function self = run(self)
		%run()
		%Run simulation syncronously
		
        try
            self.sim_engine.RunSimSet();
        catch Exception                        
            if isempty(regexp(Exception.message, 'SimulationRun: No simulation set inputs have been specified.  AddInput must be called first.','once'))
                rethrow(Exception)
            end
            self.sim_engine.SingleRun();
        end
    end
    
    function data = get_data(self, metro_plane, ID)
		%data = get_data(metro_plane, ID)
        %Returns data from the last run simulation
		%
		%metro_plane - string
		%Name of metrology plane
		%
		%ID - integer
		%PROLITH output ID for desired output 
        
        num_sim = self.sim_engine.NumResultsRecords;

        if num_sim > 0
            data = zeros(num_sim-1,1);
            
            for i=0:num_sim-1
                data(i+1) = self.sim_engine.GetMetrologyPlaneSimSetResult(metro_plane, ID, i);
            end
        else
            data = self.sim_engine.GetMetrologyPlaneSingleRunResult(metro_plane, ID);
        end
    end
   
    function self = close(self)
       %closes PROLITH link
       
       %close file
       %self.document.Close();
       
       %set object status
       self.status='CLOSED';
       
       %release PROLITH objects
       delete(self.application);
    end        
    
    %-----------%
    %| Setters |%
    %-----------%
    function self = set.focus(self, value)
        self.focus.value = value;
        self.img_system.focus.value = value;
    end
        
    function self = set.wavelength(self, value)
        self.wavelength.value = value;
        self.img_system.wavelength.value = value;
    end
    
    function self = set.threshold(self, value)        
        self.threshold.value = value;
        invoke(self.sim_engine, 'SetInput', self.threshold.id, 0, value);
    end   
        
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

    %-----------%
    %| Getters |%
    %-----------%
    function value = get.final_step(self)
        value = self.final_step;
    end
    
    function value = get.mask_rotation(self)
       value = self.mask.IsRotated; 
    end
   end
end