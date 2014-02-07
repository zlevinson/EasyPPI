%--------------------------------%
%File name: PROLITH.m
%Author: Zac Levinson
%Email: zal2186@rit.edu
%
%Description:
%Object oriented interface with
%PROLITH lithography simulator
%--------------------------------%

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
        self.document = self.document.New();
        self = self.getLowLevelObjects();
    end
    
    function self = open_file(self, path)
        self.document = self.document.get('Open', path);            
        self = self.getLowLevelObjects();
    end
    
    function self = close_file(self)
        self.document.Close();
        self = self.getLowLevelObjects();
    end
    
    function self = import(self,varargin)
        if numel(varargin) < 2
            self.database.Import(varargin{1},0);
        else
            self.database.Import(varargin{1},1);
        end
    end
    
    function rotate_mask(self)
        if self.mask.IsRotated
            self.mask.IsRotated = false; 
        else
            self.mask.IsRotated = true;
        end
    end
   
    function self = set_source(self, varargin)
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
        source = self.img_system.GetSource();
        source_type = lower(source.Name);
        
        if strcmp(source_type,'conventional--partially coherent')
            source.get('Radius').set('Value', varargin{1});
        elseif strcmp(source_type,'dipole')            
            source.get('Center').set('Value', varargin{1});
            source.get('Radius').set('Value', varargin{2});
        elseif strcmp(source_type,'annular')
            source.get('OuterRadius').set('Value', varargin{1});
            source.get('InnerRadius').set('Value', varargin{2});
        elseif strcmp(source_type,'quadrupole')
            source.get('Center').set('Value', varargin{1});
            source.get('Radius').set('Value', varargin{2});
        end
    end
    
    function set_target_cd(self, metro_plane, target)
        self.sim_engine.SetMetrologyPlaneLithoTargetCD(metro_plane, target);
    end

    function add_aberration(self,index, value)
        %Add aberration to imaging system  
        self.aberration.get('Coefficient', index).set('Value', value);
    end
    
    function add_set_aberration(self,index, values)
        %Add aberration to simulation set
        ID = self.aberration.get('Coefficient', index).get('ID');
        self.add_set_input(ID, values);
    end

    function add_set_input(self, ID, values)
       self.sim_engine.AddInput(ID, values(1), values(2), values(3), 0);
    end
    
    function add_output(self, ID)
        self.sim_engine.AddOutput(ID);
    end

    function self = clear_inputs(self)
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
        try
            self.sim_engine.RunSimSet();
        catch Exception                        
            if isempty(regexp(Exception.message, 'SimulationRun: No simulation set inputs have been specified.  AddInput must be called first.','once'))
                rethrow(Exception)
            end
            self.sim_engine.SingleRun();
        end
    end
    
    function data = get_data(self, metroPlane, output)
        %gets data from PROLITH simulation
        
        numSim = self.sim_engine.NumResultsRecords;

        if numSim > 0
            data = zeros(numSim-1,1);
            
            for i=0:numSim-1
                data(i+1) = self.sim_engine.GetMetrologyPlaneSimSetResult(metroPlane, output, i);
            end
        else
            data = self.sim_engine.GetMetrologyPlaneSingleRunResult(metroPlane, output);
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