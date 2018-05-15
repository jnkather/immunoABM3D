%CLASS_INTERFACE Example MATLAB class wrapper to an underlying C++ class
classdef SIMengine_interface < handle
    properties (SetAccess = public, Hidden = true)
        objectHandle; % Handle to the underlying C++ class instance
    end
    methods
        %% Constructor - Create a new C++ class instance
        function this = SIMengine_interface(varargin)
            this.objectHandle = SIMengine_interface_mex('new', varargin{:});
        end
        
        %% Destructor - Destroy the C++ class instance
        function delete(this)
            SIMengine_interface_mex('delete', this.objectHandle);
        end
        
        function varargout = initialize(this, varargin)
            [varargout{1:nargout}] = SIMengine_interface_mex('initialize', this.objectHandle, varargin{:});
        end
        
        function varargout = initializeFromState(this, varargin)
            [varargout{1:nargout}] = SIMengine_interface_mex('initializeFromState', this.objectHandle, varargin{:});
        end
        
        
        function varargout = getState(this, varargin)
            [varargout{1:nargout}] = SIMengine_interface_mex('getState', this.objectHandle, varargin{:});
        end
        
        function varargout = readState(this, varargin)
            [varargout{1:nargout}] = SIMengine_interface_mex('readState', this.objectHandle, varargin{:});
        end
        
        function TU_go_grow_die(this, varargin)
            SIMengine_interface_mex('TU_go_grow_die', this.objectHandle, varargin{:});
        end
        
        function modulateAdjuvanticity(this, varargin)
            SIMengine_interface_mex('modulateAdjuvanticity', this.objectHandle, varargin{:});
        end
        
        function decayAdjuvanicityMap(this, varargin)
            SIMengine_interface_mex('decayAdjuvanicityMap', this.objectHandle, varargin{:});
        end
        
        function IMinflux(this, varargin)
            SIMengine_interface_mex('IMinflux', this.objectHandle, varargin{:});
        end
        
        function MPinflux(this, varargin)
            SIMengine_interface_mex('MPinflux', this.objectHandle, varargin{:});
        end
        
        function lymphocytesAct(this, varargin)
            SIMengine_interface_mex('lymphocytesAct', this.objectHandle, varargin{:});
        end
        
        function macrophagesAct(this, varargin)
            SIMengine_interface_mex('macrophagesAct', this.objectHandle, varargin{:});
        end
        
        function updateNecroMap(this, varargin)
            SIMengine_interface_mex('updateNecroMap', this.objectHandle, varargin{:});
        end
        
        function updateChemoMap(this, varargin)
            SIMengine_interface_mex('updateChemoMap', this.objectHandle, varargin{:});
        end
        
        function seedFibrosis(this, varargin)
            SIMengine_interface_mex('seedFibrosis', this.objectHandle, varargin{:});
        end
		
		function varargout = TUcellsNum(this, varargin)
            [varargout{1:nargout}] = SIMengine_interface_mex('TUcellsNum', this.objectHandle, varargin{:});
        end
        
    end
end