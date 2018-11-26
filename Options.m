classdef Options < handle
    % Also known as: I wish MATLAB had named arguments like Python
    %
    %
    % (c) 2011-2012 F. Poloni <poloni@math.tu-berlin.de> and others 
    % see AUTHORS.txt and COPYING.txt for details
    % https://bitbucket.org/fph/pgdoubling

    properties (Access=private)
        values_ %values.key = something
        used_  %used.key = boolean
    end
    methods
        
        function o=Options(varargin)
            %class constructor for Options
            %
            % Options('key1',val1,'key2',val2,...) or
            % Options(previousOptions,'keyN',valN,...)
            %
            % Construct an Options structure that can be used and passed
            % around.
            %

            o.values_=struct();
            o.used_=struct();
            
            if nargin==0
                return
            end
            index=1;            
            if isa(varargin{index},'Options')
                o=varargin{index};
                index=index+1;
            else
                o.values_=struct(varargin{index:end});
                
                %populates o.used_
                keys=fieldnames(o.values_);
                for i=1:length(keys)
                    o.used_.(keys{i})=false;
                end 
            end
        end
        
        function set(o,key,value)
            % sets an option
            %
            % set(o,key,value);
            
            o.values_.(key)=value;
            o.used_.(key)=false;
        end
            
        function b=isSet(o,key)
            % returns true if an option is set
            %
            % bool = o.isSet('key')
            
            b=isfield(o.values_,key);
        end
        
        function v=get(o,key,varargin)
            % returns an option registered in an Options structure
            %
            % value=o.get('key',defaultValue)
            % returns option named 'key', if it exists, or defaultValue
            % otherwise
            %
            % value o.get('key') returns option named 'key', if it exists,
            % otherwise throws an error
            %
            % see also: Options.look
            
            v=o.look(key,varargin{:});
            o.used_.(key)=true;
        end
        
        function v=look(o,key,default)
            % returns an option without registering it as "used"
            %
            % works like Option.get, but does not mark the key as "used" --
            % you will still get a warning in Options.delete if you don't
            % get() it.
            %
            % see also: Options.get, Options.delete
            if(not(ischar(key)))
                error('Options:keyMustBeAString','key name is of type "%s" instead of string',class(key));
            end
            if(isfield(o.values_,key))
                v=o.values_.(key);
            else 
                if(nargin<=2)
                    error('Options:noSuchKey','key %s does not exist in the Options and no default value specified',key);
                else
                    v=default;
                end
            end            
        end
        
        function delete(o)
            % Options destructor
            %
            % apart from freeing structures, returns a warning if some
            % options were defined but never used (probably a typo)

            keys=fieldnames(o.values_);
            for i=1:length(keys)
                if o.used_.(keys{i})==false
                    warning('Options:unusedKey','option %s was defined but never used. Maybe you mis-spelled it?', keys{i});
                end
            end
            
        end
        
    end
end
