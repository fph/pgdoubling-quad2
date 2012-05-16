function assertWarningThrown(f,id)
% warning test to match assertExceptionThrown in XUnit
%
% assertWarningThrown(f,id)
%
% throws an error unless f throws a warning with id "id".
% well, in fact, unless it throws *some* warnings, whose last one has id
% "id". I don't see how to make it work in the other way.

s=warning('query',id);
warning('off',id);
lastwarn('dummy warning test','dummy:warning:id');

f();

[~,id_obtained]=lastwarn;
if strcmp(id_obtained,'dummy:warning:id')
   error('matgic:assertWarningThrown:noWarning','expected warning with message_id %s, but none thrown',id); 
elseif not(strcmp(id_obtained,id))
   error('matgic:assertWarningThrown','expected warning %s, but got warning %s instead',id,id_obtained);
end

warning(s);

end
