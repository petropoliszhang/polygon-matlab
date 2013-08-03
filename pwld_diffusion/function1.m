function function1()

function2(1);

a2=function2(2);
a2
return
end
%---------------------

function varargout=function2(i)

if(i<1.5)
    varargout{1} = 0;
else
    varargout{1} = 100;
end

return
end
%---------------------