function [opt, isdefault]= set_defaults(opt, varargin)
%[opt, isdefault]= set_defaults(opt, defopt)
%[opt, isdefault]= set_defaults(opt, field/value list)
%
% This functions fills in the given struct opt some new fields with
% default values, but only when these fields DO NOT exist before in opt.
% Existing fields are kept with their original values.
% There are two forms in which you can can specify the default values,
% (1) as struct, 
%   opt= set_defaults(opt, struct('color','g', 'linewidth',3));
%
% (2) as property/value list, e.g.,
%   opt= set_defaults(opt, 'color','g', 'linewidth',3);
%
% The second output argument isdefault is a struct with the same fields
% as the returned opt, where each field has a boolean value indicating
% whether or not the default value was inserted in opt for that field.
%
% The default values should be given for ALL VALID property names, i.e. the
% set of fields in 'opt' should be a subset of 'defopt' or the field/value
% list. A warning will be issued for all fields in 'opt' that are not present
% in 'defopt', thus possibly avoiding a silent setting of options that are
% not understood by the receiving functions. 
%
% $Id$
% 
% Copyright (C) Fraunhofer FIRST
% Authors: Frank Meinecke (meinecke@first.fhg.de)
%          Benjamin Blankertz (blanker@first.fhg.de)
%          Pavel Laskov (laskov@first.fhg.de)

if length(opt)>1,
  error('first argument must be a 1x1 struct');
end

% Set 'isdefault' to ones for the field already present in 'opt'
isdefault= [];
if ~isempty(opt),
  for Fld=fieldnames(opt)',
    isdefault= setfield(isdefault, Fld{1}, 0);
  end
end

% Check if we have a  field/value list
if length(varargin) > 1
  
  % If the target is a propertylist structure use propertylist2struct to
  % convert the property list to a defopt structure.
  if (ispropertystruct(opt))
    defopt = propertylist2struct(varargin{:});
      
  else  % otherwise construct defopt from scratch
    
    
    % Create a dummy defopt structure: a terrible Matlab hack to overcome
    % impossibility of incremental update of an empty structure.
    defopt = struct('matlabsucks','foo');
  
    % Check consistency of a field/value list: even number of arguments
    nArgs= length(varargin)/2;
    if nArgs~=round(nArgs) & length(varargin~=1),
      error('inconsistent field/value list');
    end
    
    % Write a temporary defopt structure
    for ii= 1:nArgs,
      defopt= setfield(defopt, varargin{ii*2-1}, varargin{ii*2});
    end
    
    % Remove the dummy field from defopt
    defopt = rmfield(defopt,'matlabsucks');
  end
  
else  
  
  % If varargin has only one element, it must be a defopt structure.
  defopt = varargin{1};
  
end
  
% Replace the missing fields in 'opt' from their 'defopt' counterparts. 
for Fld=fieldnames(defopt)',
  fld= Fld{1};
  if ~isfield(opt, fld),
    opt= setfield(opt, fld, getfield(defopt, fld));
    isdefault= setfield(isdefault, fld, 1);
  end
end

% Check if some fields in 'opt' are missing in 'defopt': possibly wrong
% options.
for Fld=fieldnames(opt)',
  fld= Fld{1};
  if ~isfield(defopt,fld)
%    warning('set_defaults:DEFAULT_FLD',['field ''' fld ''' does not have a valid default option']);
  end
end
