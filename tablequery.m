function [ newdatatable, criteriaStrings ] = tablequery( DataTable, varargin )
%TABLEQUERY Query from a matlab DataTable, using match & range criteria
%   TABLEQUERY operates as a logical AND query to a table, i.e.
%        extracts the part of DataTable that satisfies ALL query criteria
%   TABLEQUERY extracts data from a table based on combinations of
%        Variable Names + Query Criteria
%   TABLEQUERY supports two methods of specifying query criteria:
%     1) Using separate input cell arrays for 'match' and/or 'range'
%         w/cells consisting of paired {VarName,QryCriteria,VN2,QC2,...}
%     2) Using an input structure with 'match' and/or 'range' fields
%
%   SYNTAX
%
%   TABLEQUERY(DataTable)
%   TABLEQUERY(... 'match',matchCriteriaCell)
%   TABLEQUERY(... 'range',rangeCriteriaCell)
%   TABLEQUERY(... 'qryStruct',qryCriteriaStructure)
%   TABLEQUERY(... 'varnameq',{'din','phy'})
%   TABLEQUERY(... 'suppressoutput',false)
%
%     INPUT (REQUIRED):
%         - DataTable - The table from which to query a portion.
%                Each var. column in DataTable has only one data type
%                e.g. replace missing value chars w/NaN in numeric columns
%
%     INPUT (OPTIONAL VARARGIN PARAMETER/VALUE PAIRS):
%         - 'match' - a cell in the form:
%                {'qryVarName', qryVarCriteria, 'qryVarName2', etc.}
%         - 'range' - a cell in the form:
%                {'qryVarName', qryVarCriteria, 'qryVarName2', etc.}
%                NOTE: range is an INCLUSIVE range, & must contain two-
%                         element numeric vectors in a cell, such as
%                         {[ rangeMin rangeMax ]}
%         - 'qryStruct' - an alternative method of specifying 
%                         query criteria. A structure that contains two
%                         fields, each of which are structures themselves
%                Default form: struct('match',struct(),'range',struct())
%                NOTE: if match/range fields are empty then they are unused
%                      each criteria in the structure must be within a cell
%         - 'varnameq' / default: {'all'}
%                either - a char name of one of the table variables or
%                           - a cell containing multiple variable names or
%                           - a cell containing the string 'all'
%         - 'suppressoutput' / default: false
%                if set to true, then this function will  display
%                progress and matching in the MATLAB command window
%               
%     USAGE EXAMPLES:
%         USAGE EXAMPLE 1:
%                 (using a 'qryStruct' parameter)
%            >> varsToRetrieve = {'site','t','Concentration','z'}
%            >> qCriteria.match.site = {'j171'}
%            >> qCriteria.range.t = {[302.5 307]}
%            >> qCriterua.range.z = {[100 150]} 
%
%            >> tablequery(DataTable, qCriteria, 'varnameq',varsToRetrieve)
%           
%         USAGE EXAMPLE 2:
%                 (using 'match' & 'range' parameters)
%            >> Mqry = {'site', 'j171'}
%            >> Rqry = {'time',[305 320], 'z',[100 150]}
%
%            >> tablequery(DataTable, 'match', Mqry, 'range', Rqry);
%
%         USAGE EXAMPLE 3:
%                 (using 'match' & 'range' parameters w/empty criteria)
%            >> Mqry = {'time',[307.5 309.5], 'depth',40}
%            >> Rqry = {'',''}  % <<- signifies no criteria; could
%                                       alternatively use {} or {'',[]'}
%
%            >> tablequery(DataTable, 'match', Mqry, 'range', Rqry);
%
%   MFILE:   tablequery.m
%   MATLAB:  8.4.0.150421 (R2014b)
%   AUTHOR:  Daniel Edward Kaufman (USA)
%            @ The Virginia Institute of Marine Science
%   CONTACT: dkauf42@gmail.com
%   REVISION HISTORY:
%   - Bug-fix (Feb, 2016)
%   - Added support for matches and ranges of class datetime. (Nov, 2015)
%   - Added transfer of Data Table Units and Descriptions. (Apr, 2015)
%   - Removed extraneous code and reordered inputs. (Mar, 2015)
%   - Updated input options, increased speed and bug-fixes. (Feb, 2015)
%   - Initial generation. (Jan, 2015)
%

blankQryStruct = struct('match',struct(),'range',struct());
blankMatch     = {'',''};
blankRange     = {'',[]};
%% Initial Function Housekeeping 
% Parse Function Input
    p=inputParser;
    addRequired(p,'DataTable',@istable);
    addParameter(p,'match',{'',''},@iscell);
    addParameter(p,'range',{'',''},@iscell);
    addParameter(p,'qryStruct',blankQryStruct,@chkQCriteria);
    addParameter(p,'varnameq',{'all'},@iscellstr);
    addParameter(p,'suppressoutput',false,@islogical);
parse(p,DataTable, varargin{:})
matchQ=p.Results.match;
rangeQ=p.Results.range;
qryStruct=p.Results.qryStruct;
varnameq=p.Results.varnameq;
sprsoutput=p.Results.suppressoutput;

% Additional Input Checks
passedStruct = ~isequal(qryStruct,blankQryStruct);
passedMatch  = ~isequal(matchQ,blankMatch);
passedRange  = ~isequal(rangeQ,blankRange);
if passedStruct && (passedMatch || passedRange)
    error(sprintf(['Ambiguous, overconstrained input.\n ',...
            'Please specify query criteria as either ',...
            '(''match'' &/or ''range'') ',...
            'cell arrays or a (''qryStruct'') structure, not both.']));
end
if isempty(varnameq); error('varnameq input cannot be empty'); end;
if ~sprsoutput; fprintf('tablequery processing query criteria... \n'); end;
clear 'blankQryStruct' 'blankMatch' 'blankRange'

% Get Info About Queried Variable(s)
allVarNames = DataTable.Properties.VariableNames;
numOfTableVars = numel(allVarNames);
initialdatalength = size(DataTable,1);
if numel(varnameq)==1 && strcmpi(varnameq{1},'all')
    varnameq = allVarNames;
end
numOfVarQ = numel(varnameq);

%% Parse 'match' and/or 'range' input
if passedMatch
    qryStruct = parseMorRcriteriaCell(matchQ,qryStruct,'match');
end
if passedRange
    qryStruct = parseMorRcriteriaCell(rangeQ,qryStruct,'range');
end

%% Investigate/Validate the Query Criteria
ctns=fieldnames(qryStruct); % ctns stands for Criteria Type Names, i.e. 'match' and 'range'
types = struct();
criteriaStrings = cell(0); critStrii=0;
for cii = 1:numel(ctns) % CHECK BOTH TYPES OF CRITERIA (MATCH and RANGE)
    mORr = ctns{cii}; % mORr = Match or Range
    qcnames.(mORr)= fieldnames(qryStruct.(mORr));
    qcn.(mORr) = numel(qcnames.(mORr)); % number of query criteria within match and/or range
    
    types.(mORr) = struct();
    for iiqc=1:qcn.(mORr)
        thisCrit = qcnames.(mORr){iiqc};
        % Each query criterion must have a corresponding queried variable
        % Ensure qCriteria variables match at least one queried variable(s)
        validatestring( thisCrit, varnameq );
        % Get the matlab type of each match criterion
        types.(mORr).(thisCrit)=class(qryStruct.(mORr).(thisCrit){1});
        % Ensure queried vars. have corresponding vars. in DataTable
        validatestring( thisCrit, allVarNames );
        if strcmp(mORr,'range')
            if ~isnumeric(qryStruct.(mORr).(thisCrit){1}) && ~isdatetime(qryStruct.(mORr).(thisCrit){1})
                error('Range Criteria Must be Numeric')
            elseif numel(qryStruct.(mORr).(thisCrit){1}) ~= 2 && ~isempty(qryStruct.(mORr).(thisCrit){1})
                error('Range Criterion Must Be A Two Element Numeric Vector')
            end
        end
    end
    
    % Display Query Criteria in Command Window
    if ~sprsoutput
        if qcn.(mORr) > 0
            fprintf(['-- Query (', mORr, ') Criteria = \n']);
            for iiqc=1:qcn.(mORr)
                critStrii = critStrii+1;
                switch types.(mORr).(qcnames.(mORr){iiqc})
                    case 'double'
                        css = num2str( qryStruct.(mORr).(qcnames.(mORr){iiqc}){:} );
                    case 'cell'
                        error('Criteria themselves should not be in cell format');
                    case 'datetime'
                        css = datestr( qryStruct.(mORr).(qcnames.(mORr){iiqc}){:} );
                        if ~isvector(css)
                            css = css.';
                            css(end+1,:) = repmat(' ',size(css,2),1);
                            css = css(:).';
                        end
                    otherwise
                        css = num2str( qryStruct.(mORr).(qcnames.(mORr){iiqc}){:} );
                end
                critStr = [' ',qcnames.(mORr){iiqc},': ',css];
                fprintf(critStr);
                fprintf('\n');
                criteriaStrings{critStrii} = critStr;
                critStrii = critStrii+1;
            end
        else
            fprintf(['-- No (', mORr, ') Query Criteria \n']);
        end
    end
end

%% Part I : Identify data section of (Match) Criteria
matches = ones(initialdatalength,1);
if any(strcmp(ctns,'match'));
    for iiqc = 1:qcn.match
        thisCrit = qcnames.match{iiqc};
        thisCritData = DataTable.(thisCrit);
        switch types.match.(thisCrit)
            case 'char'
                chartempmatches = zeros(initialdatalength,1);
                for iicrit=1:numel(qryStruct.match.(thisCrit))
                    chartempmatches = chartempmatches | ...
                        strcmp(qryStruct.match.(thisCrit){iicrit}, thisCritData);
                end
                matches = matches & chartempmatches;
            case {'double', 'datetime'}
                matches = matches & ismember( thisCritData, ...
                    qryStruct.match.(thisCrit){:} );
            otherwise
                error('unknown type')
        end
    end
    if ~any(matches)
        warning(sprintf(['\n**** No Matches Were Found ****\n',...
                           '**** Output Table Is Empty ****']));
    end
end

%% Part II : Identify data section of (Range) Criteria
inRanges = ones(initialdatalength,1);
if any(strcmp(ctns,'range'));
    for iiqc = 1:qcn.range
        thisCrit = qcnames.range{iiqc};
        thisCritData = DataTable.(thisCrit);
        if isempty(qryStruct.range.(thisCrit){1})
            continue
        end
        switch types.range.(thisCrit)
            case {'double','datetime'}
                inRanges = inRanges & ...
                  ( thisCritData >= qryStruct.range.(thisCrit){1}(1) & ...
                    thisCritData <= qryStruct.range.(thisCrit){1}(2) );
            otherwise
                error('unknown type')
        end
    end
    if ~any(inRanges)
        warning('*** No Table Elements Within Queried Range Were Found ***')
    end
end

%% Part III : COMBINE MATCHES AND INRANGES
okElmnts = matches & inRanges;

%% Part IV : Extract data section of interest
%       (delete variables not queried from data structure and use
%        criteria for queried variables that are being retained)
% Unused Code:: % varsToRemove = [];
newdatatable=table();
for iiav = 1:numOfTableVars
    if ~any(strcmp(varnameq,allVarNames(iiav)))
        % Unused Code:: % varsToRemove = [varsToRemove,iiav];
    else
        switch class(DataTable.(allVarNames{iiav})(1))
            case {'double' , 'single'}
                % Unused Code:: % scratch = num2cell( DataTable.(allVarNames{iiav})(okElmnts) );
                % Unused Code:: % newdatatable.(allVarNames{iiav}) = cell2mat(scratch);
                newdatatable.(allVarNames{iiav}) = DataTable.(allVarNames{iiav})(okElmnts);
            case 'cell'
                newdatatable.(allVarNames{iiav}) = DataTable.(allVarNames{iiav})(okElmnts);
            case 'char'
                newdatatable.(allVarNames{iiav}) = DataTable.(allVarNames{iiav})(okElmnts,:);
            case 'datetime'
                newdatatable.(allVarNames{iiav}) = DataTable.(allVarNames{iiav})(okElmnts);
            otherwise
                error(['unknown type: ',class(DataTable.(allVarNames{iiav})(1))])
        end
        
        % Transfer Old DataTable Units and Descriptions to New DataTable
        varqTableIdx=strcmp(DataTable.Properties.VariableNames,(allVarNames{iiav}));
        newVidx = strcmp(newdatatable.Properties.VariableNames,(allVarNames{iiav}));
        if numel(DataTable.Properties.VariableNames) == numel(DataTable.Properties.VariableUnits)
            newdatatable.Properties.VariableUnits{newVidx} = DataTable.Properties.VariableUnits{varqTableIdx};
        end
        if numel(DataTable.Properties.VariableNames) == numel(DataTable.Properties.VariableDescriptions)
            newdatatable.Properties.VariableDescriptions{newVidx} = DataTable.Properties.VariableDescriptions{varqTableIdx};
        end
    end
end

end
%% END OF MAIN FUNCTION
% Local subfunctions follow...

function ok = chkQCriteria(x)
if ~isstruct(x)
    error('qCriteriaStruct must be a structure')
end

validFieldNames = {'match','range'};
fields = fieldnames(x);
if numel(fields) > 2
    error('First level fields (Query Type Names) in qryStruct can only be ''match'' and/or ''range''')
end
for i=1:numel(fields)
  if all(~strcmpi(fields{i}, validFieldNames))
      error('First level fields (Query Type Names) in qryStruct can only be ''match'' and/or ''range''')
  end
end
ok = 1;
end

function qryStruct = parseMorRcriteriaCell(x,qryStruct,qryType)
    nCrit = length(x); % count criteria
    if round(nCrit/2)~=nCrit/2
        error('%s input needs criteriaName/criteriaValue pairs',qryType)
    end

    for pair = reshape(x,2,[]) % pair is {criteriaName;criteriaValue}
        qryStruct.(qryType).(pair{1}) = pair(2);
    end
end

%#ok<*SPWRN>
%#ok<*SPERR>
