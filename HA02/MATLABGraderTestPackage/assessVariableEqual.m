function assessVariableEqual(VariableName, ExpectedValue, varargin)
% assessVariableEqual checks whether variable properties are equal to expected values
%
%     assessVariableEqual(variableName,expectedValue) determines whether the
%     value of variableName is equal to expectedValue.
%
%     assessVariableEqual(variableName,expectedValue,Name,Value) uses
%     additional options specified by one or more Name,Value pair arguments.
%
%     - Feedback - Specify additional feedback to display to the learner.
%
%     - RelativeTolerance - Specify the relative tolerance to apply to the submitted variable value. The default relative tolerance is 0.1%.
%
%     - Absolute Tolerance - Specify the absolute tolerance to apply to the submitted variable value. The default absolute tolerance is 1e-4.
%
%     If no tolerance is specified, the function applies both tolerances with their default values. If either tolerance passes, the variable is considered equal.
%
%     Examples
%     --------
%
%     % Assess value of variable studentValue. The expected
%     % value is 2. The submitted value is 1.
%     assessVariableEqual('studentValue', '2')
%
%        Variable studentValue has an incorrect value: 1
%
%
%     % Assess value of variable studentArray and provide
%     % additional feedback. The expected value is 'two'. The
%     % submitted value is 'four'.
%     assessVariableEqual('studentArray', 'two', 'Feedback', 'Refer to the week2 handout on Prime Numbers')
%
%        Variable studentString has an incorrect value: four.
%        Refer to the week2 handout on Prime Numbers.
%
%
%     % Assess value of variable studentArray using a relative
%     % tolerance of 0.5%. The variable has an expected value of [1 2 3 4 6 4 3 4 5]. The submitted value is [1 2 3 4 6 4 3 4 1].
%     assessVariableEqual('studentArray',[1 2 3 4 6 4 3 4 5],'RelativeTolerance',0.5)
%
%        Variable studentValue has an incorrect value: 1 2 3 4 6 4 3 4 1
%
%
%     See also assessFunctionPresence, assessFunctionAbsence, matlab.unittest.constraints.RelativeTolerance, matlab.unittest.constraints.AbsoluteTolerance

%% Default Values
    DefaultRelativeTolerance = 0.001;
    DefaultAbsoluteTolerance = 0.0001;
    DefaultUserFeedback      = '';

    %% Input Handling
    inputHandler = inputParser;
    inputHandler.PartialMatching = false;
    addRequired(inputHandler,'VariableName',@ischar);
    addRequired(inputHandler,'ExpectedValue');
    addParameter(inputHandler,'RelativeTolerance', DefaultRelativeTolerance, @isnumeric);
    addParameter(inputHandler,'AbsoluteTolerance', DefaultAbsoluteTolerance, @isnumeric);
    addParameter(inputHandler,'Feedback', DefaultUserFeedback, @ischar);
    parse(inputHandler, VariableName, ExpectedValue, varargin{:});

    %% Variable Existence Check
    cmd = sprintf('exist(''%s'', ''var'')', VariableName);
    variableExists = isequal(evalin('caller', cmd), 1);
    if ~variableExists
         variableNotPresentMsg = sprintf('%s%s%s\n%s', 'The submission must contain a variable named ', ...
                     VariableName, ...
                     '.', inputHandler.Results.Feedback);
        ME = MException('AssessmentToolbox:Feedback:VariableNotPresent', variableNotPresentMsg);
        throwAsCaller(ME);
    end

    %% Get Value from VariableName
    ActualValue = evalin('caller', VariableName);

    %% Calling IsValueEqualTo Constraint

    variableEqualConstraint = IsValueEqualTo(ExpectedValue, varargin{:}, 'VariableName', VariableName);

    if ~variableEqualConstraint.satisfiedBy(ActualValue)
        diagnostics = variableEqualConstraint.getDiagnosticFor(ActualValue);
        diagnostics.diagnose;
        ME = MException(diagnostics.DiagnosticIdentifier, diagnostics.DiagnosticResult);
        throwAsCaller(ME)
    end
end