classdef IsValueEqualTo < matlab.unittest.constraints.Constraint

    properties(SetAccess=immutable)
        Expected
    end

    properties(Constant)
        DefaultRelativeTolerance = 0.001;
        DefaultAbsoluteTolerance = 0.0001;
        DefaultUserFeedback = '';
        DefaultVariableName = '';
    end

    properties (Access=private)
        ActualVariableName;
        Feedback;
        RelativeToleranceValue;
        AbsoluteToleranceValue;
        AbsoluteToleranceSpecified = false;
        RelativeToleranceSpecified = false;
    end

    methods
        function constraint = IsValueEqualTo(ExpectedValue, varargin)
            inputHandler = inputParser;
            inputHandler.PartialMatching = false;
            addRequired(inputHandler,'ExpectedValue');
            addParameter(inputHandler,'RelativeTolerance', constraint.DefaultRelativeTolerance, @isnumeric);
            addParameter(inputHandler,'AbsoluteTolerance', constraint.DefaultAbsoluteTolerance, @isnumeric);
            addParameter(inputHandler,'Feedback', constraint.DefaultUserFeedback, @ischar);
            addParameter(inputHandler,'VariableName', constraint.DefaultVariableName, @ischar);

            parse(inputHandler, ExpectedValue, varargin{:})

            if ~isempty(inputHandler.Results.VariableName)
                constraint.ActualVariableName = [ ' ' inputHandler.Results.VariableName];
            end

            constraint.AbsoluteToleranceSpecified = ~any(strcmp('AbsoluteTolerance', inputHandler.UsingDefaults));
            constraint.RelativeToleranceSpecified = ~any(strcmp('RelativeTolerance', inputHandler.UsingDefaults));

            constraint.Expected = inputHandler.Results.ExpectedValue;
            constraint.RelativeToleranceValue = inputHandler.Results.RelativeTolerance;
            constraint.AbsoluteToleranceValue = inputHandler.Results.AbsoluteTolerance;
            constraint.Feedback = inputHandler.Results.Feedback;
        end

        function bool = satisfiedBy(constraint, ActualValue) 
            bool = constraint.hasSameSize(ActualValue) && ...
                   constraint.hasSameDataType(ActualValue) && ...
                   constraint.hasSameValue(ActualValue);
        end

        function diag = getDiagnosticFor(constraint, ActualValue)
            import matlab.unittest.internal.diagnostics.ConstraintDiagnosticFactory;
            import matlab.unittest.internal.diagnostics.DiagnosticSense;
            constraintDiagnostic = matlab.unittest.diagnostics.ConstraintDiagnostic;

            if ~constraint.hasSameDataType(ActualValue)
                datatypeMismatchMsg = sprintf('%s%s%s%s%s%s%s', ...
                    'Variable', constraint.ActualVariableName, ' must be of data type ', class(constraint.Expected), '.', ' It is currently of type ', class(ActualValue), '. Check where the variable is assigned a value.');

                diag = constraint.generateDiagnostics('AssessmentToolbox:Feedback:DataTypeMismatch', datatypeMismatchMsg);
            elseif ~constraint.hasSameSize(ActualValue)
                sizeMismatchMsg = sprintf('%s%s%s%s%s%s%s', ...
                    'Variable', constraint.ActualVariableName, ' must be of size ', mat2str(size(constraint.Expected)), '.', ' It is currently of size ', mat2str(size(ActualValue)), '. Check where the variable is assigned a value.');

                diag = constraint.generateDiagnostics('AssessmentToolbox:Feedback:SizeMismatch', sizeMismatchMsg);
            elseif ~constraint.hasSameValue(ActualValue)

                displayPattern = '%s%s%s';
                valueMismatchMsg = sprintf(displayPattern, ...
                    'Variable', constraint.ActualVariableName, ' has an incorrect value.');

                diag = constraint.generateDiagnostics('AssessmentToolbox:Feedback:ValueMismatch', valueMismatchMsg);
            else
                diag = ConstraintDiagnosticFactory.generatePassingDiagnostic(constraint, ...
                    DiagnosticSense.Positive);
            end
        end
    end

    methods(Access=private)
        function bool = hasSameDataType(constraint, actual)
            bool = isequal(class(actual), class(constraint.Expected));
        end


        function bool = hasSameSize(constraint, actual)
            if  ~ischar(actual)
                bool = isequal(size(actual), size(constraint.Expected));
            else
                bool = true;
            end
        end

        function bool = hasSameValue(constraint, actual)
            import matlab.unittest.constraints.IsEqualTo;
            import matlab.unittest.constraints.AbsoluteTolerance;
            import matlab.unittest.constraints.RelativeTolerance;
            
            isEqualToConstraint = IsEqualTo(constraint.Expected, 'Using', [SymbolicExpressionComparator, IsEqualTo.DefaultComparator]);
            
            if (constraint.AbsoluteToleranceSpecified && constraint.RelativeToleranceSpecified) || (~constraint.AbsoluteToleranceSpecified && ~constraint.RelativeToleranceSpecified)
                absoluteTolerance = AbsoluteTolerance(single(constraint.AbsoluteToleranceValue), double(constraint.AbsoluteToleranceValue));
                relativeTolerance = RelativeTolerance(single(constraint.RelativeToleranceValue), double(constraint.RelativeToleranceValue));
                isEqualToConstraint = isEqualToConstraint.within(absoluteTolerance | relativeTolerance);
            elseif constraint.RelativeToleranceSpecified
                isEqualToConstraint = isEqualToConstraint.within(RelativeTolerance(constraint.RelativeToleranceValue));
            elseif constraint.AbsoluteToleranceSpecified
                isEqualToConstraint = isEqualToConstraint.within(AbsoluteTolerance(constraint.AbsoluteToleranceValue));
            end
            bool = isEqualToConstraint.satisfiedBy(actual);
        end

        function stringDiagnostic = generateDiagnostics(constraint, identifier, errorMessage)
            if ~isempty(constraint.Feedback)
                stringDiagnostic = AssessmentDiagnostic(identifier, sprintf('%s\n%s', errorMessage, constraint.Feedback));
            else
                stringDiagnostic = AssessmentDiagnostic(identifier, errorMessage);
            end
        end
    end
end