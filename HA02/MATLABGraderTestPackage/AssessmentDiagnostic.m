classdef AssessmentDiagnostic < matlab.unittest.diagnostics.Diagnostic
    properties(SetAccess=private, GetAccess=public)
        DiagnosticIdentifier
    end
    
    methods
        function diag = AssessmentDiagnostic(identifier, value)
            % AssessmentDiagnostic - Class constructor
            %
            %   AssessmentDiagnostic(IDENTIFIER, VALUE) creates a new AssessmentDiagnostic instance using
            %   the VALUE provided.
            %
            %   Examples:
            %
            %       AssessmentDiagnostic('identifier', 'Diagnostic text');
            %       AssessmentDiagnostic('identifier', sprintf('This text is first created using %s', 'sprintf'));
            %
            validateattributes(identifier, {'char'}, {'2d'}, '', 'identifier');
            validateattributes(value, {'char'}, {'2d'}, '', 'value');
            diag.DiagnosticResult = strjoin(cellstr(value),'\n');
            diag.DiagnosticIdentifier = identifier;
        end
        function diagnose(~)
        end
    end
    
end