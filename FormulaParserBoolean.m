classdef FormulaParserBoolean < handle
    % This class definition is based on 'FormulaParsed' of Thomas Pfau
    % (2016), which parses logic formulas.
    % It has been modified in order to fulfill the De Morgan's laws for 
    % the inhibitory rules presented in Boolean rules. 
    %
    % Morgan's Laws say that in Boolean algebra 
    % not (A or B) = (not A) and (not B) /// 
    % not (A and B) = (not A) or (not B)
    
    % Two functions: createNodeStructure and createNodeStructure_not. This
    % last one is used when a ! is present in a node and during all its
    % layers until other ! is present. The function used at that time is
    % createNodeStructure, because the negation of the negation is the
    % original one. 
    
    properties
        literalpat = 'x\([0-9]+\)';
        pat = '(?!\(.*\(.*?\))(\((.*?)\))'
        subformulas;
    end
    
    methods
        function obj = FormulaParserBoolean()
        % Default FormulaParserBoolean constructor.        
        % USAGE:
        %    obj = FormulaParserBoolean()                
        %    
        % OUTPUTS:
        %    obj:    The FormulaParser Object
        %
        obj.subformulas = containers.Map();
        end
        
        function Head = parseFormulaBoolean_on(self,formula)            
        % Parse a Formula in the COBRA rules format.
        % USAGE:
        %    Head = FormulaParserBoolean.parseFormulaBoolean(formula)                
        %    
        % INPUTS:
        %    formula:   A String of a Boolean rules formula in rules format 
        %                ( &/| as operators, x(1) as literal symbols
        %
        % OUTPUTS:
        %    Head:       The Head of a Tree representing the formula
        %  
            id = 1;
            % Lets replace our rules (x([0-9])) by the corresponding number
            formula = regexprep(formula,'x\(([0-9]+)\)','$1');
            % For workability with grRules, we will also replace any version
            % of " and " by & and any version of " or " by |
            formula = regexprep(formula,'([ \)])and([ \(])','$1&$2','ignorecase');
            formula = regexprep(formula,'([ \)])or([ \(])','$1|$2','ignorecase');
            exp = regexp(formula,self.pat,'match');
            newf = formula;
            if not(isempty(exp))
                self.subformulas(['$' num2str(id)]) = exp;
                newf = strrep(formula,exp{1},['$' num2str(id)]);
                finalid = ['$' num2str(id)];
                id = id +1;
            else
                if not(length(formula) == 0)
                    exp = formula;
                    self.subformulas(['$' num2str(id)]) = {exp};
                    newf = strrep(formula,exp,['$' num2str(id)]);
                    finalid = ['$' num2str(id)];
                    id = id+1;
                end
            end
            while not(strcmp(newf,formula))
                exp = regexp(newf,self.pat,'match');
                formula = newf;
                if not(isempty(exp))
                    self.subformulas(['$' num2str(id)]) = exp;
                    newf = strrep(newf,exp{1},['$' num2str(id)]);
                    finalid = ['$' num2str(id)];
                    id = id +1;
                else
                    if length(formula) > length(['$' num2str(id)])
                        exp = formula;
                        self.subformulas(['$' num2str(id)]) = {exp};
                        newf = strrep(formula,exp,['$' num2str(id)]);
                        finalid = ['$' num2str(id)];
                        id = id+1;
                    end
                end
                
            end
            Head = self.createNodeStructure(finalid);
        end
        function Head = parseFormulaBoolean_off(self,formula)            
        % Parse a Formula in the COBRA rules format (as detailed above).
        % USAGE:
        %    Head = FormulaParser.parseFormula(formula)                
        %    
        % INPUTS:
        %    formula:   A String of a GPR formula in rules format ( &/| as
        %               operators, x(1) as literal symbols
        %
        % OUTPUTS:
        %    Head:       The Head of a Tree representing the formula
        %  
            id = 1;
            %Lets replace our rules (x([0-9])) by the corresponding number
            formula = regexprep(formula,'x\(([0-9]+)\)','$1');
            %For workability with grRules, we will also replace any version
            %of " and " by & and any version of " or " by |
            formula = regexprep(formula,'([ \)])and([ \(])','$1&$2','ignorecase');
            formula = regexprep(formula,'([ \)])or([ \(])','$1|$2','ignorecase');
            exp = regexp(formula,self.pat,'match');
            newf = formula;
            if not(isempty(exp))
                self.subformulas(['$' num2str(id)]) = exp;
                newf = strrep(formula,exp{1},['$' num2str(id)]);
                finalid = ['$' num2str(id)];
                id = id +1;
            else
                if not(length(formula) == 0)
                    exp = formula;
                    self.subformulas(['$' num2str(id)]) = {exp};
                    newf = strrep(formula,exp,['$' num2str(id)]);
                    finalid = ['$' num2str(id)];
                    id = id+1;
                end
            end
            while not(strcmp(newf,formula))
                exp = regexp(newf,self.pat,'match');
                formula = newf;
                if not(isempty(exp))
                    self.subformulas(['$' num2str(id)]) = exp;
                    newf = strrep(newf,exp{1},['$' num2str(id)]);
                    finalid = ['$' num2str(id)];
                    id = id +1;
                else
                    
                    if length(formula) > length(['$' num2str(id)])
                        exp = formula;
                        self.subformulas(['$' num2str(id)]) = {exp};
                        newf = strrep(formula,exp,['$' num2str(id)]);
                        finalid = ['$' num2str(id)];
                        id = id+1;
                    end
                end
                
            end
            Head = self.createNodeStructure_not(finalid);
            
        end
                        
        function HeadNode = createNodeStructure(self,finalid)
            if self.subformulas.isKey(finalid)
                currentstring = self.subformulas(finalid);
                currentstring = currentstring{1};
                pos = find(regexp(currentstring,'( or )|( OR )|( Or )| ?\| ?| ?\|\|'));
                if not(isempty(pos))
                    HeadNode = OrNode();
                    literals = strsplit(regexprep(currentstring,'( or )|( OR )|( Or )|( ?\| ?)|( ?\|\| ?)', '$$'),'$$');
                    for i=1:numel(literals)
                        pos = find(regexp(currentstring,'( and )|( AND )|( And )|( ?\& ?)|( ?\&\& ?)'));
                        if not(isempty(pos))
                            NewNode = AndNode();
                            HeadNode.addChild(NewNode);
                            andliterals = strsplit(regexprep(literals{i},'( and )|( AND )|( And )|( ?\& ?)|( ?\&\& ?)', '$$'),'$$');
                            for j = 1:numel(andliterals)
                                literal = andliterals{j};
                                literal = regexprep(literal,'\[|\]|\{|\}|\(|\)|\s','');
                                if find(regexp(literal,'(!)|(not)'))
                                    literal = regexprep(literal, '(!)|(not)', '');
                                    NewNode.addChild(self.createNodeStructure_not(literal));
                                else
                                    NewNode.addChild(self.createNodeStructure(literal));
                                end
                            end
                        else
                            literal = literals{i};
                            literal = regexprep(literal,'\[|\]|\{|\}|\(|\)|\s','');
                            if find(regexp(literal,'(!)|(not)'))
                                literal = regexprep(literal, '(!)|(not)', '');
                                HeadNode.addChild(self.createNodeStructure_not(literal));
                            else
                                HeadNode.addChild(self.createNodeStructure(literal));
                            end
                            
                        end
                        
                    end
                else
                    pos = find(regexp(currentstring,'( and )|( AND )|( And )|( ?\& ?)|( ?\&\& ?)'));
                    if not(isempty(pos))
                        % there are only Ands in this node
                        HeadNode = AndNode();
                        literals = strsplit(regexprep(currentstring,'( and )|( AND )|( And )|( ?\& ?)|( ?\&\& ?)', '$$'),'$$');
                        for i=1:numel(literals)
                            literal = literals{i};
                            literal = regexprep(literal,'\[|\]|\{|\}|\(|\)|\s','');
                            if find(regexp(literal,'(!)|(not)'))
                                literal = regexprep(literal, '(!)|(not)', '');
                                HeadNode.addChild(self.createNodeStructure_not(literal));
                            else
                                HeadNode.addChild(self.createNodeStructure(literal));
                            end
                        end
                    else
                        % We have literals which were seperated by brackets
                        % (or additional brackets). so just parse the
                        % current node ignoring the brackets.
                        nodestring = regexprep(currentstring,'\[|\]|\{|\}|\(|\)|\s','');
                        if find(regexp(nodestring,'(!)|(not)'))
                        	nodestring = regexprep(nodestring, '(!)|(not)', '');
                        	HeadNode = self.createNodeStructure_not(nodestring);
                        else
                        	HeadNode = self.createNodeStructure(nodestring);
                        end
                    end
                end
                
            else
                % this either represents a literal OR the Head node without
                % any brackets
                currentstring = finalid;
                pos = find(regexp(currentstring,'( or )|( OR )|( Or )|( \| )|( \|\| )'));
                if not(isempty(pos))
                    HeadNode = OrNode();
                    literals = strsplit(regexprep(currentstring,'( or )|( OR )|( Or )|( \| )|( \|\| )', '$$'),'$$');
                    for i=1:numel(literals)
                        pos = find(regexp(currentstring,'( and )|( AND )|( And )|( \& )|( \&\& )'));
                        if not(isempty(pos))
                            NewNode = AndNode();
                            HeadNode.addChild(NewNode);
                            andliterals = strsplit(regexprep(literals{i},'( and )|( AND )|( And )|( \& )|( \&\& )', '$$'),'$$');
                            for j = 1:numel(andliterals)
                                literal = andliterals{j};
                                literal = regexprep(literal,'\[|\]|\{|\}|\(|\)|\s','');
                                if find(regexp(literal,'(!)|(not)'))
                                    literal = regexprep(literal, '(!)|(not)', '');
                                    NewNode.addChild(self.createNodeStructure_not(literal));
                                else
                                    NewNode.addChild(self.createNodeStructure(literal));
                                end
                                
                            end
                        else
                            literal = literals{i};
                            literal = regexprep(literal,'\[|\]|\{|\}|\(|\)|\s','');
                            if find(regexp(literal,'(!)|(not)'))
                                literal = regexprep(literal, '(!)|(not)', '');
                                HeadNode.addChild(self.createNodeStructure_not(literal));
                            else
                                HeadNode.addChild(self.createNodeStructure(literal));
                            end
                            
                        end
                        
                    end
                else
                    pos = find(regexp(currentstring,'( and )|( AND )|( And )|( \& )|( \&\& )'));
                    if not(isempty(pos))
                        %there are only Ands in this node
                        HeadNode = AndNode();
                        literals = strsplit(regexprep(currentstring,'( and )|( AND )|( And )|( \& )|( \&\& )', '$$'),'$$');
                        for i=1:numel(literals)
                            literal = andliterals{j};
                            literal = regexprep(literal,'\[|\]|\{|\}|\(|\)|\s','');
                            if find(regexp(literal,'(!)|(not)'))
                                literal = regexprep(literal, '(!)|(not)', '');
                                HeadNode.addChild(self.createNodeStructure_not(literal));
                            else
                                HeadNode.addChild(self.createNodeStructure(literal));
                            end
                        end
                    else
                        %Now there are no ORs nor ANDs in the node, so it
                        %is a literal!
                        %This SHOULD be a node representing a single literal.
                        %thus we will remove any whitespace, and any brackets and
                        %create a literal node.
                        literalstring = regexprep(finalid,'\[|\]|\{|\}|\(|\)|\s','');
                        %We will create an Or Node with a single literal
                        %child.
                        HeadNode = OrNode();
                        if find(regexp(literalstring,'(!)|(not)'))
                            literalstring = regexprep(literalstring, '(!)|(not)', '');
                            HeadNode.addChild(LiteralNode([literalstring '_off']));
                        else
                            HeadNode.addChild(LiteralNode([literalstring '_on']));
                        end
                                                
                    end
                    
                end
                
            end            
            
        end
        
        function HeadNode_not = createNodeStructure_not(self,finalid)
        % When there is an inhibitory operator (!)|(not), Morgan's laws are
        % fulfilled.
            if self.subformulas.isKey(finalid)
                currentstring = self.subformulas(finalid);
                currentstring = currentstring{1};
                pos = find(regexp(currentstring,'( or )|( OR )|( Or )| ?\| ?| ?\|\|'));
                if not(isempty(pos))
                    HeadNode_not = AndNode();
                    literals = strsplit(regexprep(currentstring,'( or )|( OR )|( Or )|( ?\| ?)|( ?\|\| ?)', '$$'),'$$');
                    for i=1:numel(literals)
                        pos = find(regexp(currentstring,'( and )|( AND )|( And )|( ?\& ?)|( ?\&\& ?)'));
                        if not(isempty(pos))
                            NewNode = OrNode();
                            HeadNode_not.addChild(NewNode);
                            andliterals = strsplit(regexprep(literals{i},'( and )|( AND )|( And )|( ?\& ?)|( ?\&\& ?)', '$$'),'$$');
                            for j = 1:numel(andliterals)
                                literal = andliterals{j};
                                literal = regexprep(literal,'\[|\]|\{|\}|\(|\)|\s','');
                                if find(regexp(literal,'(!)|(not)'))
                                    literal = regexprep(literal, '(!)|(not)', '');
                                    NewNode.addChild(self.createNodeStructure(literal));
                                else
                                    NewNode.addChild(self.createNodeStructure_not(literal));
                                end                                
                            end
                        else
                            literal = literals{i};
                            literal = regexprep(literal,'\[|\]|\{|\}|\(|\)|\s','');
                            if find(regexp(literal,'(!)|(not)'))
                                literal = regexprep(literal, '(!)|(not)', '');
                                HeadNode_not.addChild(self.createNodeStructure(literal));
                            else
                                HeadNode_not.addChild(self.createNodeStructure_not(literal));
                            end
                        end
                        
                    end
                else
                    pos = find(regexp(currentstring,'( and )|( AND )|( And )|( ?\& ?)|( ?\&\& ?)'));
                    if not(isempty(pos))
                        %there are only Ands in this node, that will become
                        %Ors
                        HeadNode_not = OrNode();
                        literals = strsplit(regexprep(currentstring,'( and )|( AND )|( And )|( ?\& ?)|( ?\&\& ?)', '$$'),'$$');
                        for i=1:numel(literals)
                            literal = literals{i};
                            literal = regexprep(literal,'\[|\]|\{|\}|\(|\)|\s','');
                            if find(regexp(literal,'(!)|(not)'))
                                literal = regexprep(literal, '(!)|(not)', '');
                                HeadNode_not.addChild(self.createNodeStructure(literal));
                            else
                                HeadNode_not.addChild(self.createNodeStructure_not(literal));
                            end
                        end
                    else
                        %We have literals which were seperated by brackets
                        %(or additional brackets). so just parse the
                        %current node ignoring the brackets.
                        nodestring = regexprep(currentstring,'\[|\]|\{|\}|\(|\)|\s','');
                        if find(regexp(nodestring,'(!)|(not)'))
                        	nodestring = regexprep(nodestring, '(!)|(not)', '');
                        	HeadNode_not = self.createNodeStructure(nodestring);
                        else
                        	HeadNode_not = self.createNodeStructure_not(nodestring);
                        end
                    end
                end
                
            else
                %this either represents a literal OR the Head node without
                %any brackets
                currentstring = finalid;
                pos = find(regexp(currentstring,'( or )|( OR )|( Or )|( \| )|( \|\| )'));
                if not(isempty(pos))
                    HeadNode_not = AndNode();
                    literals = strsplit(regexprep(currentstring,'( or )|( OR )|( Or )|( \| )|( \|\| )', '$$'),'$$');
                    for i=1:numel(literals)
                        pos = find(regexp(currentstring,'( and )|( AND )|( And )|( \& )|( \&\& )'));
                        if not(isempty(pos))
                            NewNode = OrNode();
                            HeadNode_not.addChild(NewNode);
                            andliterals = strsplit(regexprep(literals{i},'( and )|( AND )|( And )|( \& )|( \&\& )', '$$'),'$$');
                            for j = 1:numel(andliterals)
                                literal = andliterals{j};
                                literal = regexprep(literal,'\[|\]|\{|\}|\(|\)|\s','');
                                if find(regexp(literal,'(!)|(not)'))
                                    literal = regexprep(literal, '(!)|(not)', '');
                                    NewNode.addChild(self.createNodeStructure(literal));
                                else
                                    NewNode.addChild(self.createNodeStructure_not(literal));
                                end
                            end
                        else
                            literal = literals{i};
                            literal = regexprep(literal,'\[|\]|\{|\}|\(|\)|\s','');
                            if find(regexp(literal,'(!)|(not)'))
                                literal = regexprep(literal, '(!)|(not)', '');
                                HeadNode_not.addChild(self.createNodeStructure(literal));
                            else
                                HeadNode_not.addChild(self.createNodeStructure_not(literal));
                            end
                        end
                        
                    end
                else
                    pos = find(regexp(currentstring,'( and )|( AND )|( And )|( \& )|( \&\& )'));
                    if not(isempty(pos))
                        %there are only Ands in this node, that will become
                        %Ors.
                        HeadNode_not = OrNode();
                        literals = strsplit(regexprep(currentstring,'( and )|( AND )|( And )|( \& )|( \&\& )', '$$'),'$$');
                        for i=1:numel(literals)
                            literal = andliterals{j};
                            literal = regexprep(literal,'\[|\]|\{|\}|\(|\)|\s','');
                            %                                literal = literal{1};
                            if find(regexp(literal,'(!)|(not)'))
                                literal = regexprep(literal, '(!)|(not)', '');
                                HeadNode_not.addChild(self.createNodeStructure(literal));
                            else
                                HeadNode_not.addChild(self.createNodeStructure_not(literal));
                            end
                        end
                    else
                        %Now there are no ORs nor ANDs in the node, so it
                        %is a literal!
                        %This SHOULD be a node representing a single literal.
                        %thus we will remove any whitespace, and any brackets and
                        %create a literal node.
                        literalstring = regexprep(finalid,'\[|\]|\{|\}|\(|\)|\s','');
                        %We will create an Or Node with a single literal
                        %child.
                        HeadNode_not = OrNode();
                        if find(regexp(literalstring,'(!)|(not)'))
                            literalstring = regexprep(literalstring, '(!)|(not)', '');
                            HeadNode_not.addChild(LiteralNode([literalstring '_on']));
                        else
                            HeadNode_not.addChild(LiteralNode([literalstring '_off']));
                        end                      
                    end
                    
                end
                
            end            
            
        end
        
    end
    
end

