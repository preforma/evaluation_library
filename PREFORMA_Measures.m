%% PREFORMA_Measures
% 
% Computes the PREFORMA Measures matrixes.

%% Synopsis
%
%   [accuracy, AUC, LAM, consistency] = PREFORMA_Measures(pool, runSet)
%  
% Note that average precision will be NaN when there are no relevant
% documents for a given topic in the pool (this may happen due to the way
% in which relevance degrees are mapped to binary relevance).
%
% *Parameters*
%
% * *|pool|* - the pool to be used to assess the run(s). It is a table in the
% same format returned by <../io/importPoolFromFileTRECFormat.html 
% importPoolFromFileTRECFormat>;
% * *|runSet|* - the run(s) to be assessed. It is a table in the same format
% returned by <../io/importRunFromFileTRECFormat.html 
% importRunFromFileTRECFormat> or by <../io/importRunsFromDirectoryTRECFormat.html 
% importRunsFromDirectoryTRECFormat>;
%
% *Returns*
%
% * |accuracy|  - measures the overall effectiveness of a conformance checker.
% * *|AUC|* - area under the curve (AUC) measures the ability of a conformance checker to avoid false classification.
% * *|LAM|* - logistic average misclassification rate (LAM) is the geometric mean of the odds of compliance and 
%               not-compliance misclassification, converted back to a proportion .
% * *|consistency|* - assesses the ability of a conformance checker to adhere to some constraint of separation of 
%                       C0 from the other classes.

%% Information
% 
% * *Author*: <mailto:silvello@dei.unipd.it Gianmaria Silvello>
% * *Version*: 1.00
% * *Since*: 1.00
% * *Requirements*: Matlab 2013b or higher
% * *Copyright:* (C) 2013-2017 <http://ims.dei.unipd.it/ Information 
% Management Systems> (IMS) research group, <http://www.dei.unipd.it/ 
% Department of Information Engineering> (DEI), <http://www.unipd.it/ 
% University of Padua>, Italy
% * *License:* <http://www.apache.org/licenses/LICENSE-2.0 Apache License, 
% Version 2.0>

%%
%{
 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at
 
      http://www.apache.org/licenses/LICENSE-2.0
 
 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
%}

%%
function [accuracy, AUC, LAM, consistency] = PREFORMA_Measures(pool, runSet)
    
    % check that we have the correct number of input arguments. 
    narginchk(2, 2);
    
    % assess the runs
    [assessedRunSet, poolStats, runSetStats, ~] = assess(pool, runSet, 'MapToBinaryRelevance', 'Lenient');
    
     % the topic currently under processing
    ct = 1;
        
    % the run currently under processing
    cr = 1;
    
    docs = pool{1,1}{1,1}{:,1};
    for k = 2 : size(pool, 1)
        docs = union(pool{k,1}{1,1}{:,1}, docs);
    end
    
        
    % compute the measure topic-by-topic
    [confusionMatrix] = rowfun(@processTopic, assessedRunSet, 'OutputVariableNames', runSet.Properties.VariableNames, 'OutputFormat', 'table', 'ExtractCellContents', true, 'SeparateInputs', false);
    confusionMatrix.Properties.UserData.identifier = assessedRunSet.Properties.UserData.identifier;
    confusionMatrix.Properties.UserData.pool = pool.Properties.UserData.identifier;
    
    confusionMatrix.Properties.UserData.name = 'confusion';
    confusionMatrix.Properties.UserData.shortName = 'CM';
    
    accuracy = cell(size(confusionMatrix, 1), size(confusionMatrix, 2));
    accuracy = cell2table(accuracy);
    
    accuracy.Properties.UserData.identifier = assessedRunSet.Properties.UserData.identifier;
    accuracy.Properties.UserData.pool = pool.Properties.UserData.identifier;
    accuracy.Properties.VariableNames = confusionMatrix.Properties.VariableNames;
    accuracy.Properties.RowNames = confusionMatrix.Properties.RowNames;
    
    accuracy.Properties.UserData.name = 'accuracy';
    accuracy.Properties.UserData.shortName = 'A';
    

    
    AUC = cell(size(confusionMatrix, 1), size(confusionMatrix, 2));
    AUC = cell2table(AUC);
    
    AUC.Properties.UserData.identifier = assessedRunSet.Properties.UserData.identifier;
    AUC.Properties.UserData.pool = pool.Properties.UserData.identifier;
    AUC.Properties.VariableNames = confusionMatrix.Properties.VariableNames;
    AUC.Properties.RowNames = confusionMatrix.Properties.RowNames;
    
    AUC.Properties.UserData.name = 'AUC';
    AUC.Properties.UserData.shortName = 'AUC';
    
    
    %% Accuracy, AUC and LAM
    
    LAM = cell(size(confusionMatrix, 1), size(confusionMatrix, 2));
    LAM = cell2table(LAM);
    
    LAM.Properties.UserData.identifier = assessedRunSet.Properties.UserData.identifier;
    LAM.Properties.UserData.pool = pool.Properties.UserData.identifier;
    LAM.Properties.VariableNames = confusionMatrix.Properties.VariableNames;
    LAM.Properties.RowNames = confusionMatrix.Properties.RowNames;
    
    LAM.Properties.UserData.name = 'LAM';
    LAM.Properties.UserData.shortName = 'LAM';
    
    
    % classes
    for i = 1 : size(confusionMatrix, 1)
        %runs
        for j = 1 : size(confusionMatrix, 2)
            tmp = confusionMatrix{i,j};
            TP = tmp{1, 'TP'};
            FP = tmp{1, 'FP'};
            TN = tmp{1, 'TN'};
            FN = tmp{1, 'FN'};
            
            accuracy{i,j} = {(abs(TP) + abs(TN)) / (abs(TP) + abs(TN) + abs(FP) + abs(FN))};
            AUC{i,j} = {((abs(TP) / (abs(TP) + abs(FN))) +  (abs(TN) / (abs(TN) + abs(FP))))/2};
            
            fpr = (abs(FP) / (abs(FP) + abs(TN)));
            fnr = (abs(FN) / (abs(FN) + abs(TP)));
            
            a = log(fpr/(1-fpr));
            b = log(fnr/(1-fnr));
            c = (a+b)/2;
            
            LAM{i,j} = {exp(c)/(1+exp(c))};
            
        end;
    end;
    
    
    %% consistency 
    consistency = zeros(size(runSet, 1), size(runSet, 2));
    consistency = array2table(consistency);
    consistency.Properties.VariableNames = runSet.Properties.VariableNames;
    consistency.Properties.RowNames = runSet.Properties.RowNames;
    
    consistency{1,1} = 1;
    
    for i = 1 : size(runSet, 2)
        iterDocs = runSet{1, i}{1, 1}{:, 1};
        for t1 = 2 : size(runSet, 1)
            currDocs = runSet{t1, i}{1, 1}{:, 1};
            numerator = 0;
            %denominator = 0;
            %for t2 = 1 : size(runSet, 1)
             %   if (t1 ~= t2)
                    if(~isempty(intersect(currDocs, iterDocs)))
                        numerator = size(intersect(currDocs, iterDocs), 1);
                    end;
                    denominator = size(currDocs, 1);
                %end;
            %end;
            consistency{t1,i} = 1 - (numerator / denominator);
        end;
    end;
    
       
    %%
    
    % compute the measure for a given topic over all the runs
    function [varargout] = processTopic(topic)
        
        % reset the index of the run under processing for each topic
        cr = 1;
         
        % compute the measure only on those column which contain the
        % actual runs
        [varargout] = cellfun(@processRun, topic);
         
        % increment the index of the current topic under processing
        ct = ct + 1;    

        
        %% 
        
        % compute the measure for a given topic of a given run
        function [measure] = processRun(runTopic)

            % the total number of positives
            recallBase = poolStats{ct, 'RecallBase'};
            
            % the run true positives: all the relevant doc retrieved
            TP = runSetStats{ct, cr}.Relevant;
            
            % the run false positives: all the not relevant doc retrieved
            FP = size(runSet{ct, cr}{:},1) - runSetStats{ct, cr}.Relevant;
            
            % the run false negatives: all the relevant doc (in the pool)
            % which are not retrieved by the run
            FN = recallBase - runSetStats{ct, cr}.relevantRetrieved;
            
            % the true negatives are the number of unique documents in the
            % collection setminus the docs in the run
            TN = size(setdiff(docs, runSet{ct, cr}{1,1}{:, 1}), 1);
            
            
            A = table(TP, FP, TN, FN);
            A.Properties.VariableNames = {'TP', 'FP', 'TN', 'FN'};
           
            measure = {A};
            
            % increment the index of the current run under processing
            cr = cr + 1;
        end
    end
    
end



