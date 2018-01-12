%% common_parameters
% 
% Sets up parameters common to the different scripts.
%
%% Information
% 
% * *Author*: <mailto:silvello@dei.unipd.it Gianmaria Silvello>
% * *Version*: 1.00
% * *Since*: 1.00
% * *Requirements*: MATTERS 1.0 or higher; Matlab 2015b or higher
% * *Copyright:* (C) 2016 <http://ims.dei.unipd.it/ Information 
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

diary off;
warning off;

%% Path Configuration


EXPERIMENT.path.base = '/Users/silvello/Dropbox/Documents/Projects/active/PREFORMA/matlab/';


% The path for the measures, i.e. the GoP
EXPERIMENT.path.measure = sprintf('%1$s%2$s%3$s', EXPERIMENT.path.base, 'measure', filesep);

% The path for data of runs
EXPERIMENT.path.runDir = sprintf('%1$s%2$s%3$s%4$s%5$s%6$s%7$s', EXPERIMENT.path.base, 'data', filesep, 'run_test', filesep, 'text', filesep);

% The pool file name
EXPERIMENT.poolFileName = 'preforma-text-test-groundtruth.txt';

% The path for data of pools
EXPERIMENT.path.pool = sprintf('%1$s%2$s%3$s%4$s%5$s%6$s%7$s%8$s', EXPERIMENT.path.base, 'data', filesep, 'pool', filesep, 'text', filesep,  EXPERIMENT.poolFileName);


%% General Configuration

EXPERIMENT.run.id = 'PREFORMARuns';

EXPERIMENT.pool.id = 'PREFORMAPool';

EXPERIMENT.pool.RelevanceGrades = 0:1;
EXPERIMENT.pool.RelevanceDegrees = {'NotRelevant', 'Relevant'};
EXPERIMENT.pool.delimiter = 'tab';
