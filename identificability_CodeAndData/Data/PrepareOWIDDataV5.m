% Filter archived OWID ECDC US data
clear; clc;

projectDir = pwd;

while ~exist(fullfile(projectDir,'Data'),'dir')
    parentDir = fileparts(projectDir);

    if strcmp(parentDir,projectDir)
        error('Could not locate project root.');
    end

    projectDir = parentDir;
end

dataDir = fullfile(projectDir,'Data');

inputFile = fullfile(dataDir,'full_raw_data.csv');

if ~isfile(inputFile)
    inputFile = fullfile(projectDir,'full_raw_data.csv');
end

raw = readtable(inputFile,'TextType','string');

if isdatetime(raw.date)
    rawDates = raw.date;
else
    rawDates = datetime(raw.date,'InputFormat','yyyy-MM-dd');
end

locations = string(raw.location);
newCases = double(raw.new_cases);


startFull = datetime(2019,12,31);
endFull   = datetime(2020,6,30);

keep = locations == "United States" & ...
       rawDates >= startFull & rawDates <= endFull;

usDates = rawDates(keep);
usCases = newCases(keep);

[usDates,order] = sort(usDates);
usCases = usCases(order);

dateText = string(cellstr(datestr(usDates,'yyyy-mm-dd')));

fullUS = table(dateText(:),usCases(:), ...
    'VariableNames',{'date','new_cases'});


startWindow = datetime(2020,1,21);
endWindow   = datetime(2020,3,30);

windowKeep = usDates >= startWindow & usDates <= endWindow;

windowDates = usDates(windowKeep);
windowCases = usCases(windowKeep);

windowDateText = string(cellstr(datestr(windowDates,'yyyy-mm-dd')));

window70 = table((1:numel(windowDates))',windowDateText(:),windowCases(:), ...
    'VariableNames',{'day','date','new_cases'});
%% save

fullOutput = fullfile(dataDir, ...
    'owid_archived_ecdc_us_daily_cases.csv');

windowOutput = fullfile(dataDir, ...
    'real_data_window.csv');

writetable(fullUS,fullOutput);
writetable(window70,windowOutput);

fprintf('\nFiltered datasets saved:\n');
fprintf('%s\n',fullOutput);
fprintf('%s\n',windowOutput);