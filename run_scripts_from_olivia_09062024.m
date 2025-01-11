% Try running scripts that Olivias sent me on 09/06/24 meeting
% 1.) Create_Parameter_Set.m
% 2.) nParasiteGroups_ContRuptFunc.m
% pkg install -forge statistics ?
% pkg load statistics ?
work_dir='G:/My Drive/mariani_systems/malaria_ctmc_project/MS1ReviewCode/MS1ReviewCode'
addpath(work_dir);
par = Create_Parameter_Set;
%out = CTMC_ContRuptFunc(par);
nParasiteGroups_ContRuptFunc;

%Run R script from MATLAB/OCTAVE: --------------------------------------------

%https://www.mathworks.com/matlabcentral/answers/231000-how-to-run-an-r-function-from-matlab
%A couple of things if you are trying to run an R script using system(),
%that could be helpful to other fellow MATLABians out there.
%You might need to add Rscript to your system's PATH environment variables.
%If you are on Windows, to do that, in your Windows search box type
%‘Edit the system environment variables’ > click on the button
%‘Environment Variables…’ > click on the line reading ‘Path’,
%and then on button ‘Edit’ > on the window that opens, press ‘New’
%and type in the path of the Rscript.exe:  ‘C:\Program Files\R\R-4.3.3\bin\x64’.
%By default R will set your MATLAB’s pwd as the R’s working directory when
%calling the R script via MATLAB. If the R script you want to run is located elsewhere,
%simply type the full path. For example:

system('Rscript C:\Users\x\Documents\myscript.R')

%Also, you will need to enclose the path inside double quotation marks
%if you are specifying a path that contains spaces (-thats a cmd thing):

system('Rscript "C:\Users\x\Documents\Another Folder\myscript.R"')

