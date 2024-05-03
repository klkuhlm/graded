% ------------- DREAM with sampling from past and snooker updates: DREAM_ZS --------------------%
%                                                                                               %
% The code presented herein is a Markov Chain Monte Carlo algorithm that runs multiple chains   %
% in parallel for efficient posterior exploration. The algorithm, entitled DREAM_(ZS) is        %
% based on the original DREAM sampling scheme, but uses sampling from an archive of past        %
% states to generate candidate points in each individual chain. Theoy and numerical examples of %
% DREAM_(ZS) have been presented in Vrugt et al. (2009). Details can also be found in           %
% Ter Braak and Vrugt (2008)                                                                    %
%                                                                                               %
% Sampling from past has three main advantages:                                                 %
% (1) Circumvents the requirement of using N = d for posterior exploration. This will speed-up  %
% convergence to a limiting distribution, especially for high-dimensional problems (large d).   %
% (2) Outlier chains do not need explicit consideration. By sampling historical states,         %
% aberrant trajectories an jump directly to the modal region at any time during the             %
% simulation. The N path ways simulated with DREAM_(ZS) therefore maintain detailed balance at  %
% every singe step in the chain.                                                                %
% (3) The transition kernel defining the jumps in each of the chains does not require           %
% information about the current states of the chains. This is of great advantage in a           %
% multi-processor environment where the N candidate points can be generated simultaneously so   %
% that each chain can evolve most efficiently on a different computer. Details of this will be  %
% given in a later publication, which should be ready within the next few months.               %
%                                                                                               %
% DREAM_(ZS) also contains a snooker updater to maximize the diversity of candidate points      %
% and generate jumps beyond parallel direction updates. Finally, DREAM_(ZS) contains subspace   %
% learning in a similar way as DREAM, to maximize the squared jumping distance between two      %
% subsequent points in each chain. This idea has been presented in Vrugt et al. (2008) and      %
% shown to significantly increase the efficiency of posterior exploration. All these options    %
% can be activated from the input file.                                                         %
%                                                                                               %
% DREAM_(ZS) developed by Jasper A. Vrugt and Cajo ter Braak                                    %
%                                                                                               %
% This algorithm has been described in:                                                         %
%                                                                                               %
%   C.J.F. ter Braak, and J.A. Vrugt, Differential Evolution Markov Chain with snooker updater  %
%       and fewer chains, Statistics and Computing, 10.1007/s11222-008-9104-9, 2008             %
%                                                                                               %
%   Vrugt, J.A., and C.J.F. ter Braak, DiffeRential Evolution Adaptive Metropolis with Sampling %
%       from the Past and Subspace Updating, SIAM journal on Optimization                       %
%                                                                                               %
%   Vrugt, J.A., and C.J.F. ter Braak, Multiple Try DiffeRential Evolution Adaptive Metropolis  %
%       for High Performance Computing, SIAM Journal on Distributed Computing                   %
%                                                                                               %
% For more information please read:                                                             %
%                                                                                               %
%   Vrugt J.A., H.V. Gupta, W. Bouten and S. Sorooshian, A Shuffled Complex Evolution           %
%       Metropolis algorithm for optimization and uncertainty assessment of hydrologic model    %
%       parameters, Water Resour. Res., 39 (8), 1201, doi:10.1029/2002WR001642, 2003.           %
%                                                                                               %
%   ter Braak, C.J.F., A Markov Chain Monte Carlo version of the genetic algorithm Differential %
%       Evolution: easy Bayesian computing for real parameter spaces, Stat. Comput., 16,        %
%       239 - 249, doi:10.1007/s11222-006-8769-1, 2006.                                         %
%                                                                                               %
%   Vrugt, J.A., C.J.F. ter Braak, M.P. Clark, J.M. Hyman, and B.A. Robinson, Treatment of      %
%       input uncertainty in hydrologic modeling: Doing hydrology backward using Markov         %
%       chain Monte Carlo, Water Resour. Res., 44, W00B09, doi:10.1029/2007WR006720, 2008.      %
%                                                                                               %
%   Vrugt, J.A., C.J.F. ter Braak, C.G.H. Diks, D. Higdon, B.A. Robinson, and J.M. Hyman,       %
%       Accelerating Markov chain Monte Carlo simulation by self adaptive differential          %
%       evolution with randomized subspace sampling, International Journal of Nonlinear         %
%       Sciences and Numerical Simulation, In Press. 2009.                                      %
%                                                                                               %
% Copyright (c) 2008, Los Alamos National Security, LLC                                         %
%                                                                                               %
% All rights reserved.                                                                          %
%                                                                                               %
% Copyright 2008. Los Alamos National Security, LLC. This software was produced under U.S.      %
% Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory (LANL), which is     %
% operated by Los Alamos National Security, LLC for the U.S. Department of Energy. The U.S.     %
% Government has rights to use, reproduce, and distribute this software.                        %
%                                                                                               %
% NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES A NY WARRANTY, EXPRESS OR  %
% IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE. If software is modified to    %
% produce derivative works, such modified software should be clearly marked, so as not to       %
% confuse it with the version available from LANL.                                              %
%                                                                                               %
% Additionally, redistribution and use in source and binary forms, with or without              %
% modification, are permitted provided that the following conditions are met:                   %
% * Redistributions of source code must retain the above copyright notice, this list of         %
%   conditions and the following disclaimer.                                                    %
% * Redistributions in binary form must reproduce the above copyright notice, this list of      %
%   conditions and the following disclaimer in the documentation and/or other materials         %
%   provided with the distribution.                                                             %
% * Neither the name of Los Alamos National Security, LLC, Los Alamos National Laboratory, LANL %
%   the U.S. Government, nor the names of its contributors may be used to endorse or promote    %
%   products derived from this software without specific prior written permission.              %
%                                                                                               %
% THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS "AS IS" AND   %
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES      %
% OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS %
% ALAMOS NATIONAL SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, %
% SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF   %
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)        %
% HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT %
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,       %
% EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                                            %
%                                                                                               %
%                                                                                               %
% Copyright (c) 2008, Los Alamos National Security, LLC                                         %
%                                                                                               %
% Written by Jasper A. Vrugt: vrugt@lanl.gov                                                    %
%                                                                                               %
% Version 0.5: January 2009                                                                     %
% Version 1.0: April 2011         Maintenance update, explicit treatment of prior distribution  %
% Version 1.1: August 2011        Whittle likelihood function (SPECTRAL ANALYSIS !!)            %
%                                                                                               %
% --------------------------------------------------------------------------------------------- %

% Different test examples from SIAM paper
% example 1: n-dimensional Gaussian distribution
% example 2: multivariate student t distribution
% example 3: n-dimensional banana shaped Gaussian distribution
% example 4: n-dimensional multimodal mixture distribution
% example 5: real-world example using hymod rainfall - runoff model (HYMOD code in MATLAB)
% example 6: real-world example using hymod rainfall - runoff model (HYMOD code in FORTRAN)
% example 7: rainfall-runoff model with generalized log-likelihood function
% example 8: HYDRUS-1D soil hydraulic model: using prior information on soil hydraulic parameters

	
% -------------------------------- Check the following paper ------------------------------ %
%                                                                                           %
%   B. Scharnagl, J.A. Vrugt, H. Vereecken, and M. Herbst (2011), Bayesian inverse          % 
%	modeling of soil water dynamics at the field scale: using prior information             % 
%	on soil hydraulic properties, Hydrology and Earth System Sciences.                      %  
%                                                                                           %
% ----------------------------------------------------------------------------------------- % 

clear;

% Problem specific parameter settings
MCMCPar.n = 6;                          % Dimension of the problem
MCMCPar.ndraw = 400000;                   % Maximum number of function evaluations
MCMCPar.parallelUpdate = 0.9;           % Fraction of parallel direction updates

% Recommended parameter settings
MCMCPar.seq = 3;                        % Number of Markov Chains / sequences
MCMCPar.DEpairs = 1;                    % Number of chain pairs to generate candidate points
MCMCPar.nCR = 3;                        % Number of crossover values used
MCMCPar.m0 = 10 * MCMCPar.n;            % Initial size of Z
MCMCPar.k = 10;                         % Thinning parameter for appending X to Z
MCMCPar.eps = 5e-2;                     % Perturbation for ergodicity
MCMCPar.steps = 100;                     % Number of steps before calculating convergence diagnostics
MCMCPar.pJumpRate_one = 0.20;           % Probability of selecting a jumprate of 1 --> jump between modes

% --------------------------------------------------------------------------------------------
Extra.pCR = 'Update';                   % Adaptive tuning of crossover values
% --------------------------------------------------------------------------------------------

% --------------------------------------- Added for reduced sample storage -------------------
Extra.reduced_sample_collection = 'No'; % Thinned sample collection?
Extra.T = 1000;                         % Every Tth sample is collected
% --------------------------------------------------------------------------------------------

% Define the boundary handling
% None: allows values outside specified range
% Bound: sets values to boundary
% Reflect: 
% Fold: "maintains detailed balance"

%Extra.BoundHandling = 'Reflect';
Extra.BoundHandling = 'Fold';

% Save in memory or not
Extra.save_in_memory = 'Yes';

global idx;
idx = 1;

% Define feasible parameter space (minimum and maximum values)
% k0 is permeability at borehole wall
% n0 is porosity at borehole wall
% eta is power of r/r_w by which porosity drops off with distance
% tau is multiple relating eta and kappa (perm power-law factor)
% m is dimensionality (m=0 cartesian, m=1 cylindrical, m=2 spherical)
% c is formation compressibility 

cm = log10(2.69E-11);

% analytical solution parameters
%	             1      2           3     4       5      6
%                k0     n0         eta    tau     cm      m
ParRange.minn= [-24.0, -4.0,        0.001,  1.0, cm - 4.0, 0.5]; 
ParRange.maxn= [-16.0, log10(0.25), 5.0,  5.0, cm + 4.0,   3.0];

% data tables are time (seconds) and change in pressure (MPa)
data53 = readmatrix('53-1.25r-drawdown.txt');
data54 = readmatrix('54-1.5r-drawdown.txt');
data56 = readmatrix('56-2r-drawdown.txt');
data58 = readmatrix('58-3r-drawdown.txt');

data = [data53; data54; data56; data58]; % concatenate all together

Extra.ndata = [size(data53,1),size(data54,1), ...
               size(data56,1),size(data58,1)]; % size of individual vectors
Extra.time = data(:,1);
Extra.rvec = [1.25, 1.5, 2.0, 3.0];

Measurement.MeasData = data(:,2);
Measurement.N = size(Measurement.MeasData,1);

Measurement.Weights = ones(Measurement.N,1); % unit weight
%Measurement.Weights(int(Measurement.N/2):end) = 0.25;

% Define model name
ModelName = 'powerlaw';

% Define option (Model computes output simulation)
option = 3;

% Indicate the use prior information
%%Extra.InitPopulation = 'PRIOR';
Extra.InitPopulation = 'LHS_BASED';

% No restart -- just running for the first time. Restart can be used if run
% is termined while running; Use Restart = 'Yes'; to reinitialize the code
% from the saved files.
Restart = 'No';

more off;
format compact;

% Run the distributed DREAM algorithm with sampling from past
[Sequences,Reduced_Seq,X,Z,output] = dream_zs(MCMCPar,ParRange,Measurement,ModelName,Extra,option,Restart);

save('powerlaw_stormont_uniform_results.mat','-V7');

