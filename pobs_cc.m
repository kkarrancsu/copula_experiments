function [ U ] = pobs_cc( X )
%PSEUDOOBS Generates Pseudo-Observations from a given multivariate vector
%          See https://en.wikipedia.org/wiki/Copula_(probability_theory)#Empirical_copulas
%          for details
% Inputs:
%  X - an M x D multivariate input matrix, where M is the # of samples of
%      the multivariate joint distribution, and D is the dimensionality
% Optional Inputs:
%  varargin{1} - a string indicating which method to use.  By default, we
%                use ranks to generate the pseudo-observations, but an
%                empirical CDF method could also be used.  Valid values of
%                varargin{1} are: a.) rank
%                                 b.) ecdf
%  varargin{2} - if varargin{1} is ecdf, then varargin{2} must specify the
%                number of points to use in the ECDF
% Outputs:
%  U - an M x D matrix of the pseudo-observations
%
%**************************************************************************
%* 
%* Copyright (C) 2016  Kiran Karra <kiran.karra@gmail.com>
%*
%* This program is free software: you can redistribute it and/or modify
%* it under the terms of the GNU General Public License as published by
%* the Free Software Foundation, either version 3 of the License, or
%* (at your option) any later version.
%*
%* This program is distributed in the hope that it will be useful,
%* but WITHOUT ANY WARRANTY; without even the implied warranty of
%* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%* GNU General Public License for more details.
%*
%* You should have received a copy of the GNU General Public License
%* along with this program.  If not, see <http://www.gnu.org/licenses/>.
%* 
%**************************************************************************

M = size(X);
U = tiedrank(X)/(M+1);      % scale by M+1 to mitigate boundary errors

end



