%% Archetype_Computer test script.  Demonstrates how users should format their data.

% This script computes the archetype scores for the GTEx data
% The test output is saved in /Outputs/Test_Archetype_Scores.csv

% General use assumes the data to be TPM normalized (and not log-transformed)

% Generalization to other normalization schemes requires retraining the
% archetypes with the provided data.

Data=readtable('GTEx_by_tissue.txt');
Entrez_Test=Data.Name;
Data=Data(:,3:end);
Data=table2array(Data);
Entrez_Test=extractBefore(Entrez_Test,'.');

% Load the Entrez of the archetype genes
% The associated gene names are given in /Outputs/Genes.csv for reference

Entrez=readcell('Outputs/Entrez.csv');

[~,ia,ib]=intersect(Entrez,Entrez_Test); % Aligns the genes according to Entrez

% Users should make sure that they have all of the required (780) genes!

Data=Data(ib,:);
Data(ia,:)=Data;

H=Archetype_Computer(Data);

% If run correctly, this should be identical to
% /Outputs/Test_Archetype_Scores.csv

