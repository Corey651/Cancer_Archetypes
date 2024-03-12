function H=Archetype_Computer(Data)

% Function for computing archetype scores on new data

W=readmatrix('Outputs/Archetypes.csv');
S=readmatrix('Outputs/S.csv');
DD=Data./S;
DD=DD./sum(DD,1);
H=NMF_New_Weights(DD,W,6);
