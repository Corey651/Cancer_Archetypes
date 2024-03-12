%% Loads Data

Ann_Gene=readtable('genesets.v7.5.1.txt');

Apop=split(Ann_Gene{19,2},',');
DNA_Repair=split(Ann_Gene{39,2},',');
Glycolysis=split(Ann_Gene{59,2},',');
Hypoxia=split(Ann_Gene{79,2},',');
OxPhos=split(Ann_Gene{99,2},',');

clear Ann_Gene

Genes_Ann=unique([Apop;DNA_Repair;Glycolysis;Hypoxia;OxPhos]);


Data=readtable('GTEx_by_tissue.txt');
Genes=Data.Description;
Entrez=Data.Name;
Cell_Type=Data.Properties.VariableNames;
Cell_Type=Cell_Type(3:end);
Data=Data(:,3:end);
Data=table2array(Data);

S=std(Data,0,2);
Ind=find(S~=0);
Data=Data(Ind,:);
Genes=Genes(Ind);
Entrez=Entrez(Ind);

[C,ia,ib]=intersect(string(Genes),string(Genes_Ann));
Data_Ann=Data(ia,:);
Genes_Ann=Genes_Ann(ib);
Entrez=Entrez(ia);
Entrez=extractBefore(Entrez,'.');
%writecell(Entrez,'Entrez.csv')
%writecell(Genes_Ann,'Genes.csv')

Phen_Table=zeros(length(Genes_Ann),5);

Phen_Table(ismember(string(Genes_Ann),string(Apop)),1)=1;
Phen_Table(ismember(string(Genes_Ann),string(DNA_Repair)),2)=1;
Phen_Table(ismember(string(Genes_Ann),string(Glycolysis)),3)=1;
Phen_Table(ismember(string(Genes_Ann),string(Hypoxia)),4)=1;
Phen_Table(ismember(string(Genes_Ann),string(OxPhos)),5)=1;

clear Apop DNA_Repair Glycolysis Hypoxia OxPhos ia ib C Ind S

%% Set up data and Compute relative expression across each pathway in each cell type

S=std(Data_Ann,0,2);

%writematrix(S,'S.csv')
DD=Data_Ann./S;

DD=DD./sum(DD,1);


R=DD'*Phen_Table;


%% Fig SI1

% Automatic dimensionality selection from the scree plot ... Zhu, Ghodsi,
% 2005

[u,s,v]=svd(DD);
s=diag(s);
s=s(2:30);
EE=zeros(1,length(s));
for i=1:length(EE)
    m1=mean(s(1:i));
    m2=mean(s(i+1:end));
    V=(i-1)*var(s(1:i))+(length(s)-i-1)*var(s(i+1:end));
    V=V/(length(s)-2);
    F=sum((s(1:i)-m1).^2)+sum((s(i+1:end)-m2).^2);
    F=-F/(2*V)-.5*length(s)*log(V);
    EE(i)=F;
end

figure

hold on
plot(2:length(EE)+1,EE,'LineWidth',4)
xlabel('Number of Components')
ylabel('Profile Log Likelihood')
axis tight
plot([6,6],[108.9396,EE(5)],'k--','LineWidth',2)
hold off



%% N-NMF

[W6,H6]=NMF(DD,6);


T6=table('Size',[54 7],'VariableTypes',{'string','double','double','double','double','double','double'},'VariableNames',{'Cell','Factor1','Factor2','Factor3','Factor4','Factor5','Factor6'});
T6.Cell=string(Cell_Type');
T6{:,2:7}=H6';

G6=table('Size',[length(Genes_Ann) 12],'VariableTypes',{'string','double','double','double','double','double','double','double','double','double','double','double'},'VariableNames',{'Gene','Factor1','Factor2','Factor3','Factor4','Factor5','Factor6','Apop','DNA_Repair','Glycolysis','Hypoxia','OxPhos'});
G6.Gene=string(Genes_Ann);
G6{:,2:7}=W6;
G6{:,8:12}=Phen_Table;

% To save the archetype information

%writetable(G6,'Gene_Weights.csv')
%writetable(T6,'GTEx_Weights.csv')
%writematrix(W6,'Archetypes.csv')


%% Correlation between each trained gene and each archetype

PC=corr([H6;Data_Ann]','Type','Spearman');
PC=PC(1:6,7:end);
Gene_Tab=table('Size',[length(Genes_Ann) 7],'VariableTypes',{'string','double','double','double','double','double','double'},'VariableNames',{'Gene','Factor1','Factor2','Factor3','Factor4','Factor5','Factor6'});
Gene_Tab.Gene=string(Genes_Ann);
Gene_Tab{:,2:7}=PC';

PCC=PC;
for i=1:780
    for j=1:6
        Ind=setdiff(1:6,j);
        PCC(j,i)=PC(j,i)-max(PC(Ind,i));
    end
end
Gene_Tab_Check=table('Size',[length(Genes_Ann) 7],'VariableTypes',{'string','double','double','double','double','double','double'},'VariableNames',{'Gene','Factor1','Factor2','Factor3','Factor4','Factor5','Factor6'});
Gene_Tab_Check.Gene=string(Genes_Ann);
Gene_Tab_Check{:,2:7}=PCC';


%% Figure Computation


HH=H6([2,3,1,5,4,6],:);



Direc=[0,1,2,3,4,5];
Direc1=cos(pi/2-2*pi.*Direc/6);
Direc2=sin(pi/2-2*pi.*Direc/6);
Coor1=Direc1*HH;
Coor2=Direc2*HH;


%% Correlations between pathways and components


B=corrcoef([R,H6']);
B6=B(1:5,6:11);

Corr_Path_Tab=table('Size',[6 6],'VariableTypes',{'string','double','double','double','double','double'},'VariableNames',{'Archetype','Apoptosis','DNA Repair','Glycolysis','Hypoxia','OxPhos'});
Corr_Path_Tab.Archetype=['1';'2';'3';'4';'5';'6'];
Corr_Path_Tab{:,2:6}=B6';
Spearman_6=zeros(5,6);
p_6=zeros(5,6);
for i=1:5
    for j=1:6
        [rho,pval]=corr(R(:,i),H6(j,:)','Type','Spearman');
        Spearman_6(i,j)=rho;
        p_6(i,j)=pval;
    end
end

%% Benjamini-Hochberg - find largest p_6 with qq6 negative

pp6=reshape(p_6,[5*6 1]);
pp6=sort(pp6);
qq6=pp6'-((1:30)/30)*.05;  % significance threshold of .05
p_sig=pp6(max(find(qq6<0)));


%% Fig 1B

figure 

hold on
h=nsidedpoly(6);
h.Vertices=h.Vertices*[cos(2*pi/12) -sin(2*pi/12);sin(2*pi/12) cos(2*pi/12)];
plot(h)
scatter(Coor1,Coor2,'filled')
text(Coor1(4),Coor2(4),Cell_Type{4})
text(Coor1(8),Coor2(8),Cell_Type{8})
text(Coor1(11),Coor2(11),Cell_Type{11})
text(Coor1(25),Coor2(25),Cell_Type{25})
text(Coor1(32),Coor2(32),Cell_Type{32})
text(Coor1(34),Coor2(34),Cell_Type{34})
text(Coor1(36),Coor2(36),Cell_Type{36})
text(Coor1(37),Coor2(37),Cell_Type{37})
text(Coor1(50),Coor2(50),Cell_Type{50})
text(Coor1(54),Coor2(54),Cell_Type{54})
hold off

%% Fig 1C

figure  

map=[repmat([0.2157 0.4941 0.7215],4,1);[.95 .95 .95];[.95 .95 .95]; repmat([0.8941 0.1020 0.1098],4,1)];
heatmap(string(Corr_Path_Tab.Properties.VariableNames(2:6)),{'1','2','3','4','5','6'},(p_6<=p_sig)'.*sign(Spearman_6)','Colormap',map,'CellLabelColor','none')
xlabel('Pathway')
ylabel('Archetype')
axs=struct(gca);
axs.Colorbar.Label.String='Significant Trends';



%% Fig 1D


All_Gene=readtable('genesets.v2023.2.Hs.txt');
Pathway_Name=strings(1,50);
Pathway_Name(1)="Adipogenesis";
Pathway_Name(2)="Allograft Rejection";
Pathway_Name(3)="Androgen";
Pathway_Name(4)="Angiogenesis";
Pathway_Name(5)="Apical Junction";
Pathway_Name(6)="Apical Surface";
Pathway_Name(7)="Apoptosis";
Pathway_Name(8)="Bile Acid Metabolism";
Pathway_Name(9)="Cholesterol";
Pathway_Name(10)="Coagulation";
Pathway_Name(11)="Complement";
Pathway_Name(12)="DNA Repair";
Pathway_Name(13)="E2F Targets";
Pathway_Name(14)="EMT";
Pathway_Name(15)="Estrogen Early";
Pathway_Name(16)="Estrogen Late";
Pathway_Name(17)="Fatty Acid Metabolism";
Pathway_Name(18)="G2M Checkpoint";
Pathway_Name(19)="Glycolysis";
Pathway_Name(20)="HEDGEHOG Signaling";
Pathway_Name(21)="Heme Metabolism";
Pathway_Name(22)="Hypoxia";
Pathway_Name(23)="IL2/STAT5 Signaling";
Pathway_Name(24)="IL6/JAK/STAT3 Signaling";
Pathway_Name(25)="Inflammation";
Pathway_Name(26)="Interferon Alpha";
Pathway_Name(27)="Interferon Gamma";
Pathway_Name(28)="KRAS Signaling Down";
Pathway_Name(29)="KRAS Signaling Up";
Pathway_Name(30)="Mitotic Spindle";
Pathway_Name(31)="MTORC1 Signaling";
Pathway_Name(32)="MYC Targets V1";
Pathway_Name(33)="MYC Targets V2";
Pathway_Name(34)="Myogenesis";
Pathway_Name(35)="NOTCH Signaling";
Pathway_Name(36)="Oxidative Phosphorylation";
Pathway_Name(37)="P53";
Pathway_Name(38)="Pancreas Beta Cells";
Pathway_Name(39)="Peroxisome";
Pathway_Name(40)="PI3K/AKT/MTOR Signaling";
Pathway_Name(41)="Protein Secretion";
Pathway_Name(42)="Reactive Oxygen Species";
Pathway_Name(43)="Spermatogenesis";
Pathway_Name(44)="TGF Beta Signaling";
Pathway_Name(45)="TNFA/NFKB Signaling";
Pathway_Name(46)="Unfolded Protein Response";
Pathway_Name(47)="UV Response Down";
Pathway_Name(48)="UV Response Up";
Pathway_Name(49)="WNT/Beta Catenin Signaling";
Pathway_Name(50)="Drug Metabolism";

Genes_All=split(All_Gene{17,2},',');
for i=2:50
    Genes_All=[Genes_All;split(All_Gene{17+19*(i-1),2},',')];
end
Genes_All=unique(Genes_All);
Genes_All=Genes_All(2:end);
Phen_Table_All=zeros(length(Genes_All),50);
for i=1:50
    Phen_Table_All(ismember(string(Genes_All),string(split(All_Gene{17+19*(i-1),2},','))),i)=1;
end


[~,ia,ib]=intersect(string(Genes),string(Genes_All));
Data_All=Data(ia,:);
Genes_All=Genes_All(ib);
Phen_Table_All=Phen_Table_All(ib,:);

S_All=std(Data_All,0,2);
DD_All=Data_All;
DD_All=Data_All./S_All;
DD_All=DD_All./sum(DD_All,1);

R=DD_All'*Phen_Table_All;

% Correlations between pathways and components

Spearman_6=zeros(50,6);
p_6=zeros(50,6);
for i=1:50
    for j=1:6
        [rho,pval]=corr(R(:,i),H6(j,:)','Type','Spearman');
        Spearman_6(i,j)=rho;
        p_6(i,j)=pval;
    end
end


CCCC=corrcoef(Spearman_6');
[aa,bb]=Clusterfunc(CCCC,10);

lll=1:50;
figure

mappp=zeros(1002,3);
for i=1:501
    mappp(i,:)=[0.2157,0.4941,0.7215]*(1-(i-1)/500)+[.95,.95,.95]*(i-1)/500;
    mappp(i+501,:)=[.95,.95,.95]*(1-(i-1)/500)+[0.8941,0.1020,0.1098]*(i-1)/500;
end
heatmap(Pathway_Name(aa(lll)),{'Survival','Proliferation','Fibroblastic','Energy','Biomass','Senescence'},(Spearman_6(aa(lll),:))','Colormap',mappp,'CellLabelColor','none')
xlabel('Pathway')
ylabel('Archetype')
axs=struct(gca);
axs.Colorbar.Label.String='Spearman Correlation';



 % Fig 2B
 
Catlab(1)="CNS";
Catlab(2)="PNS";
Catlab(3)="Adrenal Cortex";
Catlab(4)="Cervix";
Catlab(5)="Skin";
Catlab(6)="Thyroid";
Catlab(7)="Prostate";
Catlab(8)="Colorectal";
Catlab(9)="Breast";
Catlab(10)="Uterus";
Catlab(11)="Ovary";
Catlab(12)="Gastric";
Catlab(13)="Pancreas";
Catlab(14)="Liver";
Catlab(15)="Kidney";
Catlab(16)="Esophagus";
Catlab(17)="Lung";
Catlab(18)="Blood";

HQ=H6(:,[8,40,3,24,45,51,44,26,21,52,41,49,42,36,34,28,37,54]);

mappp=zeros(1002,3);
for i=1:501
    mappp(i,:)=[0.2157,0.4941,0.7215]*(1-(i-1)/500)+[.95,.95,.95]*(i-1)/500;
    mappp(i+501,:)=[.95,.95,.95]*(1-(i-1)/500)+[0.8941,0.1020,0.1098]*(i-1)/500;
end
figure

heatmap(Catlab,{"Survival","Proliferation","Fibroblastic","Energy","Biomass","Senescence"},HQ,'Colormap',mappp)
caxis([0 1])



%% Clear old variables

clear B B6 Cell_Type Coor1 Coor2 Data_Ann DD DDD Direc Direc1 Direc2 EE F H6 HH i j m1 m2 p_6 pval Phen_Table R rho s Spearman_6 u v V

%% Set up with CCLE

fid=fopen('CCLE_expression.csv');
CCLE_Data = textscan(fid, ['%s', repmat('%f',1,19221)], 'collectoutput',true,'headerlines',1,'delimiter',',');
fclose(fid);

fid=fopen('CCLE_expression.csv');
CCLE_Genes=textscan(fid,'%s',19222,'Delimiter',',');
fclose(fid);
CCLE_Genes=CCLE_Genes{1};
CCLE_Genes=CCLE_Genes(2:end);
for i=1:length(CCLE_Genes)
    CCLE_Genes{i}=CCLE_Genes{i}(1:min(find(isspace(CCLE_Genes{i})))-1);
end
CCLE_CL=CCLE_Data{1};
CCLE_Data=CCLE_Data{2};

C_Genes_Test=CCLE_Genes;
C_Data_Test=CCLE_Data;

[~,ia,ib]=intersect(CCLE_Genes,Genes_Ann);
CCLE_Data=CCLE_Data(:,ia);
% puts the data in the same order as Genes_Ann from GTEx
CCLE_Data(:,ib)=CCLE_Data;

%% Order Annotation table and Expression Data

CCLE_Anno=readtable('sample_info_CCLE.csv');
CCLE_ID=CCLE_Anno.DepMap_ID;
[~,ia,ib]=intersect(CCLE_ID,CCLE_CL);
CCLE_Data=CCLE_Data(ib,:);
CCLE_Anno=CCLE_Anno(ia,:);
C_Data_Test=C_Data_Test(ib,:);

%% Computation of NMF coordinates

CCLE_Data=CCLE_Data';
CCLE_Data=2.^CCLE_Data-1;
DD=CCLE_Data./S;

DD=DD./sum(DD,1);


H=NMF_New_Weights(DD,W6,6);

CCLE_All_Tab=table('Size',[1405 9],'VariableTypes',{'string','string','string','double','double','double','double','double','double'},'VariableNames',{'Cell_Line','Tissue','Cancer_Subtype','1','2','3','4','5','6'});
CCLE_All_Tab.Cell_Line=CCLE_Anno.CCLE_Name;
CCLE_All_Tab.Tissue=CCLE_Anno.lineage;
CCLE_All_Tab.Cancer_Subtype=CCLE_Anno.Subtype;
CCLE_All_Tab{:,4:9}=H';
%writetable(CCLE_All_Tab,'CCLE_All_Tab.csv')


%% Clear old variables
clear ans CCLE_CL CCLE_Genes CCLE_ID Corr_Path_Tab Data fid G6 Genes i ia ib T6


%% Fig 2A

HH_canc=H([2,3,1,5,4,6],:);
Direc=[0,1,2,3,4,5];
Direc1=cos(pi/2-2*pi.*Direc/6);
Direc2=sin(pi/2-2*pi.*Direc/6);
Coor1_canc=Direc1*HH_canc;
Coor2_canc=Direc2*HH_canc;

Lins=unique(CCLE_Anno.lineage);


% Brain = 6
% Breast = 5
% Colon =  8
% Liver = 16
% Lung = 17
% Pancreas = 20
% Blood = (3,18,22)


clors=[1 0 0; 0 1 0;0 0 1; 0 1 1; 1 0 1; 1 1 0; 0 0 0];

figure 
hold on
h=nsidedpoly(6);
h.Vertices=h.Vertices*[cos(2*pi/12) -sin(2*pi/12);sin(2*pi/12) cos(2*pi/12)];
plot(h)

Indol=find(strcmp(CCLE_Anno.lineage,Lins{6}));
j=convhull(Coor1_canc(Indol)',Coor2_canc(Indol)');
hh1=plot(Coor1_canc(Indol(j)),Coor2_canc(Indol(j)),'Color',clors(1,:),'LineWidth',2)

Indol=find(strcmp(CCLE_Anno.lineage,Lins{5}));
j=convhull(Coor1_canc(Indol)',Coor2_canc(Indol)');
hh2=plot(Coor1_canc(Indol(j)),Coor2_canc(Indol(j)),'Color',clors(2,:),'LineWidth',2)

Indol=find(strcmp(CCLE_Anno.lineage,Lins{8}));
j=convhull(Coor1_canc(Indol)',Coor2_canc(Indol)');
hh3=plot(Coor1_canc(Indol(j)),Coor2_canc(Indol(j)),'Color',clors(3,:),'LineWidth',2)

Indol=find(strcmp(CCLE_Anno.lineage,Lins{16}));
j=convhull(Coor1_canc(Indol)',Coor2_canc(Indol)');
hh4=plot(Coor1_canc(Indol(j)),Coor2_canc(Indol(j)),'Color',clors(4,:),'LineWidth',2)

Indol=find(strcmp(CCLE_Anno.lineage,Lins{17}));
j=convhull(Coor1_canc(Indol)',Coor2_canc(Indol)');
hh5=plot(Coor1_canc(Indol(j)),Coor2_canc(Indol(j)),'Color',clors(5,:),'LineWidth',2)

Indol=find(strcmp(CCLE_Anno.lineage,Lins{20}));
j=convhull(Coor1_canc(Indol)',Coor2_canc(Indol)');
hh6=plot(Coor1_canc(Indol(j)),Coor2_canc(Indol(j)),'Color',clors(6,:),'LineWidth',2)

Indol=find(ismember(CCLE_Anno.lineage,Lins([3,18,22])));
j=convhull(Coor1_canc(Indol)',Coor2_canc(Indol)');
hh7=plot(Coor1_canc(Indol(j)),Coor2_canc(Indol(j)),'Color',clors(7,:),'LineWidth',2)

legend([hh1 hh2 hh3 hh4 hh5 hh6 hh7],{'Brain','Breast','Colon','Liver','Lung','Pancreas','Blood'})

hold off

%% Fig 2C

figure

Cats=categorical(CCLE_Anno.lineage);
meancat1 = groupsummary(H(1,:)',Cats,'mean');
meancat2 = groupsummary(H(2,:)',Cats,'mean');
meancat3 = groupsummary(H(3,:)',Cats,'mean');
meancat4 = groupsummary(H(4,:)',Cats,'mean');
meancat5 = groupsummary(H(5,:)',Cats,'mean');
meancat6 = groupsummary(H(6,:)',Cats,'mean');
Catts=unique(Cats);

MM=[meancat1,meancat2,meancat3,meancat4,meancat5,meancat6];
MC=mean(MM);
SC=std(MM);
MMM=(MM-MC)./SC;
CCCC=corrcoef(MMM');
[a,b]=Clusterfunc(CCCC,10);
Catlab=Catts(a);
Catlab=strrep(string(Catlab),"_"," ");

Catlab(1)="Fibroblast";
Catlab(2)="CNS";
Catlab(3)="PNS";
Catlab(4)="Bone";
Catlab(5)="Embryo";
Catlab(6)="Adrenal Cortex";
Catlab(7)="Cervix";
Catlab(8)="Unknown";
Catlab(9)="Eye";
Catlab(10)="Soft Tissue";
Catlab(11)="Skin";
Catlab(12)="Epidermis";
Catlab(13)="Thyroid";
Catlab(14)="Upper Airway";
Catlab(15)="Prostate";
Catlab(16)="Bile Duct";
Catlab(17)="Colorectal";
Catlab(18)="Breast";
Catlab(19)="Uterus";
Catlab(20)="Ovary";
Catlab(21)="Gastric";
Catlab(22)="Pancreas";
Catlab(23)="Liver";
Catlab(24)="Kidney";
Catlab(25)="Esophagus";
Catlab(26)="Urinary Tract";
Catlab(27)="Lung";
Catlab(28)="Blood";
Catlab(29)="Lymphocyte";
Catlab(30)="Plasma Cell";



mappp=zeros(1002,3);
for i=1:501
    mappp(i,:)=[0.2157,0.4941,0.7215]*(1-(i-1)/500)+[.95,.95,.95]*(i-1)/500;
    mappp(i+501,:)=[.95,.95,.95]*(1-(i-1)/500)+[0.8941,0.1020,0.1098]*(i-1)/500;
end
heatmap(Catlab,{"Survival","Proliferation","Fibroblastic","Energy","Biomass","Senescence"},MMM(a,:)','Colormap',mappp)
caxis([-1.5 1.5])



%% CCLE Drug responses

Drug_Dat=readtable('CCLE_Drugs.csv');

CQ=CCLE_Anno.CCLE_Name;
CR=Drug_Dat.CCLECellLineName;
[~,ia,ib]=intersect(CQ,CR);
Drug_CL_Names=CQ(ia);
CCLE_Anno_Drug=CCLE_Anno(ia,:);
H_Drug=H(:,ia);
Drugs=unique(Drug_Dat.Compound);
Drug_Response=10000*ones(length(Drug_CL_Names),length(Drugs)); % choosing some arbitrary number instead of NA
for i=1:length(Drug_CL_Names)
    for j=1:length(Drugs)
        Ind_Drug=find(strcmp(Drug_Dat.CCLECellLineName,Drug_CL_Names(i))&strcmp(Drug_Dat.Compound,Drugs(j)));
        if ~isempty(Ind_Drug)
            Drug_Response(i,j)=Drug_Dat.ActArea(Ind_Drug(1));
        end
    end
end

Drug_Target=Drugs;
for i=1:length(Drug_Target)
    Drug_Target(i)=Drug_Dat.Target(min(find(strcmp(Drug_Dat.Compound,Drugs(i)))));
end

%Drug Sensitivities



Drug_R=Drug_Response;
Drug_R(Drug_R>9999)=0;

%% Irinotecan  Fig 3D

% Show Blood, Breast, Lung, Others
figure

tiledlayout(1,2)
nexttile

hold on

Ind=find(strcmp(CCLE_Anno_Drug.lineage,"breast"));
scatter(H_Drug(2,Ind),Drug_R(Ind,6),50,'d','MarkerEdgeColor',[.9290,.6940,.1250],'LineWidth',1.5)
C=setdiff(1:length(H_Drug(2,:)),Ind);

Ind=find(strcmp(CCLE_Anno_Drug.lineage,"blood"));
scatter(H_Drug(2,Ind),Drug_R(Ind,6),50,'square','MarkerEdgeColor',[.6350,.0780,.1840],'LineWidth',1.5)
C=setdiff(C,Ind);

Ind=find(strcmp(CCLE_Anno_Drug.lineage,"lung"));
scatter(H_Drug(2,Ind),Drug_R(Ind,6),50,'+','MarkerEdgeColor',[0,.4470,.7410],'LineWidth',1.5)
C=setdiff(C,Ind);

scatter(H_Drug(2,C),Drug_R(C,6),30,'MarkerEdgeColor',[.5,.5,.5])


legend('Breast','Blood','Lung','Others')
xlabel('Proliferation Archetype Score')
ylabel('Irinotecan Activity Area')
hold off

%% Lapatinib Fig 3E

nexttile 
hold on

Ind=find(strcmp(CCLE_Anno_Drug.lineage,"breast"));
scatter(H_Drug(5,Ind),Drug_R(Ind,9),50,'d','MarkerEdgeColor',[.9290,.6940,.1250],'LineWidth',1.5)
C=setdiff(1:length(H_Drug(5,:)),Ind);

Ind=find(strcmp(CCLE_Anno_Drug.lineage,"colorectal"));
scatter(H_Drug(5,Ind),Drug_R(Ind,9),50,'square','MarkerEdgeColor',[.6350,.0780,.1840],'LineWidth',1.5)
C=setdiff(C,Ind);

Ind=find(strcmp(CCLE_Anno_Drug.lineage,"lung"));
scatter(H_Drug(5,Ind),Drug_R(Ind,9),50,'+','MarkerEdgeColor',[0,.4470,.7410],'LineWidth',1.5)
C=setdiff(C,Ind);

scatter(H_Drug(5,C),Drug_R(C,9),30,'MarkerEdgeColor',[.5,.5,.5])

legend('Breast','Colorectal','Lung','Others')
xlabel('Biomass Archetype Score')
ylabel('Lapatinib Activity Area')
hold off

%% Clear old variables

clear CCLE_Data CQ CR Drug_CL_Names i ia ib Ind_Drug j

%% Correlations between archetype and drugs

Corr_Drug=zeros(6,length(Drugs));
Spearman_Drug=zeros(6,length(Drugs));
p_Drug=zeros(6,length(Drugs));

for j=1:length(Drugs)
    Ind_Drug=find(Drug_Response(:,j)<10000);
    HD=H_Drug(:,Ind_Drug);
    DR=Drug_Response(Ind_Drug,j);
    for i=1:6
        CD=corrcoef(HD(i,:),DR);
        Corr_Drug(i,j)=CD(1,2);
        [rho,pval]=corr(HD(i,:)',DR,'Type','Spearman');
        Spearman_Drug(i,j)=rho;
        p_Drug(i,j)=pval;
    end
end

%% Banjamini-Hochberg

pp=reshape(p_Drug,[6*24 1]);
pp=sort(pp);
qq=pp'-((1:144)/144)*.05;  % significance threshold of .05
p_sig=pp(max(find(qq<0)));


%% Fig 3A
figure
CCCC=corrcoef(Corr_Drug);
[a,b]=Clusterfunc(CCCC,10);
aa=a([1,2,21,5,3,4,9,7,8,11,13,6,10,15,12,14,16,17,18,19,20,23,22,24]);
a=aa;
heatmap(Drugs(a),{"1","2","3","4","5","6"},sign(Spearman_Drug(:,a)).*(p_Drug(:,a)<=p_sig),'Colormap',map,'CellLabelColor','none')

Tab_Lab=['Cell_Line','Tissue','Cancer_Subtype','1','2','3','4','5','6',string(Drugs(a))'];
CCLE_Drug_Tab=table('Size',[469 33],'VariableTypes',['string','string','string',string(repmat({'double'},1,30))],'VariableNames',Tab_Lab);
CCLE_Drug_Tab.Cell_Line=CCLE_Anno_Drug.CCLE_Name;
CCLE_Drug_Tab.Tissue=CCLE_Anno_Drug.lineage;
CCLE_Drug_Tab.Cancer_Subtype=CCLE_Anno_Drug.Subtype;
CCLE_Drug_Tab{:,4:9}=H_Drug';
CCLE_Drug_Tab{:,10:33}=Drug_Response(:,a);
%writetable(CCLE_Drug_Tab,'CCLE_Drug_Tab.csv')

Drug_Table=table('Size',[24 2],'VariableType',{'string','string'},'VariableNames',{'Drug_Name','Drug_Target'});
Drug_Table.Drug_Name=Drugs(a);
Drug_Table.Drug_Target=Drug_Target(a);
%writetable(Drug_Table,'Drug_Info.csv')


%% Fig 3B

T=readtable('GSE155800_GIST_gene_TPM_matrix.txt');
Entrez_G=T.Gene;
Entrez_G=extractBefore(Entrez_G,'.');
DD=T(:,2:end);
DD=table2array(DD);

% Aligning and reordering the genes

[~,ia,ib]=intersect(Entrez,Entrez_G);
DD=DD(ib,:);
DD(ia,:)=DD;
H_GIST=Archetype_Computer(DD);

figure

h1=boxplot([H_GIST(2,6:10)',H_GIST(2,16:20)'],'Labels',{'Sensitive','Resistant'});
ylim([0 .3])
ylabel('Proliferation')
title('GIST')
set(h1,'LineWidth',2)

[h,p]=ttest2(H_GIST(2,6:10)',H_GIST(2,16:20)');

%% Fig 3C

figure

T=readtable('GSE186901_palbo_geo_tpm_mat.txt');

TG=T.Var2;

% Aligining gene names

TG{16063}='ELOA';
TG{5113}='ERO1A';
TG{8619}='LARGE1';
TG{2671}='NOCT';
TG{8910}='PLPPR4';
TG{16814}='TOMM70';
[~,ia,ib]=intersect(Genes_Ann,TG);
DD=T{:,4:end};
DDD=DD(ib,:);
DQ=zeros(780,90);
DQ(ia,:)=DDD;

H_Palbo=Archetype_Computer(DQ);
H_Pairs=H_Palbo(:,[2,3,7,8,9,10,11,12,13,14,15,16,17,18,19,20,36,37,40,41,47,48,52,53,55,56,59,60,61,62,64,65,67,68,77,78,79,80,81,82,83,84,85,86,89,90]);
boxplot([H_Pairs(6,2:2:46)'-H_Pairs(6,1:2:45)'],'Labels',{'Post-Pre'})
[h,p]=ttest(H_Pairs(6,1:2:45)',H_Pairs(6,2:2:46)');
ylabel('Senescence')



%% Clear old variables

clear a aa b CCCC CD DR Drug_Target Drug_Dat Drugs HD i Ind_Drug j pval rho Tab_Lab
clear All_Gene bb C C_Data_Test C_Genes_Test Catlab Cats Catts Coor1_canc Coor2_canc
clear Data_All DDD HH_canc H_Drug H_GIST H_Palbo HQ Indol meancat1 meancat2 meancat3 meancat4
clear meancat5 meancat6 MM MMM p p_Drug p_sig Pathway_Name PC PCC Phen_Table_All pp pp6 qq qq6 S_All
clear SC Spearman_Drug T TG CCLE_Drug_Tab Corr_Drug Drug_R Drug_Response Drug_Table H_Pairs Lins

%% Mutation Biomarkers

Mut_Data=readtable('CCLE_mutations.csv');
Mut_Data=Mut_Data(strcmp(Mut_Data.isTCGAhotspot,'True'),:);
Mut_Data=Mut_Data(Mut_Data.TCGAhsCnt>=5,:);  % Require mutation to be present in at least 5 TCGA samples to minimize false positives
Muts=unique(Mut_Data.Hugo_Symbol);
CQ=CCLE_Anno.DepMap_ID;
CR=Mut_Data.DepMap_ID;
[IDD,ia,ib]=intersect(CQ,CR);
CCLE_Anno_Mut=CCLE_Anno(ia,:);

Mut_Data=Mut_Data(ismember(Mut_Data.DepMap_ID,IDD),:);

Mut_Mat=zeros(length(Muts),length(ia));
for j=1:length(ia)
    Mut_Mat(ismember(Muts,unique(Mut_Data.Hugo_Symbol(string(Mut_Data.DepMap_ID)==IDD(j)))),j)=1;
end

S_Mut=sum(Mut_Mat,2);
Muts=Muts(S_Mut>=10);
Mut_Mat=Mut_Mat(S_Mut>=10,:);
H_Mut=H(:,ia);


%% Make table of significant correlations with mutations for each archetype

Spearman_Mut=zeros(length(Muts),6);
Mut_Table=table('Size',[length(Muts) 7],'VariableTypes',{'string','double','double','double','double','double','double'},'VariableNames',{'Gene','p_1','p_2','p_3','p_4','p_5','p_6'});
p_mut=zeros(73,6);

Mut_Table.Gene=Muts;
for i=1:length(Muts)
    MM=Mut_Mat(i,:);
    for j=1:6
        [rho,pval]=corr(H_Mut(j,:)',MM','Type','Spearman');
        Spearman_Mut(i,j)=rho;
        p_mut(i,j)=pval;
        Mut_Table{i,j+1}=sign(rho)*pval;
    end
end

CCCC=corrcoef(Spearman_Mut');
[a,b]=Clusterfunc(CCCC,100);

%% Benjamini-Hochberg

pp=reshape(p_mut,[73*6 1]);
pp=sort(pp);
qq=pp'-((1:438)/438)*.05;  % significance threshold of .05
p_sig=pp(max(find(qq<0)));

%% Plot significant mutations  Fig 5A

Mut_Ind=find(min(p_mut,[],2)<=p_sig);
aa=[2,8,9,15,16,18,13,14,4,7,5,10,17,6,1,11,12,3];
figure
P_Mat=(p_mut<=p_sig)'.*sign(Spearman_Mut)';
P_Mat=P_Mat(:,Mut_Ind);
map=[repmat([0.2157 0.4941 0.7215],4,1);[.95 .95 .95];[.95 .95 .95]; repmat([0.8941 0.1020 0.1098],4,1)];
heatmap(string(Muts(Mut_Ind(aa))),{'1','2','3','4','5','6'},P_Mat(:,aa),'Colormap',map,'CellLabelColor','none')
xlabel('Mutated Gene')
ylabel('Archetype')
axs=struct(gca);
axs.Colorbar.Label.String='Significant Trends';


%% Fig 5B

figure
hold on
histogram(H_Mut(2,Mut_Mat(58,:)==0),30,'Normalization','probability')
histogram(H_Mut(2,Mut_Mat(58,:)==1),30,'Normalization','probability')
xlabel('Proliferation Archetype Score')
ylabel('Frequency')
legend('TP53 WT','TP53 Mut')
hold off


%% Copy number analysis


CNA_Data=readtable('CCLE_gene_cn.csv');
CNA=table2array(CNA_Data(:,2:end));
CNA=(2.^CNA-1)*2;
CNA=min(CNA,4);
CNA=CNA(:,1:24210); % remove sex chromosomes

CNA_Genes=string(CNA_Data.Properties.VariableNames);
CNA_Genes=CNA_Genes(2:end);
CNA_Genes=CNA_Genes(1:24210);
CNA_Genes=extractBefore(CNA_Genes,'_');

CR=CNA_Data.Var1;
CQ=CCLE_Anno.DepMap_ID;
[~,ia,ib]=intersect(CQ,CR);
CCLE_Anno_CNA=CCLE_Anno(ia,:);
CNA=CNA(ib,:);
H_CNA=H(:,ia);
DD=[CNA];
[Spearman_CNA,p_CNA]=corr(DD,H_CNA','Type','Spearman');

pp=reshape(p_CNA,[24210*6 1]);
pp=sort(pp);
qq=pp'-((1:145260)/145260)*.05;  % significance threshold of .05
p_sig=pp(max(find(qq<0)));
Cmin=max(Spearman_CNA(p_CNA<p_sig&Spearman_CNA<0));
Cmax=min(Spearman_CNA(p_CNA<p_sig&Spearman_CNA>0));


CNA_Table=table('Size',[24210 7],'VariableTypes',{'string','double','double','double','double','double','double'},'VariableNames',{'Gene','p_1','p_2','p_3','p_4','p_5','p_6'});
CNA_Table.Gene=[CNA_Genes'];
CNA_Table{:,2:7}=Spearman_CNA;

Gene_Data=readtable('gene_data.csv');
GeneN_Name=Gene_Data.Var1;
Bin=Gene_Data.bin_id;

CNA_Table_Bin=zeros(6,20000);
Lst=[];
for i=1:20000
    [Q,ira,irb]=intersect(CNA_Genes,GeneN_Name(Bin==i));
    if isempty(Q)==0
        CNA_Table_Bin(:,i)=mean(Spearman_CNA(ira,:),1);
        if isempty(Lst)==0
            CNA_Table_Bin(:,Lst)=repmat(((CNA_Table_Bin(:,i)+CNA_Table_Bin(:,max(min(Lst)-1,1)))/2),1,length(Lst));
            Lst=[];
        end
    else
        Lst=[Lst,i];
    end
end

Impact=readtable('impact.csv');
[ccc,aaa,bbb]=intersect(Impact.HugoSymbol,CNA_Genes);
C_Impact=CNA(:,bbb);

onco=Impact.HugoSymbol(strcmp(Impact.MSK_IMPACT,'Yes'));
Impact_overlap=cell(1,6);
for i=1:6
    R=find(CNA_Table_Bin(i,:)<=quantile(CNA_Table_Bin(i,:),.005)|CNA_Table_Bin(i,:)>=quantile(CNA_Table_Bin(i,:),.995));
    MM=intersect(Gene_Data.hgnc_symbol(ismember(Gene_Data.bin_id,R)),onco);
    Impact_overlap{i}=MM;
end
Imp_ov=Impact_overlap{1};
for i=2:6
    Imp_ov=union(Imp_ov,Impact_overlap{i});
end

Imp_bin=zeros(1,length(Imp_ov));
for i=1:length(Imp_bin)
    Imp_bin(i)=Gene_Data.bin_id(strcmp(Imp_ov(i),Gene_Data.hgnc_symbol));
end

% writematrix(CNA_Table_Bin,'CNA_Table_Bin.txt')
% Fig 5C is then plotted in R, see CNA_Heatmap.R

%% Fig 5D,E,F are the Predicted vs Actual plots of the below analyses

% For all below, use the regression learner application with Data Set
% Variable = Regr... and Response = outp... along with the remaining
% instructions to replicate the figures.



RegrD=[Mut_Mat;ones(1,1250)]';
outpD=H_Mut(2,:);

% Regression learner settings - Model: Interactions Linear, No validation

%

RegrE=[C_Impact,ones(1387,1)];
outpE=H_CNA(2,:);

%Regression learner settings - use PCA, 500 PCs, Model: Linear, No
%validation

%

% Load mutation data
Mut_Data=readtable('CCLE_mutations.csv');
Mut_Data=Mut_Data(strcmp(Mut_Data.isTCGAhotspot,'True'),:);
Mut_Data=Mut_Data(Mut_Data.TCGAhsCnt>=5,:);  % Require mutation to be present in at least 5 TCGA samples to minimize false positives
Muts=unique(Mut_Data.Hugo_Symbol);
CQ=CCLE_Anno.DepMap_ID;
CR=Mut_Data.DepMap_ID;
[IDD,ia,ib]=intersect(CQ,CR);
CCLE_Anno_Mut=CCLE_Anno(ia,:);

Mut_Data=Mut_Data(ismember(Mut_Data.DepMap_ID,IDD),:);

Mut_Mat=zeros(length(Muts),length(ia));
for j=1:length(ia)
    Mut_Mat(ismember(Muts,unique(Mut_Data.Hugo_Symbol(string(Mut_Data.DepMap_ID)==IDD(j)))),j)=1;
end

S_Mut=sum(Mut_Mat,2);
Mut_Mat=Mut_Mat(S_Mut>=10,:);

H_Mut=H(:,ia);

% Load CNA Data

CNA_Data=readtable('CCLE_gene_cn.csv');
CNA=table2array(CNA_Data(:,2:end));
CNA=(2.^CNA-1)*2;
CNA=min(CNA,4);
CNA=CNA(:,1:24210); % remove sex chromosomes

CNA_Genes=string(CNA_Data.Properties.VariableNames);
CNA_Genes=CNA_Genes(2:end);
CNA_Genes=CNA_Genes(1:24210);
CNA_Genes=extractBefore(CNA_Genes,'_');

CRR=CNA_Data.Var1;
[~,iaa,ibb]=intersect(IDD,CRR);

CNA=CNA(ibb,:);
Mut_Mat=Mut_Mat(:,iaa);
Mut_Mat=Mut_Mat';
H_Comb=H_Mut(:,iaa);

Impact=readtable('impact.csv');
[ccc,aaa,bbb]=intersect(Impact.HugoSymbol,CNA_Genes);
C_Impact=CNA(:,bbb);

% Pre-PCA

[uu,ss,vv]=svd(C_Impact);
C_Impact_pca=uu(:,1:500);

RegrF=[C_Impact_pca,Mut_Mat,ones(1237,1)];
outpF=H_Comb(2,:);

% Regression learner settings - Model: Linear, No validation

clear CNA_Table_Bin C_Impact C_Impact_pca CCCC CNA CNA_Data CNA_Genes CNA_Table
clear CNA_Table_Bin DD DD_All Gene_Data H_CNA H_Comb H_Mut Impact Mut_Data Mut_Mat
clear p_CNA p_Mat p_mut pp qq Regr RegrD RegrE RegrF Spearman_CNA ss uu vv Spearman_Mut
clear DQ GeneN_Name CCLE_Anno_Mut CCLE_Anno_Drug CCLE_Anno_CNA Bin Entrez_G

%% Analysis of TCGA Breast

TCGA=readtable('TCGA_BRCA_TPM.csv');
Entrez_TCGA=TCGA.Gene_ensembl;
Entrez_TCGA=extractBefore(Entrez_TCGA,'.');
TCGA_Data=TCGA(:,4:end);
TCGA_Data=table2array(TCGA_Data);
V=sum(TCGA_Data,2);
TCGA_Data=TCGA_Data(V>0,:);
Entrez_TCGA=Entrez_TCGA(V>0);

[~,ia,ib]=intersect(Entrez,Entrez_TCGA);
TCGA_Data=TCGA_Data(ib,:);
TCGA_Data(ia,:)=TCGA_Data;

DD=TCGA_Data./S;
DD=DD./sum(DD,1);

H_TCGA=NMF_New_Weights(DD,W6,6);

%writematrix(H_TCGA,'Breast_Archetype_Scores.csv') % For testing

TCGA_clin=readtable('TCGA_BRCA_clin.csv');
PAM50=TCGA_clin.paper_BRCA_Subtype_PAM50;
death=str2double(TCGA_clin.days_to_death);
followup=(TCGA_clin.days_to_last_follow_up);
death=max([followup,death],[],2);
death=round(death./30);
status=TCGA_clin.paper_vital_status;


ind=find(string(PAM50)=='LumA'|string(PAM50)=='Basal');
indl=find(string(PAM50)=='LumA');
indb=find(string(PAM50)=='Basal');

GG=cell(1111,1);
for i=1:1111
    if H_TCGA(2,i)<=.2
        GG{i}='G1';
    elseif H_TCGA(2,i)<=.4
        GG{i}='G2';
    else
        GG{i}='G3';
    end
end


%% Figs 4A,B


%CCLE

Breast_Subs=find(strcmp(CCLE_Anno.lineage,'breast'));
HB=H(:,Breast_Subs);
CCLE_B_Anno=CCLE_Anno(Breast_Subs,:);
Ind1=find(string(CCLE_B_Anno.lineage_molecular_subtype)=='luminal');
CCLE_B_Anno.lineage_molecular_subtype(Ind1)={'Lum'};
Ind3=find(string(CCLE_B_Anno.lineage_molecular_subtype)=='basal_A');
CCLE_B_Anno.lineage_molecular_subtype(Ind3)={'BasalA'};
Ind4=find(string(CCLE_B_Anno.lineage_molecular_subtype)=='basal_B');
CCLE_B_Anno.lineage_molecular_subtype(Ind4)={'BasalB'};
Ind5=find(string(CCLE_B_Anno.lineage_molecular_subtype)=='HER2_amp');
CCLE_B_Anno.lineage_molecular_subtype(Ind5)={'HER2+'};
Indd=[Ind1;Ind3;Ind4;Ind5];

% TCGA

Ind1=find(string(PAM50)=='Normal');
Ind2=find(string(PAM50)=='LumA');
Ind3=find(string(PAM50)=='LumB');
Ind4=find(string(PAM50)=='Basal');
Ind5=find(string(PAM50)=='Her2');
PAM50(Ind5)={'HER2+'};

Inddd=[Ind1;Ind2;Ind3;Ind4;Ind5];

figure
tiledlayout(2,2)

nexttile
boxplot(HB(2,Indd),CCLE_B_Anno.lineage_molecular_subtype(Indd))
ylim([0 .9])
ylabel('Proliferation Archetype Score')
nexttile
boxplot(H_TCGA(2,Inddd),PAM50(Inddd))
ylim([0 .9])
ylabel('Proliferation Archetype Score')
nexttile
boxplot(HB(5,Indd),CCLE_B_Anno.lineage_molecular_subtype(Indd))
ylim([0 .3])
ylabel('Biomass Archetype Score')
nexttile
boxplot(H_TCGA(5,Inddd),PAM50(Inddd))
ylim([0 .3])
ylabel('Biomass Archetype Score')



%% Fig 4C, SI2A

[p2,fh2,stats2]=MatSurv(death(:), status(:), GG,'Xstep',60,'TimeMax',240,'DispHR',false,'legend',false,'RT_XAxis',false,'GroupsToUse',{{'A_{2}\leq0.2','G1'},{'0.2<A_{2}\leq0.4','G2'},{'A_{2}>0.4','G3'}});

%% Fig 4D

[p,fh,stats]=MatSurv(death(:), status(:),PAM50(:),'Xstep',60,'TimeMax',240,'DispHR',false,'legend',false,'RT_XAxis',false,'GroupsToUse',{{'LumA'},{'LumB'},{'Basal'}});


%% Analysis of TCGA Colon

TCGA=readtable('TCGA_COAD_TPM.csv');
Entrez_TCGA=TCGA.Gene_ensembl;
Entrez_TCGA=extractBefore(Entrez_TCGA,'.');
TCGA_Data=TCGA(:,4:end);
TCGA_Data=table2array(TCGA_Data);
V=sum(TCGA_Data,2);
TCGA_Data=TCGA_Data(V>0,:);
Entrez_TCGA=Entrez_TCGA(V>0);

[~,ia,ib]=intersect(Entrez,Entrez_TCGA);
TCGA_Data=TCGA_Data(ib,:);
TCGA_Data(ia,:)=TCGA_Data;

DD=TCGA_Data./S;
DD=DD./sum(DD,1);

% New archetype computation
H_TCGA=NMF_New_Weights(DD,W6,6);

TCGA_clin=readtable('TCGA_COAD_clin.csv');
death=str2double(TCGA_clin.days_to_death);
followup=(TCGA_clin.days_to_last_follow_up);
death=max([followup,death],[],2);
death=round(death./30);
status=TCGA_clin.vital_status;

GG=cell(481,1);
for i=1:481
    if H_TCGA(4,i)<=.1
        GG{i}='G1';
    else
        GG{i}='G2';
    end
end

% Fig SI2B
[p4,fh4,stats4]=MatSurv(death(:), status(:),GG,'Xstep',30,'TimeMax',120,'DispHR',false,'legend',false,'RT_XAxis',false,'GroupsToUse',{{'A_{4}\leq0.1','G1'},{'A_{4}>0.1','G2'}});

%% Pancreatic

TCGA=readtable('TCGA_PAAD_TPM.csv');
Entrez_TCGA=TCGA.Gene_ensembl;
Entrez_TCGA=extractBefore(Entrez_TCGA,'.');
TCGA_Data=TCGA(:,4:end);
TCGA_Data=table2array(TCGA_Data);
V=sum(TCGA_Data,2);
TCGA_Data=TCGA_Data(V>0,:);
Entrez_TCGA=Entrez_TCGA(V>0);

[~,ia,ib]=intersect(Entrez,Entrez_TCGA);
TCGA_Data=TCGA_Data(ib,:);
TCGA_Data(ia,:)=TCGA_Data;

DD=TCGA_Data./S;
DD=DD./sum(DD,1);

% New archetype computation
H_TCGA=NMF_New_Weights(DD,W6,6);

TCGA_clin=readtable('TCGA_PAAD_clin.csv');
death=(TCGA_clin.days_to_death);
followup=(TCGA_clin.days_to_last_follow_up);
death=max([followup,death],[],2);
death=round(death./30);
status=TCGA_clin.vital_status;

GG=cell(178,1);
for i=1:178
    if H_TCGA(6,i)<=.05
        GG{i}='G1';
    else
        GG{i}='G2';
    end
end

% Fig SI2C
[p5,fh5,stats5]=MatSurv(death(:), status(:),GG,'Xstep',30,'TimeMax',90,'DispHR',false,'legend',false,'RT_XAxis',false,'GroupsToUse',{{'A_{6}\leq0.05','G1'},{'A_{6}>0.05','G2'}});


%% Metastatic breast cancer analysis

T=readtable('QN_TPM_bC_nonlog_median.csv');

Ann_Gene=readtable('genesets.v7.5.1.txt');

Apop=split(Ann_Gene{19,2},',');
DNA_Repair=split(Ann_Gene{39,2},',');
Glycolysis=split(Ann_Gene{59,2},',');
Hypoxia=split(Ann_Gene{79,2},',');
OxPhos=split(Ann_Gene{99,2},',');

clear Ann_Gene

Genes_Ann=unique([Apop;DNA_Repair;Glycolysis;Hypoxia;OxPhos]);


Data=readtable('GTEx_by_tissue.txt');
Genes=Data.Description;
Entrez=Data.Name;
Cell_Type=Data.Properties.VariableNames;
Cell_Type=Cell_Type(3:end);
Data=Data(:,3:end);
Data=table2array(Data);

S=std(Data,0,2);
Ind=find(S~=0);
Data=Data(Ind,:);
Genes=Genes(Ind);
Entrez=Entrez(Ind);

[C,ia,ib]=intersect(string(Genes),string(Genes_Ann));
Data_Ann=Data(ia,:);
Genes_Ann=Genes_Ann(ib);
Entrez=Entrez(ia);
Entrez=extractBefore(Entrez,'.');


Clin=T.Properties.VariableNames;
Clin=Clin(2:42);
[~,ia,ib]=intersect(Entrez,T.Var43);
DD=T{:,2:42};
DDD=DD(ib,:);
DQ=zeros(780,41);
DQ(ia,:)=DDD;
H=Archetype_Computer(DQ);

HQ=H(:,[21,1,9,23,2,13,25,7,16,22,3,12,39,4,11,24,5,8]);

Catlab(1)="Normal Breast";
Catlab(2)="Primary BBC Breast";
Catlab(3)="Primary BBC TCGA";
Catlab(4)="Normal Lung";
Catlab(5)="Lung Metastasis";
Catlab(6)="Primary Lung TCGA";
Catlab(7)="Normal Skin";
Catlab(8)="Skin Metastasis";
Catlab(9)="Primary Skin TCGA";
Catlab(10)="Normal Liver";
Catlab(11)="Liver Metastasis";
Catlab(12)="Primary Liver TCGA";
Catlab(13)="Normal Brain Cortex";
Catlab(14)="Brain Metastasis";
Catlab(15)="Primary Brain TCGA";
Catlab(16)="Normal Adrenal Gland";
Catlab(17)="Adrenal Metastasis";
Catlab(18)="Primary Adrenal TCGA";

mappp=zeros(1002,3);
for i=1:501
    mappp(i,:)=[0.2157,0.4941,0.7215]*(1-(i-1)/500)+[.95,.95,.95]*(i-1)/500;
    mappp(i+501,:)=[.95,.95,.95]*(1-(i-1)/500)+[0.8941,0.1020,0.1098]*(i-1)/500;
end

%% Fig 4F

figure

heatmap(Catlab,{"Survival","Proliferation","Fibroblastic","Energy","Biomass","Senescence"},HQ,'Colormap',mappp)
caxis([0 1])


%% Fig 4E

HHH=HQ([2,3,1,5,4,6],:);

Direc=[0,1,2,3,4,5];
Direc1=cos(pi/2-2*pi.*Direc/6);
Direc2=sin(pi/2-2*pi.*Direc/6);
Coor1=Direc1*HHH;
Coor2=Direc2*HHH;


figure 

hold on
h=nsidedpoly(6);
h.Vertices=h.Vertices*[cos(2*pi/12) -sin(2*pi/12);sin(2*pi/12) cos(2*pi/12)];
plot(h)

% Breast
hh1=scatter(Coor1(1),Coor2(1),'filled','MarkerFaceColor','r');
hh2=scatter(Coor1(2),Coor2(2),'filled','MarkerFaceColor','k');

% Lung
scatter(Coor1(4),Coor2(4),'filled','MarkerFaceColor','r')
scatter(Coor1(5),Coor2(5),'filled','MarkerFaceColor','k')

% Skin
scatter(Coor1(7),Coor2(7),'filled','MarkerFaceColor','r')
scatter(Coor1(8),Coor2(8),'filled','MarkerFaceColor','k')

% Liver
scatter(Coor1(10),Coor2(10),'filled','MarkerFaceColor','r')
scatter(Coor1(11),Coor2(11),'filled','MarkerFaceColor','k')

% Brain
scatter(Coor1(13),Coor2(13),'filled','MarkerFaceColor','r')
scatter(Coor1(14),Coor2(14),'filled','MarkerFaceColor','k')

% Adrenal
scatter(Coor1(16),Coor2(16),'filled','MarkerFaceColor','r')
scatter(Coor1(17),Coor2(17),'filled','MarkerFaceColor','k')

legend([hh1 hh2],{'Normal Tissue','Metastatic BBC'})

hold off



