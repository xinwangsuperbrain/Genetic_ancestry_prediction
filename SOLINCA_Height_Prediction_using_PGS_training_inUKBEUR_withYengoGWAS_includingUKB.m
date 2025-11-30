%Using the codes below:
addpath(genpath('/home/dale/matlab'))
addpath(genpath('/space/broce-syn01/1/data/GWAS/PHS/Test_PHS/scripts'))

load('/space/broce-syn01/1/data/GWAS/xinwang/SOLINCA/results_revision/UKB_EUR_train_use_sGWASUKB/UKB_height_betas_iter1_pthreshi7_nfraci5_regtype1_alpha=1.0.mat')
% get genetic with all required SNPs
bfile_format = '/space/broce-syn01/1/data/GWAS/TOPmed/postimputaion/MEGA_HCHS_SOL.MERGED.QC.SNPSonly';
bfile = sprintf(bfile_format);
fileID = fopen(sprintf('%s.bim', bfile));
bim_file = textscan(fileID,'%s %s %s %s %s %s');
fclose(fileID);
snpidvec_sol = bim_file{2};
A1vec = bim_file{5};
A2vec = bim_file{6};
fileID = fopen(sprintf('%s.fam', bfile));
fam_file = textscan(fileID,'%s %s %s %s %s %s');
fclose(fileID);
IID = fam_file{1};
FID = fam_file{2};
nsubj=length(fam_file{1});
[dummy IA_sol IB_ukb] = intersect(snpidvec_sol ,snpids_tmp,'stable');
tic
geno_int8 = PlinkRead_binary2_amd(nsubj, IA_sol, bfile, 1:nsubj); 
geno = single(geno_int8);
toc

snpid_sol = snpidvec_sol(IA_sol); A1vec_sol = A1vec(IA_sol);A2vec_sol = A2vec(IA_sol);
snpids_tmp = snpids_tmp(IB_ukb);
disp([snpid_sol(1:10) snpids_tmp(1:10)])

% load pheno
pheno_all = readtable("/space/broce-syn01/1/data/GWAS/TOPmed/local_ancestry/height_n12690.csv", 'Delimiter', ',');

%% change the variable to specific ancestry: 
% 'All': all people; 'Europe': European; 'America': Amerindian; 'Africa': African
anc = {'All'}

if isequal(anc,{'All'});
    pheno = pheno_all;
else
    anc_ind = find(strcmp(anc, pheno_all.Kmax)); %% change if select ancestry
    pheno = pheno_all(anc_ind,:);
end

ID = pheno(:,2);
IDs=table2array(ID);
age = pheno.AGE;
gender = pheno.GENDERNUM;
Height = pheno.HEIGHT;
Center_vec = pheno.CENTERNUM-1;
Center_vec = categorical(Center_vec);
kmax_vec = pheno.Kmax;

PC_data = readtable("/space/broce-syn01/1/data/GWAS/TOPmed/postimputaion/tmp/step4_pruned/pcs_xinwang/projections.txt");
PCID = PC_data(:,2);
PCID =table2array(PCID);
PCs = PC_data(:,3:42); %

[C,ia,ib,ic] = intersect3(FID,IDs,PCID,'stable'); % find overlap subjects with genetic and pheno data
FID_new = FID(ia); IDs_new =  IDs(ib); PCID_new = PCID(ic);
% check the overlap subjects
disp([FID_new(1:10) IDs_new(1:10) PCID_new(1:10)])

geno = single(geno_int8);
geno = geno(ia,:);
ages = age(ib);
genders = gender(ib);
phenovec = Height(ib,:);
Center = Center_vec(ib);
PCs = PCs(ic,:);
PCs =table2array(PCs);

covariates = [ages PCs];

% delete nan
phe_cov_i = find(isfinite(sum(covariates,2)+phenovec));
FID_new_sol = FID_new(phe_cov_i);
phenovec_sol = phenovec(phe_cov_i);
cova_sol = covariates(phe_cov_i,:);
genders_sol = genders(phe_cov_i);
genomat_sol = geno(phe_cov_i,:);
PCs_sol = PCs(phe_cov_i,:);
Center_sol = Center(phe_cov_i);
%%%%%%%%%%%%%%%%%%%%%%%%%% PGS only prediction %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% phenotype was regressing out of age and sex (catergorical)
% Residualize phenotype (age/sex) % change
cova = cova_sol(:,1); %% age
cova= cat(2,cova,ones(size(phenovec_sol)));
phenovec_res = phenovec_sol - cova*(pinv(cova)*phenovec_sol);
genders_sol = categorical(genders_sol); %% sex
mdl = fitlm(genders_sol,phenovec_res);
phenovec_res = mdl.Residuals.Standardized;

defvec = isfinite(phenovec_res+sum(genomat_sol,2));
if 1 % Eliminate outliers?
  ivec_cens = find(defvec&(abs(phenovec_res-mean(phenovec_res(defvec)))/std(phenovec_res(defvec))<4));
else
  ivec_cens = find(defvec);
end

FID_new_sol1 = FID_new_sol(ivec_cens);
phenovec_res = phenovec_res(ivec_cens);
genomat = genomat_sol(ivec_cens,:);
phenovec_tmp = phenovec_res;
phenovec_tmp = (phenovec_tmp-mean(phenovec_tmp))/std(phenovec_tmp);

nfracvec = [0.01 0.03 0.1 0.3 0.9]; pthreshlist = logspace(-16,-4,4);  alphalist = [0.0 1.0]
iter = 1; niter_val = 1000;
pthreshlist = [10^-16,10^-12, 10^-8, 10^-4, 10^-3, 10^-2, 10^-1];
%pthreshlist = logspace(-16,-4,4);
Type = {}; Size = {}; NumberSNP = {};Method = {}; Summary_r = []; Regression = {};  
Reg_Method = {}; iteration_validation = [];

for pthreshi = 1:length(pthreshlist)
    for alphai = 2 
        alpha = alphalist(alphai);
        for regtype = 1
            for nfraci = 5
                fname = sprintf('/space/broce-syn01/1/data/GWAS/xinwang/SOLINCA/results_revision/UKB_EUR_train_use_sGWASUKB/UKB_height_betas_iter%d_pthreshi%d_nfraci%d_regtype%d_alpha=%0.1f.mat',iter,pthreshi,nfraci,regtype,alpha);
                load(fname)
                [snp snp_ind UKB_snp_ind] = intersect(snpid_sol, snpids_tmp,"stable");
                UKB_snp = snpids_tmp(UKB_snp_ind);
                UKB_snp_A1 = A1vec_tmp(UKB_snp_ind);
                UKB_snp_A2 = A2vec_tmp(UKB_snp_ind);
                A1vec_solinca_tmp = A1vec_solinca_tmp(UKB_snp_ind);
                A2vec_solinca_tmp = A2vec_solinca_tmp(UKB_snp_ind);
                snpflipvec_tmp = snpflipvec_tmp(UKB_snp_ind);
                beta_snps = beta(UKB_snp_ind,:);
                SOL_snp = snpid_sol(snp_ind);
                SOL_snp_A1 = A1vec_sol(snp_ind);
                SOL_snp_A2 = A2vec_sol(snp_ind);
                % check the intersect SNPs and A1/A2
                disp([UKB_snp(1:10) SOL_snp(1:10)])
                disp([UKB_snp_A1(1:10) UKB_snp_A2(1:10) A1vec_solinca_tmp(1:10) A2vec_solinca_tmp(1:10) SOL_snp_A1(1:10) SOL_snp_A2(1:10)])
                genomat_se = genomat(:,snp_ind);
                flip_ind = find(snpflipvec_tmp == -1);
                for i = 1:numel(flip_ind)
                    genomat_se(:,flip_ind(i)) = 2 - genomat_se(:,flip_ind(i));
                end
                genomat_se(genomat_se<0) = 0;
                genomat_se(genomat_se>2) = 0;
                clear A; A = genomat_se; A = (A-mean(A));
                prediction = zeros(size(genomat_se,1),size(beta_snps,2));
                for n = 1:size(beta_snps,2)
                    b = beta_snps(:,n);
                    prediction(:,n) = A*b;
                end
                r = corr(phenovec_res,prediction);
                R2 = r.^2;
                [maxr2 maxr2_ind]=max(R2(1,:));
                maxr = r(maxr2_ind)
                Summary_r = [Summary_r; maxr];
                m =sprintf('reesi-pc')
                t = sprintf('regtype%d_alpha=%0.1f.csv',regtype,alpha);
                p = sprintf('pthreshi%d',pthreshi);
                s = sprintf('nfraci%d',nfraci);
                if alpha == 0
                    re = sprintf('Ridge');
                else
                    re = sprintf('Lasso');
                end
                m_re = m + "_" +re;
                Type = [Type;{t}];
                Size = [Size;{s}];
                Method =[Method;{m}];
                NumberSNP = [NumberSNP;{p}];
                Regression = [Regression;{re}];
                Reg_Method = [Reg_Method;{m_re}];
                iteration_validation = [iteration_validation;1];
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%% PGS + PCs prediction %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% phenotype was regressing out of age and sex (catergorical)

% Residualize phenotype (age/sex) % change
cova = cova_sol(:,1); %% age
cova= cat(2,cova,ones(size(phenovec_sol)));
phenovec_res = phenovec_sol - cova*(pinv(cova)*phenovec_sol);
%genders_sol2 = categorical(genders_sol); 
mdl = fitlm(genders_sol,phenovec_res); %residulize sex
phenovec_res = mdl.Residuals.Standardized;  

defvec = isfinite(phenovec_res+sum(genomat_sol,2));
if 1 
  ivec_cens = find(defvec&(abs(phenovec_res-mean(phenovec_res(defvec)))/std(phenovec_res(defvec))<4));
else
  ivec_cens = find(defvec);
end

FID_new_sol2 = FID_new_sol(ivec_cens);
phenovec_res = phenovec_res(ivec_cens);
genomat = genomat_sol(ivec_cens,:);
cova_sol = cova_sol(ivec_cens,:);
PCs_sol =cova_sol(:,2:end);
PCs_sol = array2table(PCs_sol);
Center_sol2 = Center_sol(ivec_cens);
% Center and normalize phenotype
phenovec_tmp = phenovec_res;
phenovec_tmp = (phenovec_tmp-mean(phenovec_tmp))/std(phenovec_tmp);



%%%%%%%%%%%%% 3rd degree using phi cutoff 0.0442 %%%%%%%%%%%%%%%%%%%%%%%
Relat = readtable("/space/broce-syn01/1/data/GWAS/TOPmed/postimputaion/tmp/step4_pruned/pcs_xinwang/pruned_subjects_phi_0.0442.csv", 'Delimiter', ',');
Relat_ID = unique([Relat.row;Relat.column]); % find relative ID; 
[dummy ia_rel ib_rel] =intersect(FID_new_sol2,Relat_ID); 
all_indices = 1:length(FID_new_sol2); %length(FID_new_sol2)
sol_nonrela = setdiff(all_indices', ia_rel); % sol_nonrela is the index of the subjects without relatives'

%  (31%) sample were included in the training (80%) sample, so I just need to ramdomly pick 49% subjects
train = 0.8-length(ia_rel)/length(all_indices);

train_perc = (train*length(all_indices))/length(sol_nonrela);


%%%%%%%%%%%%% Incorporating PC prediction  %%%%%%%%%%%%%%%%%%%%%%%
PC_covariates = [1,2,3,4,5,10,25,40]
for pthreshi = 1:7
for alphai = 1:2 % change lasso or ridge 
        alpha = alphalist(alphai); % alpha=0: rdige; = 1: lasso
        for regtype = 1
            for nfraci = 5
                fname = sprintf('/space/broce-syn01/1/data/GWAS/xinwang/SOLINCA/results_revision/UKB_EUR_train_use_sGWASUKB/UKB_height_betas_iter%d_pthreshi%d_nfraci%d_regtype%d_alpha=%0.1f.mat',iter,pthreshi,nfraci,regtype,alpha);
                load(fname)
                [snp snp_ind UKB_snp_ind] = intersect(snpid_sol, snpids_tmp,"stable");
                UKB_snp = snpids_tmp(UKB_snp_ind);
                UKB_snp_A1 = A1vec_tmp(UKB_snp_ind);
                UKB_snp_A2 = A2vec_tmp(UKB_snp_ind);
                A1vec_solinca_tmp = A1vec_solinca_tmp(UKB_snp_ind);
                A2vec_solinca_tmp = A2vec_solinca_tmp(UKB_snp_ind);
                snpflipvec_tmp = snpflipvec_tmp(UKB_snp_ind);
                beta_snps = beta(UKB_snp_ind,:);
                SOL_snp = snpid_sol(snp_ind);
                SOL_snp_A1 = A1vec_sol(snp_ind);
                SOL_snp_A2 = A2vec_sol(snp_ind);
                % check the intersect SNPs and A1/A2
                disp([UKB_snp(1:10) SOL_snp(1:10)])
                disp([UKB_snp_A1(1:10) UKB_snp_A2(1:10) A1vec_solinca_tmp(1:10) A2vec_solinca_tmp(1:10) SOL_snp_A1(1:10) SOL_snp_A2(1:10)])
                genomat_se = genomat(:,snp_ind);
                flip_ind = find(snpflipvec_tmp == -1);
                for i = 1:numel(flip_ind)
                    genomat_se(:,flip_ind(i)) = 2 - genomat_se(:,flip_ind(i));
                end
                genomat_se(genomat_se<0) = 0;
                genomat_se(genomat_se>2) = 0;
                clear A; A = genomat_se; A = (A-mean(A));
                prediction = zeros(size(genomat_se,1),size(beta_snps,2));
                R2 = zeros(2,size(beta_snps,2));
                for n = 1:size(beta_snps,2)
                    b = beta_snps(:,n);
                    prediction(:,n) = A*b;
                end
                r = corr(phenovec_res,prediction);
                R2 = r.^2;
                [maxr2 maxr2_ind]=max(R2(1,:));
                PHS_tmp= A * beta_snps(:,maxr2_ind);
                for iter_val = 1:niter_val
                    rand = randperm(length(sol_nonrela));
                    Sol80 = [sol_nonrela(rand(1:round(train_perc*length(sol_nonrela))));ia_rel];
                    Sol20 = sol_nonrela(rand(1+round(train_perc*length(sol_nonrela)):end));
                    fid_train_vali = intersect(FID_new_sol2(Sol80),FID_new_sol2(Sol20)) 
                    PHS = array2table(PHS_tmp);
                    for pc_i = 1:length(PC_covariates)
                        cov = PC_covariates(pc_i);
                        PCs_sol1 = PCs_sol(:,1:cov);
                        PCs_sol_80 = PCs_sol1(Sol80,:);
                        PCs_sol_20 = PCs_sol1(Sol20,:);
                        phenovec_tmp_80 = phenovec_tmp(Sol80,:); phenovec_tmp_20 = phenovec_tmp(Sol20,:);
                        PHS80 = PHS(Sol80,:);PHS20 = PHS(Sol20,:);
                        PHS_PCs_80 = cat(2,PHS80,PCs_sol_80);PHS_PCs_20 = cat(2,PHS20,PCs_sol_20);
                        PHS_PCs_80a = table2array(PHS_PCs_80);PHS_PCs_20a = table2array(PHS_PCs_20);
                        tic
                        fit1 = fitlm(PHS_PCs_80a,phenovec_tmp_80);
                        toc
                        beta= fit1.Coefficients.Estimate(2:end);
                        yy_20 = PHS_PCs_20a*beta;
                        r_20 = corr(phenovec_tmp_20,yy_20);
                        Summary_r = [Summary_r;r_20];
                        m =sprintf('adding-%dpcs',cov);
                        t = sprintf('regtype%d_alpha=%0.1f.csv',regtype,alpha);
                        p = sprintf('pthreshi%d',pthreshi);  
                        s = sprintf('nfraci%d',nfraci);
                        if alpha == 0
                            re = sprintf('Ridge');
                        else
                            re = sprintf('Lasso');
                        end
                        m_re = m + "_" +re;
                        Type = [Type;{t}];
                        Size = [Size;{s}];
                        Method =[Method;{m}];
                        NumberSNP = [NumberSNP;{p}];
                        Regression = [Regression;{re}];
                        Reg_Method = [Reg_Method;{m_re}];
                        iteration_validation = [iteration_validation;iter_val];
                    end
                end
            end
        end
    end
end

R2 = Summary_r.^2
data1 = table(Size,Type,NumberSNP,Summary_r,R2,Method,Regression,Reg_Method,iteration_validation,'VariableNames',{'Traning_Sample_Size','Type','Pthreshold','r','R2','PCs','Regression','Method','iteration_validation'})
fname_out3 = sprintf('/space/broce-syn01/1/data/GWAS/xinwang/SOLINCA/results_revision/UKB_EUR_train_use_sGWASUKB/Summary_table/SOLINCA_%s_Height_prediction_trainedUKBEUR_handlerelative_phi-0.0442_10232025_sexcatergorical.csv',anc{1});
writetable(data1, fname_out3, 'Delimiter', ',', 'QuoteStrings', false);
