
% add paths
addpath(genpath('/home/dale/matlab'))
addpath(genpath('/space/broce-syn01/1/data/GWAS/PHS/Test_PHS/scripts'))

% Get genotype
bfile_format = '/space/ceph/1/UKB/DATA/genotypes/HRC_imputed/full_sample/QCed/ukb27412_imp_chr%s_v3';

chrlist = cellfun(@(x)num2str(x),num2cell([1:22]),'UniformOutput',false); 
snpidvec_geno = {}; A1vec_geno = {}; A2vec_geno = {}; chrnumvec_geno = []; posvec_geno = [];
for chri = 1:length(chrlist)
  chr = chrlist{chri};
  bfile = sprintf(bfile_format,chr);
  fileID = fopen(sprintf('%s.bim', bfile));
  bim_file = textscan(fileID,'%s %s %s %s %s %s');
  fclose(fileID);
  snpidvec = bim_file{2};
  snps=length(bim_file{1});
  fileID = fopen(sprintf('%s.fam', bfile));
  fam_file = textscan(fileID,'%s %s %s %s %s %s');
  fclose(fileID);
  nsubj=length(fam_file{1});
  PIDs_geno = fam_file{1};
  FIDs_geno = fam_file{2}; 
  fprintf('%i snps and %i subjects detected in bfile\n', snps, nsubj);
  chrnumvec = str2double(bim_file{1});
  posvec = str2double(bim_file{4});
  A1vec = bim_file{5};
  A2vec = bim_file{6};
  snpidvec_geno = cat(1,snpidvec_geno,snpidvec);
  A1vec_geno = cat(1,A1vec_geno,A1vec);
  A2vec_geno = cat(1,A2vec_geno,A2vec);
  chrnumvec_geno = cat(1,chrnumvec_geno,chrnumvec);
  posvec_geno = cat(1,posvec_geno,posvec);
end

% Read in saturated GWAS SNPs
tbl_gwas = readtable('/space/broce-syn01/1/data/GWAS/xinwang/SOLINCA/GIANT_HEIGHT_YENGO_2022_GWAS_SUMMARY_STATS_EUR','FileType','delimitedtext','delimiter','\t');
snpidvec_gwas = tbl_gwas.RSID;
betavec_gwas = tbl_gwas.BETA;
pvec_gwas = tbl_gwas.P;

% Read in SOL-Inca SNPs
bfile = '/space/broce-syn01/1/data/GWAS/TOPmed/postimputaion/MEGA_HCHS_SOL.MERGED.QC.SNPSonly';
fileID = fopen(sprintf('%s.bim', bfile));
bim_file = textscan(fileID,'%s %s %s %s %s %s');
fclose(fileID);
snpidvec_solinca = bim_file{2};
chrnumvec_solinca = str2double(bim_file{1});
posvec_solinca = str2double(bim_file{4});
A1vec_solinca = bim_file{5};
A2vec_solinca = bim_file{6};

% Should make sure that genomic info for solinca is consistent with geno
[snpidvec_com I_geno I_solinca] = intersect(snpidvec_geno,snpidvec_solinca,'stable');
fprintf(1,'Percentage of SNPs on conflicting chromosome: %0.2f%%\n',100*mean(chrnumvec_geno(I_geno)~=chrnumvec_solinca(I_solinca)));
fprintf(1,'Percentage of SNPs with conflicting position: %0.2f%%\n',100*mean(posvec_geno(I_geno)~=posvec_solinca(I_solinca)));
%snpflipvec_conflicting = GenStats_strandflip(ones(size(snpidvec_com)),A1vec_geno(I_geno),A2vec_geno(I_geno),A1vec_solinca(I_solinca),A2vec_solinca(I_solinca),false);
%fprintf(1,'Percentage of SNPs with conflicting A1/A2 alleles: %0.2f%%\n',100*mean(~isfinite(snpflipvec_conflicting)));


% Intersect SNPs

[snpvec_intersect I_snp_geno I_snp_gwas I_snp_solinca] = intersect3(snpidvec_geno,snpidvec_gwas,snpidvec_solinca);

% Check that intersectn worked
disp([snpvec_intersect([1:10 end]) snpidvec_geno(I_snp_geno([1:10 end])) snpidvec_gwas(I_snp_gwas([1:10 end])) snpidvec_solinca(I_snp_solinca([1:10 end]))])
disp([A1vec_geno(I_snp_geno([1:10 end])) A2vec_geno(I_snp_geno([1:10 end])) A1vec_solinca(I_snp_solinca([1:10 end])) A2vec_solinca(I_snp_solinca([1:10 end]))])

% Identify and censor ambiguous or inconsistent SNPs
snpflipvec_conflicting = GenStats_strandflip(ones(size(snpvec_intersect)),A1vec_geno(I_snp_geno),A2vec_geno(I_snp_geno),A1vec_solinca(I_snp_solinca),A2vec_solinca(I_snp_solinca),false);
ivec_conflicting = find(~isfinite(snpflipvec_conflicting));
for i = 1:length(ivec_conflicting)
  snpid = snpvec_intersect{ivec_conflicting(i)};
  fprintf(1,'%d %d %s %c %c %c %c\n',i,ivec_conflicting(i),snpid,A1vec_geno{I_snp_geno(ivec_conflicting(i))},A2vec_geno{I_snp_geno(ivec_conflicting(i))},A1vec_solinca{I_snp_solinca(ivec_conflicting(i))},A2vec_solinca{I_snp_solinca(ivec_conflicting(i))});
end

snpflipvec_ambiguous = GenStats_strandflip(ones(size(snpvec_intersect)),A1vec_geno(I_snp_geno),A2vec_geno(I_snp_geno),A1vec_solinca(I_snp_solinca),A2vec_solinca(I_snp_solinca),true);
ivec_ambiguous = find(~isfinite(snpflipvec_ambiguous));
%disp([A1vec_geno(I_snp_geno(ivec_ambiguous)) A2vec_geno(I_snp_geno(ivec_ambiguous)) A1vec_solinca(I_snp_solinca(ivec_ambiguous)) A2vec_solinca(I_snp_solinca(ivec_ambiguous))]);

jvec_cens = find(isfinite(snpflipvec_ambiguous));


ivec_sub = sort(jvec_cens(find(pvec_gwas(I_snp_gwas(jvec_cens))<1e-1)),'ascend'); 
snpidvec_sub = snpvec_intersect(ivec_sub);

% Read genotypes into memory for thresholded, intersected, censored SNPs
snpidvec_cat = {}; A1vec_cat = {}; A2vec_cat = {}; chrnumvec_cat = {}; posvec_cat = {};
genomat_cat = [];
for chri = 1:22 % Ignore chromosome 23, which has inconsistent subject ordering
  fprintf(1,'chri=%d/%d (%s)\n',chri,length(chrlist),datestr(now));
  chr = chrlist{chri};
  bfile = sprintf(bfile_format,chr);
  fileID = fopen(sprintf('%s.bim', bfile));
  bim_file = textscan(fileID,'%s %s %s %s %s %s');
  fclose(fileID);
  chrnumvec = str2double(bim_file{1});
  chrnumlist = unique(chrnumvec);
  snpidvec = bim_file{2};
  posvec = str2double(bim_file{4});
  A1vec = bim_file{5};
  A2vec = bim_file{6};
  snps=length(bim_file{1});
  fileID = fopen(sprintf('%s.fam', bfile));
  fam_file = textscan(fileID,'%s %s %s %s %s %s');
  fclose(fileID);
  nsubj=length(fam_file{1});
  PIDs_geno = fam_file{1};
  FIDs_geno = fam_file{2}; 
  fprintf('%i snps and %i subjects detected in bfile\n', snps, nsubj);
  [dummy IA_tmp IB_tmp] = intersect(snpidvec,snpidvec_sub,'stable');
  tic
  geno_int8 = PlinkRead_binary2_amd(nsubj, IA_tmp, bfile, 1:nsubj); 
  toc
  geno = single(geno_int8);
%  geno(geno<0) = NaN;
  geno(geno<0) = 0; % "Impute" missing as 0 (assumes coded as minor allele copies)
  genomat_cat = cat(2,genomat_cat,geno);
  snpidvec_cat = cat(1,snpidvec_cat,snpidvec(IA_tmp)); 
  A1vec_cat = cat(1,A1vec_cat,A1vec(IA_tmp));
  A2vec_cat = cat(1,A2vec_cat,A2vec(IA_tmp));
  chrnumvec_cat = cat(1,chrnumvec_cat,chrnumvec(IA_tmp));
  posvec_cat = cat(1,posvec_cat,posvec(IA_tmp));
  PrintMemoryUsage
end
clear geno
clear geno_int8

% Read height phenotype 
tbl_pheno = readtable('/space/ceph/1/UKB/DATA/users/amdale/height.csv');
eidvec_pheno = cellfun(@(x)num2str(x),num2cell(tbl_pheno.eid),'UniformOutput',false); 
phenovec_pheno = nanmean(table2array(tbl_pheno(:,2:end)),2); 

% Read covariates - Hao to coordinate with Rob
tbl_demogShort = readtable('/space/gwas-syn1/1/data/GWAS/UKBioBank/covariates/UKB500k_demogShort_230519.txt'); 
tbl_popComps = readtable('/space/gwas-syn1/1/data/GWAS/UKBioBank/covariates/UKB500k_popComps_230519.txt'); 
tbl_covars = cat(2,tbl_demogShort,tbl_popComps(:,2:end));
tbl_covars.ID = cellfun(@(x)num2str(x),num2cell(tbl_demogShort.ID),'UniformOutput',false); 
eidvec_covars = tbl_covars.ID;
agevec_covars = tbl_covars.Age;
sexvec_covars = NaN(size(eidvec_covars)); sexvec_covars(strcmp('Male',tbl_covars.Sex)) = 1; sexvec_covars(strcmp('Female',tbl_covars.Sex)) = 0;
PCmat_covars = table2array(tbl_popComps(:,2:end));
covmat_covars = cat(2,agevec_covars,PCmat_covars);

% read european 
Cauca_data = readtable('/space/broce-syn01/1/data/GWAS/xinwang/UKB/22006.csv');
Cauca_eid = cellfun(@(x)num2str(x),num2cell(Cauca_data.eid),'UniformOutput',false); 
Cauca = Cauca_eid(find(Cauca_data.x22006_0_0_ukb52562 ==1)); 

[dummy I_subj_covars I_subj_pheno I_subj_geno I_subj_cauca] = intersectn(eidvec_covars,eidvec_pheno,PIDs_geno,Cauca); 

ivec = find(isfinite(sum(covmat_covars(I_subj_covars,:),2)+phenovec_pheno(I_subj_pheno)));  
%disp([eidvec_covars(I_subj_covars(ivec(1:10))) eidvec_pheno(I_subj_pheno(ivec(1:10))) PIDs_geno(I_subj_geno(ivec(1:10)))]) 
phenovec = double(phenovec_pheno(I_subj_pheno(ivec),:));
eidlist = eidvec_covars(I_subj_covars(ivec));
covmat = cat(2,covmat_covars(I_subj_covars(ivec),:),ones(size(phenovec))); 
agevec = agevec_covars(I_subj_covars(ivec));
sexvec = categorical(sexvec_covars(I_subj_covars(ivec)));
genomat = genomat_cat(I_subj_geno(ivec),:);
clear genomat_cat

% Residualize phenotype
phenovec_res = phenovec - covmat*(pinv(covmat)*phenovec); % Residualize by age and PCs
mdl = fitlm(sexvec,phenovec_res); % Residualize by sex as catergorical
phenovec_res = mdl.Residuals.Standardized;
defvec = isfinite(phenovec_res+sum(genomat,2));
if 1 
  ivec_cens = find(defvec&(abs(phenovec_res-mean(phenovec_res(defvec)))/std(phenovec_res(defvec))<4));
else
  ivec_cens = find(defvec);
end

phenovec = phenovec(ivec_cens);
phenovec_res = phenovec_res(ivec_cens);
genomat = genomat(ivec_cens,:);
agevec = agevec(ivec_cens);
sexvec = sexvec(ivec_cens);
covmat = covmat(ivec_cens,:);
eidlist = eidlist(ivec_cens); 
% Center and normalize phewnotype
phenovec_tmp = phenovec_res;
phenovec_tmp = (phenovec_tmp-mean(phenovec_tmp))/std(phenovec_tmp);

disp(datestr(now))

jf=java.text.DecimalFormat;

%keyboard

[snpvec_intersect_new I_snp_geno_new I_snp_gwas_new I_snp_solinca_new] = intersect3(snpidvec_cat,snpidvec_gwas, snpidvec_solinca,'stable');

snpidvec_cat_new = snpidvec_cat(I_snp_geno_new);
snpidvec_sub = snpidvec_gwas(I_snp_gwas_new);
snpidvec_solinca =snpidvec_solinca(I_snp_solinca_new);
% check the snps id are sorting again
disp([snpvec_intersect_new(1:10) snpidvec_cat_new(1:10) snpidvec_sub(1:10) snpidvec_solinca(1:10)])

genomat_new = genomat(:,I_snp_geno_new);
A1vec_geno = A1vec_cat(I_snp_geno_new);
A2vec_geno = A2vec_cat(I_snp_geno_new);
clear genomat
% check the snps id and A1 A2 are matching
disp([snpidvec_cat_new(1:10) A1vec_geno(1:10) A2vec_geno(1:10)])

pvec_gwas_sub = pvec_gwas(I_snp_gwas_new);

A1vec_solinca =A1vec_solinca(I_snp_solinca_new);
A2vec_solinca =A2vec_solinca(I_snp_solinca_new);
disp([snpidvec_cat_new(1:10) A1vec_solinca(1:10) A2vec_solinca(1:10)])

snpflipvec_ambiguous = GenStats_strandflip(ones(size(A1vec_solinca)),A1vec_geno,A2vec_geno,A1vec_solinca ,A2vec_solinca,true);

% Perform random cross-validation using glmnet
nfracvec = [0.01 0.03 0.1 0.3 0.9]; 
pthreshlist = [10^-16,10^-12, 10^-8, 10^-4, 10^-3, 10^-2, 10^-1];
nvec = round(nfracvec*length(phenovec_tmp));

PrintMemoryUsage
niter = 1;
nlambda = 10; alphalist = [0.0 1.0];
lammat = NaN(niter,length(pthreshlist),length(nfracvec),nlambda,2); betamats = cell(niter,length(pthreshlist),length(nfracvec),2);
rmat_disc = NaN(niter,length(pthreshlist),length(nfracvec),nlambda,2); rmat_repl = NaN(niter,length(pthreshlist),length(nfracvec),nlambda,2);
for iter = 1:niter
  for alphai = 1:2
    alpha = alphalist(alphai);
    for pthreshi = 1:length(pthreshlist)
    ivec_pvec = find(pvec_gwas_sub<=pthreshlist(pthreshi)); 
    clear A; PrintMemoryUsage; A = genomat_new(:,ivec_pvec); A = (A-mean(A)); 
      for nfraci = length(nfracvec)
       %clear genomat_new
        nfrac = nfracvec(nfraci);
        permvec = randperm(length(phenovec_tmp));
        ivec_disc = permvec(1:round(nfrac*length(phenovec_tmp)));
        ivec_repl = permvec(1+round(nfrac*length(phenovec_tmp)):end);
        eidlist_repl = eidlist(ivec_repl);
        A_disc = double(A(ivec_disc,:)); phenovec_disc = double(phenovec_tmp(ivec_disc));
        for regtype = 1
          fprintf(1,'alphai=%d/%d iter=%d/%d pthreshi=%d/%d nfraci=%d/%d regtype=%d %s\n',alphai,length(alphalist),iter,niter,pthreshi,length(pthreshlist),nfraci,length(nfracvec),regtype,datestr(now));
            clear options
          options = glmnetSet;
          options.nlambda = nlambda;
          options.alpha = alpha;
          PrintMemoryUsage
          tic
          fit1 = glmnet(A_disc,phenovec_disc,'gaussian',options);
          toc
          lambda = rowvec(fit1.lambda)
          lammat(iter,pthreshi,nfraci,:,regtype) = lambda;
          beta = fit1.beta;
          yy_hatmat = A*beta;
          rvec_disc = corr(phenovec_tmp(ivec_disc),yy_hatmat(ivec_disc,:))
          rvec_repl = corr(phenovec_tmp(ivec_repl),yy_hatmat(ivec_repl,:))
          rmat_disc(iter,pthreshi,nfraci,:,regtype) = rvec_disc;
          rmat_repl(iter,pthreshi,nfraci,:,regtype) = rvec_repl;
          betamats{iter,pthreshi,nfraci,regtype} = beta;
          fname_out = sprintf('/space/broce-syn01/1/data/GWAS/xinwang/SOLINCA/results_revision/UKB_EUR_train_use_sGWASUKB/UKB_height_betas_iter%d_pthreshi%d_nfraci%d_regtype%d_alpha=%0.1f.mat',iter,pthreshi,nfraci,regtype,alpha);
          snpids_tmp = snpidvec_sub(ivec_pvec);
            A1vec_tmp = A1vec_geno(ivec_pvec);
            A2vec_tmp = A2vec_geno(ivec_pvec);
          A1vec_solinca_tmp = A1vec_solinca(ivec_pvec);
          A2vec_solinca_tmp = A2vec_solinca(ivec_pvec);
          snpflipvec_tmp = snpflipvec_ambiguous(ivec_pvec);
          save(fname_out,'beta','lambda','rvec_disc','rvec_repl','snpids_tmp','A1vec_tmp','A2vec_tmp','A1vec_solinca_tmp','A2vec_solinca_tmp','snpflipvec_tmp','eidlist_repl','rvec_repl','rvec_disc');
        end
        clear A_disc;
      end
    end
  end
end


keyboard
