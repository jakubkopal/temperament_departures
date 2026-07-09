%% ========================================================================
%  LONGITUDINAL FEMA GWAS ANALYSIS
%  ========================================================================
%  Purpose: Perform FEMA-Long genome-wide association study (GWAS) 
%           with longitudinal data. This script analyzes temperament (EAS)
%           outcomes in relation to genetic variants.
%
%  Input requirements:
%    - CSV file with phenotypic data (age, sex, genetic PCs, genotyping batch)
%    - EAS outcomes (emotional, activity, social, shyness)
%    - Binary GRM (Genetic Relatedness Matrix) for genetic kinship
%    - PLINK binary files (.bed/.bim/.fam) for SNP genotypes
%    - Family structure identifiers (IID, EID, FID)
%
%  Output: Beta estimates, standard errors, z-statistics, p-values for SNP effects
%  ========================================================================

%% Add required toolboxes and functions to MATLAB path
%  These contain FEMA functions and utility scripts for genetic analysis
addpath '.../cmig_tools_latest/cmig_tools_utils/matlab'
addpath '.../cmig_tools_latest/FEMA'

%% Define file and directory paths
%  User must specify the PLINK file prefix and genetics data directory
filePLINK   = '';  % PLINK binary file prefix (.bed/.bim/.fam)
dirGenetics = '';  % Directory containing PLINK and GRM files


%% Load phenotypic data from CSV file
%  Table T contains all demographic, covariate, and outcome variables
T = readtable(strcat(path_main, '....csv'));

%% Extract and preprocess demographic variables
%  Age: continuous variable in months; z-score standardized for numerical stability
age = T{:, "age"};  % Age in months

%  Sex: binary indicator;
sex = T{:, "sex"};  

%  Genetic PCs (Principal Components): 20 PCs from genetic data used to control for 
%  population stratification and ancestry effects.
genPCs = T{:, ["PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8",...
    "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16",...
    "PC17", "PC18", "PC19", "PC20"]};  % 20 genetic principal components

%  Genotyping batch: categorical variable accounting for batch effects in genotyping process
cols = contains(T.Properties.VariableNames, 'genotyping_batch');
genBatch = table2array(T(:, cols));

%% Extract family structure identifiers (random effects)
%  These define the hierarchical structure for mixed-effects modeling:
%  IID = Individual ID (unique per person)
%  EID = Event ID (visit time point: 18, 36, or 60 months)
%  FID = Family ID (identifies family clusters)
%  eid = numeric encoding of visit time points (1, 2, or 3)
IID = cellstr(T{:, 'iid'});  % Individual identifiers
EID = cellstr(T{:, 'eid'});  % Event time point labels
FID = cellstr(T{:, 'fid'});  % Family identifiers
eid = double(categorical(EID, {'18_months', '36_months', '60_months'}, 'Ordinal', true));  % Numeric visit codes
numVisits   = length(unique(eid));  % Number of longitudinal time points
intercept   = ones(height(eid), 1);  % Intercept term for regression model

%% Extract and standardize outcome variables (EAS Temperament scales)
%  Four EAS subscales measuring child behavior across three visits:
%  - emot: Emotional reactivity
%  - act:  Activity level
%  - soc:  Sociability 
%  - shy:  Shyness
%  All outcomes z-score standardized for consistent scaling and ease of interpretation
ymat = T{:, ["emot", "act", "soc", "shy"]};  % EAS outcome scores

%% Create smooth basis functions for age effect modeling
%  Natural spline basis functions capture non-linear age-related developmental effects
%  Knot placement: placed at median age for natural anchoring
%  Two representations created:
%  - raw_basisFunction: unprocessed basis for diagnostic purposes
%  - basisFunction: SVD-processed basis for numerical stability in regression
knots = zeros(1, 1);  % Initialize knot locations
knots(1,1) = median(age);  % Place knot at median age

%  Generate basis functions: 'nsk' = natural spline with k degrees of freedom
%  Raw basis: unprocessed natural spline basis functions
raw_basisFunction = createBasisFunctions(age, knots, 'nsk', [], 'raw');

%  SVD-orthogonalized basis: numerically stable orthogonal representation
basisFunction     = createBasisFunctions(age, knots, 'nsk', [], 'svd');

%% Construct covariate matrix for FEMA model
%  X matrix includes all fixed effects:
%  - Intercept: overall mean
%  - Basis functions: age effects (smooth across visits)
%  - Sex: binary sex indicator
%  - Genetic PCs: ancestry correction (20 PCs)
%  - Genotyping batch: technical covariate controlling for batch variation
X = [intercept, basisFunction, sex, genPCs, genBatch];

%  Ensure design matrix is full rank (required for statistical inference)
%  Non-full rank indicates redundant or collinear covariates
if rank(X) ~= size(X,2)
    error('X is not full rank - recheck covariates');
end

%% ========================================================================
%  GENETIC RELATEDNESS MATRIX (GRM) CONSTRUCTION AND VALIDATION
%  ========================================================================
%  The GRM encodes pairwise genetic similarity between all individuals,
%  accounting for family relationships. Required for family-based analysis.
%  Values range: diagonal = 1, off-diagonal = kinship coefficient (0 to 1)

%% Load GRM metadata and construct kinship matrix
%  grmInfo: structure containing sample identifiers and GRM matrix dimensions
grmInfo = load(fullfile(path_main, '....mat'));

%% Read binary GRM file
%  GRM stored as binary doubles in compressed format (lower triangle only)
%  Dimensions: [n_samples × n_samples] where n_samples = number of unique individuals
fid      = fopen(fullfile(path_main, 'GRM.dat'), 'r');
pihatmat = fread(fid, [length(grmInfo.uqObservations) length(grmInfo.uqObservations)], 'double');
fclose(fid);

%% Expand GRM to full symmetric matrix
%  GRM file contains only lower triangle; create full matrix by adding its transpose
%  This produces symmetric matrix where pihatmat(i,j) = pihatmat(j,i)
pihatmat = pihatmat + pihatmat.';

%% Set diagonal to 1
%  Diagonal elements represent self-kinship. After symmetrization, diagonal may be ~2,
%  so reset to 1 (self-kinship = 1 by definition)
tmp = size(pihatmat,1);
pihatmat(1:1+tmp:tmp*tmp) = 1;  % Linear indexing to access diagonal elements

%% Convert to single precision for memory efficiency
%  Double precision necessary for stable GRM computation; convert to single for storage/computation
pihatmat = single(pihatmat);

%% Validate and reorder GRM to match phenotype sample order
%  GRM row/column order may differ from phenotype data order.
%  This section checks alignment and reorders GRM if necessary.

%  Extract unique individual IDs in order of appearance (stable sort)
[a, b, c] = unique(IID, 'stable');  % a = unique IDs, b = first occurrence indices

%  Check 1: Do we have the same number of individuals?
%  If not, GRM reordering is needed (possible missing individuals or duplicates)
if length(a) ~= length(grmInfo.uqObservations)
    reorderGRM = true;
else
    %  Check 2: Are individual IDs in the same order and spelling?
    %  strcmpi performs case-insensitive string comparison; sum counts exact matches
    if sum(strcmpi(grmInfo.uqObservations, a)) ~= length(a)
        reorderGRM = true;
    else
        reorderGRM = false;  % IDs match exactly; no reordering needed
    end
end

%% Reorder GRM rows and columns to match phenotype data order (if needed)
%  ismember finds the index position of each phenotype IID in the GRM's IID list
%  These indices used to permute GRM rows and columns accordingly
if reorderGRM
    [~, wch] = ismember(a, grmInfo.uqObservations);  % Map phenotype IIDs to GRM indices
    pihatmat = pihatmat(wch, wch);  % Reorder GRM to match phenotype order
end

%% ========================================================================
%  FEMA MODEL CONFIGURATION
%  ========================================================================
%  Parameters defining the variance components and analysis approach

%  Random Effects: variance components to estimate
%  'F' = Family effect (shared family environment)
%  'A' = Additive genetic effect (heritability)
%  'S' = Shared environment effect (common to siblings)
%  'E' = Individual-specific environment effect
RandomEffects    = {'F', 'A', 'S', 'E'};

%  returnReusable: return pre-computed components for efficient GWAS iterations
%  (reduces computation time when testing multiple genetic variants)
returnReusable  = true;

%  contrasts: custom contrasts for comparing effects (empty = default)
contrasts       = [];

%  nbins: number of bins for stratified analysis (e.g., by family size or relationship type)
nbins           = 20;

%  niter: number of optimization iterations for parameter estimation convergence
niter           = 1;

%  CovType: structure of variance-covariance matrix for random effects
%  'unstructured' = fully flexible (each variance component estimated separately)
CovType         = 'unstructured';

%  GWAS-specific settings for SNP analysis
chunkSize       = 5000;  % Process 5000 SNPs per computational chunk
splitBy         = 'snp';  % Divide computational work by SNPs (vs. by family, etc.)
precision  = 'double';  % Numerical precision for computation

%  Parallel computing configuration
numWorkers      = ;  % Number of parallel workers (CPUs) to utilize
numThreads      = ;  % Threads per worker
GWAS_outPrefix  = 'FEMA_GWAS_Aggregate';  % Prefix for output files from GWAS

%% ========================================================================
%  FIT FEMA MODEL TO ESTIMATE BASELINE VARIANCE COMPONENTS
%  ========================================================================
%  FEMA_fit: Main function fitting mixed-effects model with family structure
%  and genetic relatedness. Estimates variance components for random effects.
%  
%  Outputs:
%  - beta_hat, beta_se: fixed effect estimates and standard errors (covariates)
%  - zmat, logpmat: z-statistics and p-values for covariate effects
%  - sig2tvec, sig2mat: estimated variance components (time-varying and constant)
%  - Hessmat, logLikvec: Hessian matrix and log-likelihood (for model diagnostics)
%  - reusableVars: pre-computed terms reused in subsequent GWAS iterations

[beta_hat,      beta_se,        zmat,        logpmat,              ...
 sig2tvec,      sig2mat,        Hessmat,     logLikvec,            ...
 beta_hat_perm, beta_se_perm,   zmat_perm,   sig2tvec_perm,        ...
 sig2mat_perm,  logLikvec_perm, binvec_save, nvec_bins,            ...
 tvec_bins,     FamilyStruct,   coeffCovar,  unstructParams,       ...
 residuals_GLS, info] =                                            ...
 FEMA_fit(X, IID, eid, FID, age, ymat, niter, contrasts, nbins, ...
          pihatmat, 'RandomEffects', RandomEffects, ...
          'CovType', CovType, 'returnReusable', returnReusable, 'precision', precision);
%  
%  Input parameters:
%  - X: covariate matrix (fixed effects)
%  - IID, EID, FID: individual, event, family identifiers
%  - age: continuous age variable
%  - ymat: outcome variable
%  - pihatmat: genetic relatedness matrix from GRM
%  - RandomEffects: variance components to estimate (Familial, Additive genetic, Shared environment, Error)

%% ========================================================================
%  EXTRACT GENETIC DATA FROM PLINK FILES
%  ========================================================================
%  Parse PLINK binary files (.bed/.bim/.fam) to extract SNP information
%  and validate that all phenotyped individuals have genetic data

%  FEMA_parse_PLINK: Extract chromosome, SNP ID, base pair position from PLINK files
%  - Chr: chromosome number for each SNP
%  - SNPID: SNP identifier (rs number or chr:pos:ref:alt)
%  - BP: base pair position on chromosome
%  - genInfo: additional genetic mapping information
%  - check: error checking flag (true = successful parse)
%  - errMsg: error message if parse fails
[~, Chr, SNPID, BP, check, errMsg, genInfo] = ...
 FEMA_parse_PLINK(fullfile(dirGenetics, filePLINK), IID, [], true);

%  Error handling: abort if PLINK file parsing failed
if ~check
    error(errMsg);
end

%% ========================================================================
%  COMPILE VARIANCE COMPONENT TERMS FOR GWAS
%  ========================================================================
%  Pre-compute and compile variance components, family relationships, and
%  covariance structures needed for SNP effect estimation. These "W terms"
%  depend only on baseline model, not on individual SNP genotypes.
%  Compiling once avoids redundant recalculations across all SNPs.

[allWsTerms, tCompile] = FEMA_compileTerms(FamilyStruct.clusterinfo, binvec_save,            ...
                                           sig2mat, RandomEffects, FamilyStruct.famtypevec,  ...
                                           GroupByFamType, CovType, precision, unstructParams.visitnum);

%% ========================================================================
%  DIVIDE SNPS INTO COMPUTATIONAL CHUNKS
%  ========================================================================
%  Split SNP list into manageable chunks for parallel processing
%  Each chunk contains SNP metadata and output filename designation

[splitInfo, tDivide] = divideSNPs(fullfile(dirGenetics, filePLINK), splitBy, ...
                                 chunkSize, GWAS_outPrefix, SNPID, Chr, BP, genInfo);

%% ========================================================================
%  CONFIGURE PARALLEL COMPUTING
%  ========================================================================
%  Initialize parallel cluster for computationally intensive GWAS analysis
%  across thousands of SNPs using multiple processors

%  Create local cluster object with specified threading configuration
local            = parcluster('local');  % Use local machine resources
local.NumThreads = ;  % Threads per MATLAB worker

%  Start parallel pool: ... workers with ... idle timeout before auto-shutdown
%  This pool will be used for parfor loop below to parallelize SNP testing
pool             = local.parpool(, 'IdleTimeout', );

%% ========================================================================
%  MAIN GWAS: TEST EACH SNP FOR ASSOCIATION WITH OUTCOMES
%  ========================================================================
%  For each SNP chunk, estimate fixed effect (regression coefficient) and its
%  standard error, adjusted for covariates, family structure, and ancestry.
%  FEMA accounts for family relationships via the kinship matrix.

%  Prepare residual phenotypes from baseline FEMA fit
%  (phenotype adjusted for covariates but not SNP effects)
ymat_res_gls = residuals_GLS;

%  Basis functions for SNP analysis: intercept + age smooth trend
%  (lower-dimensional basis for computational efficiency in GWAS)
bfSNP        = [intercept, basisFunction];

%  Loop over SNP chunks in parallel
%  Each iteration tests a chunk of ~5000 SNPs and outputs results to file
parfor parts = 1:length(splitInfo)
    %  FEMA_fit_GWAS: Estimate SNP effect for each variant in this chunk
    %  Inputs:
    %    - splitInfo{parts}: metadata for this SNP chunk (SNP names, positions, etc.)
    %    - ymat_res_gls: residualized phenotypes from baseline model
    %    - binvec_save: family cluster assignments
    %    - X: covariate matrix
    %    - allWsTerms: pre-compiled variance components (from variance component model)
    %  Outputs: Beta estimates and standard errors written to files
    FEMA_fit_GWAS(splitInfo{parts}, ymat_res_gls, binvec_save, X, allWsTerms, 'outDir', dirOutput, 'outName', splitInfo{parts}.outName, 'bfSNP', bfSNP, 'doCoeffCovar', true, 'precision', precision);
end

%% ========================================================================
%  CLEANUP: SHUT DOWN PARALLEL CLUSTER
%  ========================================================================
%  Release computational resources after GWAS completes
delete(pool);  % Terminate parallel pool and free memory
