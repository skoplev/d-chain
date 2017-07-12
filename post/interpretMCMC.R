# Functions for parsing and interpretting MCMC

require(data.table)
require(compiler)
enableJIT(3)

require(reshape2)
require(stringr)
require(abind)

require(dplyr)

# Decodes vector of ; (default) separated strings. 3 parameters by default.
# Returns numeric matrix.
decodeParamVector = function(param_vec, par_names=c("K", "h", "alpha"), sep=";") {
	mat = matrix(
		unlist(
			strsplit(param_vec, sep)
		),
		ncol=3,
		byrow=TRUE
	)
	colnames(mat) = par_names
	mat = apply(mat, 2, as.numeric)  # char->numeric
	return(mat)
}

# Casts flat format into 3d array
sampleArrayCast = function(flat) {
	theta = abind(
		acast(flat, sample~variable, value.var="K"),
		acast(flat, sample~variable, value.var="h"),
		acast(flat, sample~variable, value.var="alpha"),
		along=3
	)
	dimnames(theta)[[3]] = c("K", "h", "alpha")
	return(theta)
}

# Dose-response model
response <- function(conc, par) {
	K <- par[1]
	h <- par[2]
	alpha <- par[3]
	out <- ((1 - alpha) / (1 + (K * conc)^h)) + alpha
	return(out)
}

# Parses output files from mcmcGlobal fitting run.
# folder is the directory that the .csv MCMC sample files are stored.
parseMCMC = function(folder) {
	d = list()

	# Load MCMC samples for the strain in specified folder
	# --------------------------------------------------------

	message("loading parameter samples...")
	message("\tsingle-drug parameters...")
	# Load drug indices
	d$drugs = fread(file.path(folder, "drugs.csv"),
		nrows=1,
		header=FALSE)  # drug names matching indicies
	d$drugs = as.vector(t(d$drugs))

	d$beta_residual = fread(file.path(folder, "beta_residual.csv"),
		header=FALSE)  # residual norm factor
	d$beta_residual = as.data.frame(d$beta_residual)

	d$lambda = fread(file.path(folder, "lambda.csv"),
		header=FALSE)
	d$lambda = as.data.frame(d$lambda)

	# Theta, single drugs
	theta_tmp = fread(file.path(folder, "theta.csv"),
		header=FALSE)  # baseline parameters
	theta_tmp = as.data.frame(theta_tmp)
	colnames(theta_tmp) = d$drugs
	theta_tmp$sample = 1:nrow(theta_tmp)  # add MCMC sample # for melt IDs

	theta_tmp_flat = melt(theta_tmp, id.vars="sample")

	# Expand theta (parameters) matrix.
	# Returns character matrix.
	theta_mat = decodeParamVector(theta_tmp_flat$value)

	# Combine flat 
	theta_tmp_flat = do.call(cbind, list(theta_tmp_flat, theta_mat))

	d$theta = sampleArrayCast(theta_tmp_flat)

	message("\tcombinatorial lambdas...")
	# Combined experiments
	lambda_AB_tmp = fread(file.path(folder, "lambda_AB.csv"),
		header=FALSE)
	lambda_AB_tmp = as.matrix(lambda_AB_tmp)

	d$lambda_AB = array(dim=c(nrow(lambda_AB_tmp), length(d$drugs), length(d$drugs)))
	# For loop is quick enough.
	for (i in 1:nrow(lambda_AB_tmp)) {
		d$lambda_AB[i,,] = matrix(lambda_AB_tmp[i,], ncol=length(d$drugs), nrow=length(d$drugs),
			byrow=TRUE)
	}

	dimnames(d$lambda_AB)[[2]] = d$drugs
	dimnames(d$lambda_AB)[[3]] = d$drugs
	names(dimnames(d$lambda_AB)) = c("Sample", "First", "Second")

	# Data structure transformations
	# ----------------------------------------------------
	# Convert theta_AB_tmp to array. Takes ~5min or longer
	message("\tcombinatorial thetas...")

	# Load rows of vectors of parameter sets.
	# Each row vector encodes a complete set of combinatorial parameter sets.
	theta_AB_tmp = fread(file.path(folder, "theta_AB.csv"),
		header=FALSE)

	# Parse parameters
	drug_comb = expand.grid(d$drugs, d$drugs)[, 2:1]
	colnames(drug_comb) = c("first", "second")
	colnames(theta_AB_tmp) = apply(drug_comb, 1, paste, collapse="_")  # by-row encoding of matrix
	theta_AB_tmp$sample = 1:nrow(theta_AB_tmp)

	# FLatten ; separated values by condition
	theta_AB_flat = melt(theta_AB_tmp, id.vars="sample")

	# Evaluate ; separated parameters
	theta_AB_mat = decodeParamVector(theta_AB_flat$value)

	# Separate drug combination IDs, first and second treatment.
	drug_ids = matrix(
		unlist(strsplit(as.character(theta_AB_flat$variable), "_")),
		ncol=2,
		byrow=TRUE
	)
	colnames(drug_ids) = c("first", "second")

	theta_AB_flat = bind_cols(theta_AB_flat, as.data.frame(drug_ids), as.data.frame(theta_AB_mat))  # from dplyr, fast

	# Combine samples into arrays
	d$theta_AB = abind(
		acast(theta_AB_flat, sample~first~second, value.var="K"),
		acast(theta_AB_flat, sample~first~second, value.var="h"),
		acast(theta_AB_flat, sample~first~second, value.var="alpha"),
		along=4
	)
	dimnames(d$theta_AB)[[4]] = c("K", "h", "alpha")

	return(d)
}

summaryStatisticsMCMC = function(mcmc) {
	message("Calculating MCMC summary statistics...")
	d = list()  # output data
	d$beta_residual_mean = apply(mcmc$beta_residual, 2, mean)  # pretreatment specific renormalization
	d$lambda_mean = apply(mcmc$lambda, 2, mean)  # baseline 

	# baseline parameters
	d$theta_mean = array(NA, dim=c(length(mcmc$drugs), 3))  # conditional (on lambda) parameters
	for (i in 1:length(mcmc$drugs)) {
		id = which(as.logical(mcmc$lambda[,i]))  # iterations where parameters are used
		if (length(id) > 1) {
			d$theta_mean[i,] = apply(mcmc$theta[id,i,], 2, mean)
		} else if (sum(id) == 1) {
			d$theta_mean[i,] = mcmc$theta[id,i,]
		} else {
			d$theta_mean[i,] <- c(0, 1.5, 0)  # default
		}
	}

	# Calculate average parameters, lambda_AB
	d$lambda_AB_mean = matrix(0.0, ncol=length(mcmc$drugs), nrow=length(mcmc$drugs))
	for (i in 1:dim(mcmc$lambda_AB)[1]) {
		d$lambda_AB_mean = d$lambda_AB_mean + mcmc$lambda_AB[i,,]
	}
	d$lambda_AB_mean = d$lambda_AB_mean / dim(mcmc$lambda_AB)[1]

	# theta_AB
	init_par = c(0.01, 1.5, 0.0)
	d$theta_AB_mean = array(NA, dim=c(length(mcmc$drugs), length(mcmc$drugs), 3))
	for (a in 1:length(mcmc$drugs)) {
		for (b in 1:length(mcmc$drugs)) {
			id = mcmc$lambda_AB[,a,b] == TRUE
			if (sum(id) > 1) {
				d$theta_AB_mean[a, b,] = apply(mcmc$theta_AB[id, a, b,], 2, mean)
			} else if (sum(id) == 1) {
				d$theta_AB_mean[a, b,] = mcmc$theta_AB[id, a, b,]
			} else {
				d$theta_AB_mean[a, b,] = init_par
			}
		}
	}

	# Synergy index calculation based on MCMC samples
	# -------------------------------------------------
	nsample = nrow(mcmc$lambda)
	synergy_index = array(0.0, dim=c(nsample, length(mcmc$drugs), ncol=length(mcmc$drugs)))

	s = seq(0.01, 10, length.out=10)  # evaluation points for integral
	for (a in 1:length(mcmc$drugs)) {
		if (a %% 10 == 0) {
			cat("Iteration (drug) ", a, " out of ", length(mcmc$drugs), "\n")
		}
		for (b in 1:length(mcmc$drugs)) {
			for (i in 1:nsample) {
				baseline = (1 - mcmc$lambda[i, b]) + mcmc$lambda[i, b] * (response(s, mcmc$theta[i, b,]))
				resp = (1 - mcmc$lambda_AB[i, a, b]) * baseline + mcmc$lambda_AB[i, a, b] * response(s, mcmc$theta_AB[i, a, b,])

				synergy_index[i, a, b] = mean(baseline - resp)
			}
		}
	}


	d$synergy_index_mean = matrix(NA, length(mcmc$drugs), length(mcmc$drugs))
	d$synergy_index_sd = matrix(NA, length(mcmc$drugs), length(mcmc$drugs))
	for (a in 1:length(mcmc$drugs)) {
		for (b in 1:length(mcmc$drugs)) {
			d$synergy_index_mean[a, b] = mean(synergy_index[,a,b])
			d$synergy_index_sd[a, b] = sd(synergy_index[,a,b])
		}
	}
	message("done.")
	return(d)
}
