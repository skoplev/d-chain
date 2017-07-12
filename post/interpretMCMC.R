# Functions for parsing and interpretting MCMC

# Parses output files from mcmcGlobal fitting run.
# folder is the directory that the .csv MCMC sample files are stored.
parseMCMC = function(folder) {
	require(data.table)
	require(compiler)
	enableJIT(3)

	require(reshape2)
	require(stringr)
	require(abind)

	d = list()

	# Load MCMC samples for the strain in specified folder
	# --------------------------------------------------------

	# Load drug indices
	d$drugs = fread(file.path(folder, "drugs.csv"),
		nrows=1,
		header=FALSE)  # drug names matching indicies
	d$drugs = as.vector(t(d$drugs))

	d$beta_residual = fread(file.path(folder, "beta_residual.csv"),
		header=FALSE)  # residual norm factor

	d$lambda = fread(file.path(folder, "lambda.csv"),
		header=FALSE)

	# Theta, single drugs
	theta_tmp = fread(file.path(folder, "theta.csv"),
		header=FALSE)  # baseline parameters
	theta_tmp = as.data.frame(theta_tmp)
	colnames(theta_tmp) = d$drugs
	theta_tmp$sample = 1:nrow(theta_tmp)  # add MCMC sample # for melt IDs

	theta_tmp_flat = melt(theta_tmp, id.vars="sample")

	# Expand theta (parameters) matrix.
	# Returns character matrix.
	theta_mat = matrix(
		unlist(
			strsplit(theta_tmp_flat$value, ";")
		),
		ncol=3,
		byrow=TRUE
	)
	colnames(theta_mat) = c("K", "h", "alpha")
	theta_mat = apply(theta_mat, 2, as.numeric)  # char->numeric

	# Combine flat 
	theta_tmp_flat = cbind(theta_tmp_flat, theta_mat)

	# Drop ; separated string
	# theta_tmp_flat = theta_tmp_flat[, colnames(theta_tmp_flat) != "value"]

	d$theta = abind(
		acast(theta_tmp_flat, sample~variable, value.var="K"),
		acast(theta_tmp_flat, sample~variable, value.var="h"),
		acast(theta_tmp_flat, sample~variable, value.var="alpha"),
		along=3
	)
	dimnames(d$theta)[[3]] = c("K", "h", "alpha")

	# d$theta = array(dim=c(nrow(theta_tmp), ncol(theta_tmp), 3));
	# for (i in 1:nrow(theta_tmp)) {
	# 	for (j in 1:ncol(theta_tmp)) {
	# 			row = strsplit(toString(theta_tmp[i, j]), ";")
	# 			d$theta[i, j,] = as.numeric(row[[1]])
	# 	}
	# }

	message("loading lambdas...")
	# Combined experiments
	lambda_AB_tmp = fread(file.path(folder, "lambda_AB.csv"),
		header=FALSE)
	lambda_AB_tmp = as.matrix(lambda_AB_tmp)

	message("loading lambdas...")
	d$lambda_AB = array(dim=c(nrow(lambda_AB_tmp), length(d$drugs), length(d$drugs)))
	for (i in 1:nrow(lambda_AB_tmp)) {
		d$lambda_AB[i,,] = matrix(lambda_AB_tmp[i,], ncol=length(d$drugs), nrow=length(d$drugs), byrow=TRUE)
	}

	# Data structure transformations
	# ----------------------------------------------------
	# Convert theta_AB_tmp to array. Takes ~5min or longer
	message("loading thetas...")

	theta_AB_tmp = fread(file.path(folder, "theta_AB.csv"),
		header=FALSE)

	d$theta_AB = array(dim=(c(nrow(theta_AB_tmp), ncol(theta_AB_tmp), 3)))
	for (i in 1:nrow(theta_AB_tmp)) {
		row = sapply(theta_AB_tmp[i,], toString)  # list of strings
		row = sapply(row, strsplit, ";")

		for (j in 1:length(row)) {
			d$theta_AB[i, j,] = as.numeric(row[[j]])
		}
	}

	# Convert to array for [a][b] access
	d$theta_AB_arr = array(dim=c(nrow(theta_AB_tmp), length(d$drugs), length(d$drugs), 3))
	for (i in 1:nrow(theta_AB_tmp)) {
		for (j in 1:3) {
			d$theta_AB_arr[i,,,j] = matrix(d$theta_AB[i,,j], ncol=length(d$drugs), nrow=length(d$drugs), byrow=TRUE)
		}
	}

	return(d)
}
