// mcmcGlobal.cpp
//
// Global MCMC analysis of systematic cell viability screen. Parametric curve fits using selector
// variables to fall back on baseline models. Experiment types AB, A, A0 have their own, albeit, related
// parametrizations. x_AB = beta * f_A * f_B, assuming Bliss independence in the two time intervals. 
// x_A = f_A, x_A0 = beta * f_A.
// T-statistic error model with Gamma prior on the log relative
// viability data.
// 
// Reads csv file in the format specified by data field codes. Handles experiments at different concentration
// by constructing experimental id's.
// 
// Compilation:  g++ -O mcmcGlobal.cpp -o mcmcGlobal -lboost_system -lboost_program_options
// Uses C++11 and <boost>
//
// Old data specification used for analysis
// string data_file_name = "all_data_linReg_validation.csv";  // was all_data_q75_norm.csv. Includes 2nd round validation
// string data_file_name = "all_data_linregnorm.csv";  // primary screen
// int iterations = 500000;
// int burn = 100000;
// int subsample = 200;  // store every nth parameter set.
// int init_phase = 20000;  // relies on init_lambda to be true.
//
// Authors:  Simon Koplev (skoplev@gmail.com)
// 

#include <iostream>
#include <random>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <cstring>
#include <set>
#include <cmath>
#include <array>
#include <utility>  // for pair
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

// DEBUG PRINT
// ----------------------------
#ifdef NDEBUG
#define DPRINT(value) ((void)0)  // removed by interpreter
#else
#define DPRINT(value)\
	cout << "DEBUG " << __FILE__ \
		<< ", " << __FUNCTION__ << "()"\
		<< ", line " << __LINE__ \
		<< ": " << #value << " = " << value << endl;
#endif

// Data fields code
#define FIELD_EXPERIMENT 0
#define FIELD_STRAIN 1
#define FIELD_RUN 2
#define FIELD_PLATE 3
#define FIELD_PRETREATMENT 4
#define FIELD_COMPOUND 5
#define FIELD_CONCENTRATION 6
#define FIELD_RELCELL 7  // relative cell count

#define NPAR 3  // number of response curve paramters

using namespace std;
namespace po = boost::program_options;
namespace fs = boost::filesystem;

// split() 
// ----------------------------
// Splits a string by delimiter character. Returns vector of string elements.
void split(const std::string &s, char delim,  // in
	std::vector<std::string> &elems)  // out
{
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
}
std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

// readFile()
// ----------------------------------------------
// Reads content of file and returns a vector of lines encoded as strings. 
// Input is path of file to read.
std::vector<std::string> readFile(const std::string &s) 
{
	ifstream in_file (s);
	std::vector<std::string> data;
	if (in_file.is_open()) {
		string line;
		while(getline(in_file, line)) {
			data.push_back(line);
		}
		in_file.close();
	} else {
		cout << "readFile() ERROR: unable to open file: " << s << endl;
		exit(EXIT_FAILURE);
	}
	in_file.close();
	return(data);
}

// genExpID()
// ------------------------------------------------------
// Generates unique experiment ID, a string. Used for identifying triplicate experiments.
string genExpID(const vector<string> &row) 
{
	string id =
		row[FIELD_EXPERIMENT] + 
		row[FIELD_STRAIN] + 
		row[FIELD_PRETREATMENT] + 
		row[FIELD_COMPOUND] + 
		row[FIELD_CONCENTRATION];
	return(id);
}

// CLASS: Observation
// ------------------------
// A single observation (data point) from cell viability experiment.
class Observation {
public:
	string experiment;
	string strain;
	int run;  // run id
	int plate;  // plate id
	string pretreatment;
	string compound;
	double concentration;  // drug concentration
	double log_rel_count;
};

// parseRow()
// -------------------------------------
// Parses data row into Observation object. vector<string> -> Observation object.
Observation parseRow(const vector<string> &d) {
	Observation obs;
	obs.experiment = d[FIELD_EXPERIMENT];
	obs.strain = d[FIELD_STRAIN];
	obs.run = stoi(d[FIELD_RUN]) - 1;  // to zero indexed
	obs.plate = stoi(d[FIELD_PLATE]) - 1;
	obs.pretreatment = d[FIELD_PRETREATMENT];
	obs.compound = d[FIELD_COMPOUND];
	obs.concentration = stod(d[FIELD_CONCENTRATION]);
	obs.log_rel_count = log(stod(d[FIELD_RELCELL]));
	return obs;
}

typedef vector<Observation> Replicates;  // Collection of observation that are replicates, same experimental condition.
typedef vector<Replicates> ExperimentSet;  // a collection of experiments

// CLASS: SufficientStat
// ------------------------------------------
class SufficientStat {
public:
	int n;  // number of observations
	double x_mean;  // log relative count mean
	double x_var;  // log relative count variance. 0 if only one observation.
	double conc;  // concentration
};

// PRIORS
struct GammaPrior {
	double a;
	double b;
};

struct ExpPrior {
	double lambda;
};

struct LogNormalPrior {
	double mu;
	double sd;
};

struct BetaPrior {
	double a, b;
};

struct BernPrior {
	double p;
};

// tQuotient()
// ------------------------------------------
// Calculate t distribution quotient for provided experiment and model.
double tQuotient(
	const SufficientStat &stat,  // observation statistics
	double model,  // model prediction
	const GammaPrior &var_prior)  // Gamma variance priors
{
	double a_post = var_prior.a + (stat.n - 1) / 2;
	double b_post = var_prior.b + stat.x_var * (stat.n - 1) / 2;
	double out = pow(1 + pow(stat.x_mean - model, 2) * stat.n / (2 * b_post), -a_post - 0.5);

	return(out);
}

// calcSufficientStat()
// -----------------------------
// Calculates sufficient statistics for
SufficientStat calcSufficientStat(
	const Replicates &obs,  // set of observations in experiment
	const vector<double> &beta_plate, 
	const vector<double> &beta_run)
{
	SufficientStat stat;
	stat.n = obs.size();

	// Estimate mean
	stat.x_mean = 0.0;

	for (int i = 0; i < obs.size(); i++) {
		// size_t nplate = obs[i].plate;
		// size_t nrun = obs[i].run;
		// double plate_offset = log(beta_plate[nplate]);
		// double run_offset = log(beta_run[nrun]);

		stat.x_mean += obs[i].log_rel_count;
		// stat.x_mean += plate_offset;
		// stat.x_mean += run_offset;
	}
	stat.x_mean /= obs.size();

	// Estimate variance
	stat.x_var = 0.0;

	if (obs.size() > 1) {
		for (int i = 0; i < obs.size(); i++) {
			// double plate_offset = log(beta_plate[obs[i].plate]);  // Can be optimized
			// double run_offset = log(beta_run[obs[i].run]);

			stat.x_var += pow(obs[i].log_rel_count - stat.x_mean, 2);
		}
		stat.x_var /= obs.size() - 1;
	}

	// Find concentration
	stat.conc = obs[0].concentration;  // assumes that concentration is the same for all observations

	return(stat);
}

// PRINT FUNCTIONS
// ---------------------------------------
template <typename T, size_t sz>
std::ostream& operator<<(std::ostream& os, const std::array<T, sz> &ar);  // declaration required for dependent operato<<

// operator<<(vector<T>)
// -------------------------------
// Writes content of vector as comma separated string. Requires T to have operator<< function.
// Vectors of vectors bottomns out, writing content one line at a time.
template<typename T>
std::ostream& operator<<(std::ostream& os, const vector<T> &v) {
	for (int i = 0; i < v.size() - 1; i++) {
		os << v[i] << ',';
	}
	os << v[v.size() - 1];  // last element
	return os;
}

// operator<<(array<T, sz>)
template <typename T, size_t sz>
std::ostream& operator<<(std::ostream& os, const std::array<T, sz> &ar) {
	for (int i = 0; i < sz - 1; i++) {
		os << ar[i] << ';';  // note auxilliary array print separator ';'
	}
	// last element
	os << ar[sz - 1];
	return(os);
}

std::ostream& operator<<(std::ostream& os, const SufficientStat &stat) {
	os << stat.n << ','
		<< stat.x_mean << ','
		<< stat.x_var;
	return os;
}

std::ostream& operator<<(std::ostream& os, const Observation &obs) {
	os << obs.experiment << ","
		<< obs.strain << ","
		<< obs.run << ","
		<< obs.plate << ","
		<< obs.pretreatment << ","
		<< obs.compound << ","
		<< obs.concentration << ","
		<< obs.log_rel_count;
	return os;
}

// Calculates response curve
double logResponse(double conc, 
	array<double, NPAR> params  // K, h, alpha
	)
{
	// double out = (1 - params[2]) / (1 + pow(params[0] * conc, params[1])) + params[2];
	double out = log((1 - params[2]) / (1 + pow(params[0] * conc, params[1])) + params[2]);
	return(out);
}

// Prior quotient correction of acceptance probability ratio. new->old ratio evaluation.
// Quotient correction from a Gamma prior
double qprior(double old_par, double new_par, GammaPrior gpr) {
	return(pow(new_par / old_par, gpr.a - 1) * pow(exp(old_par - new_par), gpr.b));
}

// Exponential prior.
double qprior(double old_par, double new_par, ExpPrior epr) {
	return(exp(epr.lambda * (old_par - new_par)));
}

// Log normal prior
double qprior(double old_par, double new_par, LogNormalPrior lnp) {
	double out = 
		(old_par / new_par) * 
		exp( 
			(1.0 / (2.0 * lnp.sd)) * 
			(
				pow(log(old_par) - log(lnp.mu), 2) -
				pow(log(new_par) - log(lnp.mu), 2)
			)
		);
	return(out);
}

double qprior(bool old_par, bool new_par, BernPrior bp) {
	double q;
	if (old_par && !new_par) {
		// 1 -> 0
		q = (1 - bp.p) / bp.p;
	} else if (!old_par && new_par) {
		// 0 -> 1
		q = bp.p / (1 - bp.p);
	} else {
		q = 1;
	}
	return(q);
}

// Beta prior
double qprior(double old_par, double new_par, BetaPrior bp) {
	double out = pow(new_par / old_par, bp.a - 1) * pow((1 - new_par) / (1 - old_par), bp.b - 1);
	return(out);
}

void printUsage(string program_name) {
	cout << "USAGE: " << program_name << " FILE"
		<< " <options>" << endl;
}

// options, set by <boost/program_options>
// Global variable
struct Options {
	int iterations;
	int burn;
	int subsample;  // store every nth parameter set.
	int init_phase;  // relies on init_lambda to be true.
	string strain;
	string data_file;
	string out_dir;  // output directory
} options;

int main(int argc, char const *argv[])
{
	// Parse commandline arguments
	po::options_description opts_desc("Options");

	try {
		// Configure command-line options
		po::variables_map vm;

		// Positional arguments
		opts_desc.add_options()
			("file",
				po::value<string>(&options.data_file)->required(),
				"Path to .tsv data file containing cell viability data in flat format -- one viability measurement per line.")
			;

		// Optional arguments
		opts_desc.add_options()
			("out,o",
				po::value<string>(&options.out_dir)->default_value("."),
				"Output directory to print sampled MCMC parameters.")
			("cell,c",
				po::value<std::string>(&options.strain),
				"Cell type to use from input data")
			("iterations,i",
				po::value<int>(&options.iterations)->default_value(500000),
				"Total MCMC iterations to evaluate.")
			("burn,b",
				po::value<int>(&options.burn)->default_value(100000),
				"MCMC burn iterations")
			("subsample,s",
				po::value<int>(&options.subsample)->default_value(200),
				"MCMC subsample every s iterations for output.")
			("init,n",
				po::value<int>(&options.init_phase)->default_value(20000),
				"Iterations used for initial choice of lambda.")
			("help,h", "Help screen")
			;

		po::positional_options_description p;
		p.add("file", -1);

		po::command_line_parser parser(argc, argv);
		parser.options(opts_desc).positional(p);
		po::parsed_options parsed_options = parser.run();

		po::store(parsed_options, vm);
		po::notify(vm);

		if (vm.count("help") || argc == 1) {
			printUsage(argv[0]);
			std::cout << opts_desc << std::endl;  // print options
			exit(EXIT_SUCCESS);
		}
	} catch (const po::error &e) {
		std::cerr << e.what() << std::endl;
		printUsage(argv[0]);
		cerr << opts_desc << endl;
		exit(EXIT_SUCCESS);
	} catch (const std::exception &e) {
		std::cerr << e.what() << std::endl;
		exit(EXIT_FAILURE);
	}

	// Print input arguments
	cout << argv[0] << " called with arguments:" << endl;
	cout << "\tfile: " << options.data_file << endl;
	cout << "\tout_dir: " << options.out_dir << endl;
	cout << "\tStrain: " << options.strain << endl;
	cout << "\tburn: " << options.burn << endl;
	cout << "\titerations: " << options.iterations << endl;
	cout << "\tsubsample: " << options.subsample << endl;
	cout << "\tinit_phase: " << options.init_phase << endl;

	// Initial values
	double init_beta = 1.0;
	bool init_lambda = true;
	bool init_lambda_AB = true;
	array<double, NPAR> init_theta = {0.01, 1.5, 0.5};  // K, h, alpha

	// Priors
	GammaPrior var_prior;
	var_prior.a = 0.6;  // worst case from statistics covering PANC1 and A375
	var_prior.b = 0.02;

	// Beta prior
	LogNormalPrior beta_prior;
	beta_prior.mu = 1.0;
	beta_prior.sd = 0.05;

	// K_prior;
	LogNormalPrior K_prior;
	K_prior.mu = 0.1;
	K_prior.sd = 2.0;

	LogNormalPrior h_prior;
	h_prior.mu = 1.5;  // log operation is performed subsequently
	h_prior.sd = 0.5;

	BetaPrior alpha_prior;
	alpha_prior.a = 1;
	alpha_prior.b = 3;  // was 3

	// Proposal distributions
	double switch_prop = 0.1;  // lambda switch proposals, was 0.05
	bernoulli_distribution flip(switch_prop);
	array<double, NPAR> prop_sd = {0.5, 0.1, 0.1};   // logN, N, N : K, h, alpha
	double beta_prop_sd = 0.1;

	// Parse data
	// -----------------------------------------------------
	std::cout << "Loading data..." << std::endl;

	set<string> drugs;  // Set of drugs to analyze
	map<string, Replicates> experiments;  // exp_id -> data triplicates
	int max_plate = 0;  // highest plate id
	int max_run = 0;  // highest run number
	{
		// Read whole data set into string vector
		std::vector<std::string> viability_data = readFile(options.data_file);
		cout << "observations: " << viability_data.size() - 1 << endl;

		cout << "Transforming data..." << endl;
		for (int i = 1; i < viability_data.size(); i++)  // first line is colnames
		{
			// Transform data line into vector of fields -- a data row.
			vector<string> row = split(viability_data[i], ',');
			Observation obs = parseRow(row);

			if (obs.strain != options.strain) continue;  // exclude strains not included in analysis

			// Include drugs for which AB data is available
			if (obs.experiment == "AB") {
				drugs.insert(obs.compound);
			}

			// Add data row to Experimental ID map, for calculation of 
			string exp_id = genExpID(row);
			experiments[exp_id].push_back(obs);

			// Find maximum plate and run id
			if (obs.plate > max_plate) {
				max_plate = obs.plate;  // new max
			}

			if (obs.run > max_run) {
				max_run = obs.run;
			}
		}
	}

	cout << "number of experiments: " << experiments.size() << endl;
	cout << "max run id: " << max_run << endl;
	cout << "max plate id: " << max_plate << endl;

	// Create drug indices, name->index. Create by iterating over 
	// the included drug names
	map<string, int> drug_index;
	vector<string> drugs_vec;  // drug name vector, matched to indices. index->name
	{
		int n;
		set<string>::iterator namep;
		for (namep = drugs.begin(), n = 0; namep != drugs.end(); ++namep, ++n) {
			drug_index[*namep] = n;
			drugs_vec.push_back(*namep);

		}
	}

	// Print drug list
	ofstream drugs_file(
		(fs::path(options.out_dir) / fs::path("drugs.csv"))  // constructs path
		.string());
	drugs_file << drugs_vec << endl;
	drugs_file.close();

	// Separate data into experiment specific observations.
	// Initialize data structures
	vector<ExperimentSet> obsA;  // obsA[a][e][r],  name,experiment,repeat
	obsA.resize(drugs.size());

	vector<ExperimentSet> obsA0;  // obsA0[a][e][r], although only single experiment
	obsA0.resize(drugs.size());

	vector<vector<ExperimentSet> > obsAB;  // obsAB[a][b][e][r]
	obsAB.resize(drugs.size());
	for (int i = 0; i < drugs.size(); i++) {
		obsAB[i].resize(drugs.size());
	}

	// Iterate over all experiments, separating each experiment into the structures obsA, obsA0, and obsAB.
	for (map<string, Replicates>::const_iterator replip = experiments.begin(); 
		replip != experiments.end(); 
		++replip) 
	{
		// replip->first is the experimentalID, replip->second is a vector of Observation objects.
		// Check if experiment is on list of drugs to analyze
		Observation first_obs = replip->second[0];
		string drug_name = replip->second[0].compound;
		set<string>::iterator it = drugs.find(drug_name);
		if (it == drugs.end()) continue;  // not found, next experiment

		int b = drug_index[first_obs.compound];  // drug index, assumed to be same for all Observations in set
		if (first_obs.experiment == "A") {
			obsA[b].push_back(replip->second);
		} else if (first_obs.experiment == "A0") {
			obsA0[b].push_back(replip->second);
		} else if (first_obs.experiment == "AB") {
			int a = drug_index[replip->second[0].pretreatment];
			obsAB[a][b].push_back(replip->second);
		} else {
			cerr << "ERROR: unrecognized experiment ID: " << first_obs.experiment << endl;
			exit(EXIT_FAILURE);
		}
	}

	cout << "number of drugs included in analysis: " << drugs.size() << endl;

	// MCMC analysis
	// ------------------------------------------------------------------ 
	cout << "Initializing MCMC..." << endl;

	ofstream theta_file(
		(fs::path(options.out_dir) / fs::path("theta.csv"))
		.string()
		);  // overwrites previous file by default
	ofstream theta_AB_file(
		(fs::path(options.out_dir) / fs::path("theta_AB.csv"))
		.string());
	ofstream lambda_file(
		(fs::path(options.out_dir) / fs::path("lambda.csv"))
		.string()
		);
	ofstream lambda_AB_file(
		(fs::path(options.out_dir) / fs::path("lambda_AB.csv"))
		.string()
		);
	ofstream beta_residual_file(
		(fs::path(options.out_dir) / fs::path("beta_residual.csv"))
		.string()
		);

	if (!theta_file.is_open() || 
		!theta_AB_file.is_open() || 
		!lambda_file.is_open() || 
		!lambda_AB_file.is_open() || 
		!beta_residual_file.is_open()) 
	{
		cout << "ERROR: unable to open output file." << endl;
		exit(EXIT_FAILURE);
	}

	// Model parameters
	vector<bool> lambda;
	lambda.resize(drugs.size(), init_lambda);

	vector<bool> lambda_AB;  // matrix
	lambda_AB.resize(drugs.size() * drugs.size(), init_lambda_AB);

	vector<array<double, NPAR> > theta;
	theta.resize(drugs.size(), init_theta);

	vector<array<double, NPAR> > theta_AB;
	theta_AB.resize(drugs.size() * drugs.size(), init_theta);

	vector<double> beta_residual;  // BetaA0
	beta_residual.resize(drugs.size(), init_beta);

	vector<double> beta_plate;
	beta_plate.resize(max_plate + 1, init_beta);
	vector<double> beta_run;
	beta_run.resize(max_run + 1, init_beta);

	// Random number generator
	default_random_engine generator;
	uniform_real_distribution<double> unif01(0.0, 1.0);

	cout << "Running MCMC analysis..." << endl;
	cout.precision(15);
	for (int iter = 0; iter < options.iterations; iter++) {
		if (iter % 100 == 0) {
			cout << "Iteration " << iter << " out of " << options.iterations << endl;
		}

		// Beta residual proposal
		// ------------------------------------------
		for (int a = 0; a < drugs.size(); a++) {
			normal_distribution<double> rnorm(beta_residual[a], beta_prop_sd);
			double new_beta_residual = rnorm(generator);
			if (new_beta_residual < 0.0) continue;  // reject proposal

			double ratio = 1.0;  // acceptance ratio

			// A0 experiment 
			// For each experiment involving beta_residual[a]
			for (int i = 0; i < obsA0[a].size(); i++) {
				SufficientStat stat = calcSufficientStat(obsA0[a][i], beta_plate, beta_run);

				double model_mean = log(beta_residual[a]);
				double new_model_mean = log(new_beta_residual);
				if (lambda[a]) {
					model_mean += logResponse(stat.conc, theta[a]);
					new_model_mean += logResponse(stat.conc, theta[a]);
				}  // else flat line, log(1) = 0, no change to model_mean

				double q = tQuotient(stat, new_model_mean, var_prior) / tQuotient(stat, model_mean, var_prior);
				ratio *= q;
			}

			// AB experiments involving beta_A0
			for (int b = 0; b < obsAB[a].size(); b++) {
				for (int i = 0; i < obsAB[a][b].size(); i++) {
					SufficientStat stat = calcSufficientStat(obsAB[a][b][i], beta_plate, beta_run);

					double model_mean = log(beta_residual[a]);
					double new_model_mean = log(new_beta_residual);

					if (lambda[a]) {
						model_mean += logResponse(1.0, theta[a]);
						new_model_mean += logResponse(1.0, theta[a]);
					}

					if (lambda_AB[a * drugs.size() + b]) {
						model_mean += logResponse(stat.conc, theta_AB[a * drugs.size() + b]);
						new_model_mean += logResponse(stat.conc, theta_AB[a * drugs.size() + b]);
					} else if (lambda[b]) {
						model_mean += logResponse(stat.conc, theta[b]);
						new_model_mean += logResponse(stat.conc, theta[b]);
					}

					double q = tQuotient(stat, new_model_mean, var_prior) / tQuotient(stat, model_mean, var_prior);

					ratio *= q;
				}
			}

			// Priors
			double q = qprior(beta_residual[a], new_beta_residual, beta_prior);  // take prior into acount
			ratio *= q;

			if (unif01(generator) < ratio) {
				// accept new value
				beta_residual[a] = new_beta_residual;
			}  // else keep old parameter
		}

		// Theta and lambda proposals
		// ------------------------------------------------------
		for (int a = 0; a < drugs.size(); a++) {
			// Propose lambda and theta for single drug.
			bool new_lambda;
			if (iter < options.init_phase) {
				new_lambda = true;  // Forced on in initial phase to settle starting parameters
			} else {
				if (flip(generator)) {
					// switch
					new_lambda = !lambda[a];
				} else {
					// keeo previous lambda
					new_lambda = lambda[a];
				}
			}

			// Propose new parameter set conditioned on lambda
			array<double, NPAR> new_theta;
			if ((!lambda[a] && new_lambda) ||  (lambda[a] && !new_lambda)) {
				// p2, p3 case: off-> on, on->off
				// propose initial curve
				// new_theta = init_theta;  // reset to initial
				new_theta = theta[a];  // maintain as latent, tests for lambda difference only
			} else if (!lambda[a] && !new_lambda) {
				// p4,
				continue;  // keep current, both lambda and theta
			} else {
				// p1
				// choose random parameter to change
				uniform_int_distribution<> dpar(0, NPAR - 1);
				int npar = dpar(generator);

				new_theta = theta[a];

				if (npar == 0) {
					// K, Log normal
					normal_distribution<double> rnorm(log(theta[a][0]), prop_sd[0]);
					new_theta[0] = exp(rnorm(generator));
				} else if (npar == 1) {
					// h, Normal
					normal_distribution<double> rnorm(theta[a][1], prop_sd[1]);
					new_theta[1] = rnorm(generator);
					if (new_theta[1] < 0) continue;
				} else if (npar == 2) {
					// alpha, Normal
					normal_distribution<double> rnorm(theta[a][2], prop_sd[2]);
					new_theta[2] = rnorm(generator);
					if (new_theta[2] < 0.0 || new_theta[2] > 1.0) continue;
				} else {
					cerr << "ERROR: unrecognized parameter choice" << endl;
					exit(EXIT_FAILURE);
				}
			}

			// Acceptance ratio
			// -----------------
			// Calculate ratio for (lambda, theta) -> (new_lambda, new_theta)
			double ratio = 1.0;

			// A0 experiment, likelihood contribution
			for (int i = 0; i < obsA0[a].size(); i++) {
				SufficientStat stat = calcSufficientStat(obsA0[a][i], beta_plate, beta_run);

				double model_mean = log(beta_residual[a]);
				double new_model_mean = model_mean;  // same

				if (lambda[a]) {
					model_mean += logResponse(stat.conc, theta[a]);
				}

				if (new_lambda) {
					new_model_mean += logResponse(stat.conc, new_theta);
				}

				double q = tQuotient(stat, new_model_mean, var_prior) / tQuotient(stat, model_mean, var_prior);
				ratio *= q;
			}

			// A experiment, likelihood contribution
			for (int i = 0; i < obsA[a].size(); i++) {
				SufficientStat stat = calcSufficientStat(obsA[a][i], beta_plate, beta_run);

				double model_mean = 0.0;
				double new_model_mean = 0.0;
				if (lambda[a]) {
					model_mean += logResponse(stat.conc, theta[a]);
				}

				if (new_lambda) {
					new_model_mean += logResponse(stat.conc, new_theta);
				}

				double q = tQuotient(stat, new_model_mean, var_prior) / tQuotient(stat, model_mean, var_prior);
				ratio *= q;
			}

			// Likelihood contribution from AxB combinations involving A.
			for (int b = 0; b < obsAB[a].size(); b++) {
				// AB experiment likelihood contribution, where the parameters change the pretreatment effect
				for (int i = 0; i < obsAB[a][b].size(); i++) {
					SufficientStat stat = calcSufficientStat(obsAB[a][b][i], beta_plate, beta_run);

					double model_mean = log(beta_residual[a]);
					double new_model_mean = log(beta_residual[a]);

					// first treatment
					if (lambda[a]) {
						model_mean += logResponse(1.0, theta[a]);
					}

					if (new_lambda) {
						new_model_mean += logResponse(1.0, new_theta);
					}

					// second treatment
					if (lambda_AB[a * drugs.size() + b]) {
						// both B component determined by the theta_AB parameters.
						model_mean += logResponse(stat.conc, theta_AB[a * drugs.size() + b]);
						new_model_mean += logResponse(stat.conc, theta_AB[a * drugs.size() + b]);
					} else if (lambda[b]) {
						model_mean += logResponse(stat.conc, theta[b]);
						new_model_mean += logResponse(stat.conc, theta[b]);
					}  // else baseline, log(1) = 0

					double q = tQuotient(stat, new_model_mean, var_prior) / tQuotient(stat, model_mean, var_prior);

					ratio *= q;
				}

				// BA Likelihood contributions, where the proposed parameters change the compound/post effect.
				for (int i = 0; i < obsAB[b][a].size(); i++) {
					SufficientStat stat = calcSufficientStat(obsAB[b][a][i], beta_plate, beta_run);

					double model_mean = log(beta_residual[b]);
					double new_model_mean = log(beta_residual[b]);

					// first treatment
					if (lambda[b]) {
						model_mean += logResponse(1.0, theta[b]);
						new_model_mean += logResponse(1.0, theta[b]);
					}

					// second treatment
					if (lambda_AB[b * drugs.size() + a]) {
						continue;  // changing theta[a] has no effect
					} else {
						// falls back on proposed parameters
						if (lambda[a]) {
							model_mean += logResponse(stat.conc, theta[a]);
						}

						if (new_lambda) {
							new_model_mean += logResponse(stat.conc, new_theta);
						}
					}

					double q = tQuotient(stat, new_model_mean, var_prior) / tQuotient(stat, model_mean, var_prior);

					ratio *= q;
				}
			}

			// Proposal correction
			double q = new_theta[0] / theta[a][0];  // log normal proposal correction for K.
			if (isnan(q)) {
				cerr << "WARNING: NaN quotient from log normal correction." << endl;
			} else {
				ratio *= q;
			}

			// Prior correction
			q = qprior(theta[a][0], new_theta[0], K_prior);
			q *= qprior(theta[a][1], new_theta[1], h_prior);
			q *= qprior(theta[a][2], new_theta[2], alpha_prior);

			ratio *= q;

			// Accept/reject
			if (unif01(generator) < ratio) {
				// accept
				theta[a] = new_theta;
				lambda[a] = new_lambda;
			}  // else keep old theta, lambda set
		}

		// Theta, lambda proposals for AB combined
		// --------------------------------------------------------
		for (int a = 0; a < drugs.size(); a++) {
			for (int b = 0; b < drugs.size(); b++) {
				bool new_lambda;
				if (iter < options.init_phase) {
					// new_lambda = false;  // model depends only on single drug components
					new_lambda = true;  // 
				} else {
					if (flip(generator)) {
						// switch
						new_lambda = !lambda_AB[a * drugs.size() + b];
					} else {
						// keep previous lambda
						new_lambda = lambda_AB[a * drugs.size() + b];
					}
				}

				// Propose new parameter set conditioned on lambda
				array<double, NPAR> new_theta_AB;
				if ((!lambda_AB[a * drugs.size() + b] && new_lambda) ||  (lambda_AB[a * drugs.size() + b] && !new_lambda)) {
					// p2, p3 case: off-> on, on->off
					// propose initial curve
					new_theta_AB = theta_AB[a * drugs.size() + b];  // maintain as latent, tests for lambda difference only
					// new_theta_AB = theta[b];  // revert to b
				} else if (!lambda_AB[a * drugs.size() + b] && !new_lambda) {
					// p4,
					continue;  // keep current, both lambda and theta
				} else {
					// p1
					// choose random parameter to change
					uniform_int_distribution<> dpar(0, NPAR - 1);
					// cout << dpar(generator) << endl;
					int npar = dpar(generator);

					new_theta_AB = theta_AB[a * drugs.size() + b];

					if (npar == 0) {
						// K, Log normal
						normal_distribution<double> rnorm(log(theta_AB[a * drugs.size() + b][0]), prop_sd[0]);
						new_theta_AB[0] = exp(rnorm(generator));
					} else if (npar == 1) {
						// h, Normal
						normal_distribution<double> rnorm(theta_AB[a * drugs.size() + b][1], prop_sd[1]);
						new_theta_AB[1] = rnorm(generator);
						if (new_theta_AB[1] < 0) continue;
					} else if (npar == 2) {
						// alpha, Normal
						normal_distribution<double> rnorm(theta_AB[a * drugs.size() + b][2], prop_sd[2]);
						new_theta_AB[2] = rnorm(generator);
						if (new_theta_AB[2] < 0.0 || new_theta_AB[2] > 1.0) continue;
					} else {
						cerr << "ERROR: unrecognized parameter choice" << endl;
						exit(EXIT_FAILURE);
					}
				}

				double ratio = 1.0;
				// Likelihood contribution from AB experiments
				for (int i = 0; i < obsAB[a][b].size(); i++) {
					SufficientStat stat = calcSufficientStat(obsAB[a][b][i], beta_plate, beta_run);

					double model_mean = log(beta_residual[a]);
					double new_model_mean = log(beta_residual[a]);

					// first treatment
					if (lambda[a]) {
						model_mean += logResponse(1.0, theta[a]);
						new_model_mean += logResponse(1.0, theta[a]);
					}

					// second treatment
					if (lambda_AB[a * drugs.size() + b]) {
						// both B component determined by the theta_AB parameters.
						model_mean += logResponse(stat.conc, theta_AB[a * drugs.size() + b]);
					} else if (lambda[b]) {
						model_mean += logResponse(stat.conc, theta[b]);
					}  // else baseline, log(1) = 0

					if (new_lambda) {
						new_model_mean += logResponse(stat.conc, new_theta_AB);
					} else if (lambda[b]) {
						new_model_mean += logResponse(stat.conc, theta[b]);
					}  // else baseline model

					double q = tQuotient(stat, new_model_mean, var_prior) / tQuotient(stat, model_mean, var_prior);

					ratio *= q;
				}

				// Proposal correction
				double q = new_theta_AB[0] / theta_AB[a * drugs.size() + b][0];  // log normal proposal correction for K.
				ratio *= q;

				// Prior corrections
				q = qprior(theta_AB[a * drugs.size() + b][0], new_theta_AB[0], K_prior);
				q *= qprior(theta_AB[a * drugs.size() + b][1], new_theta_AB[1], h_prior);
				q *= qprior(theta_AB[a * drugs.size() + b][2], new_theta_AB[2], alpha_prior);

				ratio *= q;

				// Accept/reject
				if (unif01(generator) < ratio) {
					// accept
					theta_AB[a * drugs.size() + b] = new_theta_AB;
					lambda_AB[a * drugs.size() + b] = new_lambda;
				}  // else keep old theta, lambda set
			}
		}

		// Store parameter set in files
		if (iter > options.burn && iter % options.subsample == 0) {
			beta_residual_file << beta_residual << endl;

			lambda_file << lambda << endl;
			lambda_AB_file << lambda_AB << endl;

			theta_file << theta << endl;
			theta_AB_file << theta_AB << endl;
		}
	}
	beta_residual_file.close();
	lambda_file.close();
	lambda_AB_file.close();
	theta_file.close();
	theta_AB_file.close();

	return 0;
}