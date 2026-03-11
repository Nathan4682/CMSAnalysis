#include <vector>
#include <algorithm>
#include "Math/DistFunc.h"
#include "TF1.h"
#include "CMSAnalysis/Analysis/interface/HiggsCompleteAnalysis.hh"
#include "CMSAnalysis/Analysis/interface/HistVariable.hh"

// Function definitions
// b is point at which you take KDE, vector is data set, h is bandwidth (you choose h!)
double boxcarKDE(double b, const std::vector<double> &samples, double h)
{
	// prevents error, returns 0 if bandwidth is invalid (0 or negative) or the data set is empty
	int n = samples.size();
	if (n == 0 || h <= 0.0)
		return 0.0;

	// counter and summation part of function
	int count = 0;
	for (double bi : samples)
	{
		if (std::abs(b - bi) <= h / 2.0)
			count++;
	}

	// 1/nh part of function before returning value
	return static_cast<double>(count) / (n * h);
}

double mu;
std::vector<double> backgroundSamples;
double h;
int Nobs;

double confusing(double *x, double *p)
{
	double b = x[0];
	double lambda = mu + b;
	// A root function that does poission cdf for you! awesome!
	double tail = 1.0 - ROOT::Math::poisson_cdf(Nobs - 1, lambda);
	return tail * boxcarKDE(b, backgroundSamples, h);
}

// function of stuff that works apparently
double marginalizedPoissonIntegral(double mu,
								   const std::vector<double> &backgroundSamples,
								   double h,
								   int Nobs)
{
	// Integration limits: min/max of background ± half bandwidth, saves some time so not calculating parts of integral where no data exists
	double bmin = *std::min_element(backgroundSamples.begin(), backgroundSamples.end()) - h / 2.0;
	double bmax = *std::max_element(backgroundSamples.begin(), backgroundSamples.end()) + h / 2.0;
	// remove lambda
	::mu = mu;
	::backgroundSamples = backgroundSamples;
	::h = h;
	::Nobs = Nobs;
	//  TF1 for integrand
	// double conf_val = confusing(mu, Nobs, h, backgroundSamples);
	TF1 integrand("integrand", confusing, bmin, bmax, 0);

	// Perform numerical integration
	return integrand.Integral(bmin, bmax);
}

double limitRootFunction(double *x, double *)
{
	double muVal = x[0];
	return marginalizedPoissonIntegral(muVal, ::backgroundSamples, ::h, ::Nobs) - 0.05;
}

void PdeWIP()
{
	HiggsCompleteAnalysis analysis;
	HistVariable histvar(HistVariable::VariableType::InvariantMass, "", true, false);
	TH1 *hist = analysis.getHist(histvar, "ZZ Background", false, "eeee");
	int nBins = hist->GetNbinsX();

	std::vector<double> backgroundSamplesLocal = {};
	double hLocal = 1.0;
	int NobsLocal = 7;

	for (int i = 1; i <= nBins; ++i)
	{
		double x_center = hist->GetBinCenter(i);
		double y_events = hist->GetBinContent(i);
		double y_error = hist->GetBinError(i);

		for (int j = 1; j <= y_events; ++j)
		{
			backgroundSamplesLocal.push_back(x_center);
		}
	}

	// Copy local data into the globals
	::backgroundSamples = backgroundSamplesLocal;
	::h = hLocal;
	::Nobs = NobsLocal;

	// Bounds of mu search is 0-2000
	TF1 f("f", limitRootFunction, 0, 2000, 0);

	double mu_limit = f.GetX(0.0, 0.0, 2000.0);

	std::cout << "Histogram " << hist->GetName()
			  << (hist->GetEntries() > 0 ? " has content." : " is empty.")
			  << " (Entries = " << hist->GetEntries() << ")"
			  << std::endl;

	std::cout << "95% CL upper limit on mu = " << mu_limit << std::endl;
}