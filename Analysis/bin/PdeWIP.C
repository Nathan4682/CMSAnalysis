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

double signalRegionMin;
double signalRegionMax;
double totalBackgroundNorm; // this is v_B

double intervalOverlap(double a1, double a2, double b1, double b2)
{
	double left = std::max(a1, b1);
	double right = std::min(a2, b2);
	return std::max(0.0, right - left);
}

// Computes int dm for the boxcar KDE
double backgroundFractionInSignalRegion(const std::vector<double> &massSamples,
										double h,
										double sMin,
										double sMax)
{
	int n = massSamples.size();
	if (n == 0 || h <= 0.0 || sMax <= sMin)
		return 0.0;

	double totalOverlap = 0.0;

	for (double m_i : massSamples)
	{
		double kernelMin = m_i - h / 2.0;
		double kernelMax = m_i + h / 2.0;
		totalOverlap += intervalOverlap(kernelMin, kernelMax, sMin, sMax);
	}

	return totalOverlap / (n * h);
}

// Computes b_s int dm
double expectedBackgroundYieldInSignalRegion(const std::vector<double> &massSamples,
											 double h,
											 double sMin,
											 double sMax,
											 double nuB)
{
	return nuB * backgroundFractionInSignalRegion(massSamples, h, sMin, sMax);
}

double integrand_call(double *x, double *p)
{
	(void)x;
	(void)p;

	double b_estimated =
		expectedBackgroundYieldInSignalRegion(backgroundSamples,
											  h,
											  signalRegionMin,
											  signalRegionMax,
											  totalBackgroundNorm);

	double lambda = mu + b_estimated;
	double tail = 1.0 - ROOT::Math::poisson_cdf(Nobs - 1, lambda);

	return tail;
}

// function of stuff that works apparently
double marginalizedPoissonIntegral(double mu,
								   const std::vector<double> &backgroundSamples,
								   double h,
								   int Nobs)
{
	::mu = mu;
	::backgroundSamples = backgroundSamples;
	::h = h;
	::Nobs = Nobs;

	return integrand_call(nullptr, nullptr);
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