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
// function of stuff that works apparently
double marginalizedPoissonIntegral(double mu,
								   const std::vector<double> &backgroundSamples,
								   double h,
								   int Nobs)
{
	// Integration limits: min/max of background ± half bandwidth, saves some time so not calculating parts of integral where no data exists
	double bmin = *std::min_element(backgroundSamples.begin(), backgroundSamples.end()) - h / 2.0;
	double bmax = *std::max_element(backgroundSamples.begin(), backgroundSamples.end()) + h / 2.0;

	// TF1 for integrand
	TF1 integrand("integrand", [&](double *x, double *p)
				  {
			double b = x[0];
			double lambda = mu + b;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
			// A root function that does poission cdf for you! awesome!
			double tail = 1.0 - ROOT::Math::poisson_cdf(Nobs - 1, lambda);
			return tail * boxcarKDE(b, backgroundSamples, h); }, bmin, bmax, 0);

	// Perform numerical integration
	return integrand.Integral(bmin, bmax);
}
void PdeWIP()
{
	HiggsCompleteAnalysis analysis;
	HistVariable histvar(HistVariable::VariableType::RecoSameSignInvariantMass, "", true, false);
	TH1 *hist = analysis.getHist(histvar, "ZZ Background", false, "");
	// count numb bins
	int nBins = hist->GetNbinsX();

	std::cout << "Hist Name:    " << hist->GetName() << std::endl;
	std::cout << "Hist Title:   " << hist->GetTitle() << std::endl;
	std::cout << "Total Bins:   " << hist->GetNbinsX() << std::endl;
	std::cout << "Total Events: " << hist->GetEntries() << std::endl;
	std::cout << "-------------------" << std::endl;

	// ENTER BACKGROUND SAMPLES HERE!!!
	std::vector<double> backgroundSamples = {};
	double h = 1.0;
	int Nobs = 7;

	for (int i = 1; i <= nBins; ++i)
	{
		double x_center = hist->GetBinCenter(i);
		double y_events = hist->GetBinContent(i);
		double y_error = hist->GetBinError(i);

		std::cout << "Bin " << i
				  << ": Center = " << x_center
				  << ", Events = " << y_events
				  << " +/- " << y_error << std::endl;
		}

	// Use TF1 to define the function f(mu) = integral - 0.05
	TF1 f("f", [&](double *x, double *)
		  {
        double mu = x[0];
        return marginalizedPoissonIntegral(mu, backgroundSamples, h, Nobs) - 0.05; }, 0, 20, 0);

	// Use GetX to find root, does bisection and other stuff for you!
	double mu_limit = f.GetX(0.0); // initial guess = 0

	std::cout << "95% CL upper limit on mu = " << mu_limit << std::endl;
}