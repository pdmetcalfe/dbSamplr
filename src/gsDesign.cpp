#include <algorithm>
#include <exception>
#include <vector>

#include <Rcpp.h>

namespace {

  class bin_gen {
    public:
      bin_gen(double prob, int size) : prob(prob), size(size), ind(0) {};
      double operator()() {
        return R::dbinom(ind++, size, prob, 0);
      }
    private:
      double prob;
      int size;
      int ind;
  };

}

//' Simulate simple multiple testing strategy for single binomial event
//'
//' @param p event probability
//' @param sizes monotonically increasing vector of sample sizes
//' @param crits vector below which sampling stops
//' @return a vector giving probability of stopping at each sample
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector gsProbs(const double p,
                            const Rcpp::IntegerVector& sizes,
			    const Rcpp::IntegerVector& crits) {
  if (sizes.size() != crits.size()) {
    throw std::runtime_error("Size mismatch");
  }

  if (0 == sizes.size()) {
    throw std::runtime_error("Huh?");
  }
  
  if (p < 0 || p > 1) {
    throw std::runtime_error("Your probability should be between 0 and 1!");
  }
  
  // fill the probabilities
  std::vector<double> probs(1 + sizes[0]);
  std::generate(probs.begin(), probs.end(), bin_gen(p, sizes[0]));
  
  // here's where the result is allocated
  Rcpp::NumericVector result(sizes.size());

  result[0] = std::accumulate(probs.begin(),
			      probs.begin() + 1 + crits[0],
			      double(0));

  for (int i=1; i < sizes.size(); ++i) {
    // how many values are we getting in this step
    std::vector<double> stepProbs(1 + sizes[i] - sizes[i - 1]);
    std::generate(stepProbs.begin(), stepProbs.end(),
		  bin_gen(p, sizes[i] - sizes[i - 1]));
    std::vector<double> newProb(sizes[i]);
    for (int j=1 + crits[i-1]; j < probs.size(); ++j) {
      const int toDo = std::min(stepProbs.size(),
				newProb.size() - j);
#pragma omp parallel for
      for (int k=0;k < toDo; ++k) {
	newProb[j + k] += probs[j] * stepProbs[k];
      }
    }
    
    result[i] = std::accumulate(newProb.begin(),
				newProb.begin() + 1 + crits[i],
				double(0));
    std::swap(newProb, probs);
  }
  
  return result;
}
