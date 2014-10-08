#ifndef DiHadronCorrelationAnalyzer_EventShape_h
#define DiHadronCorrelationAnalyzer_EventShape_h
/** \class EventShape
 *
 * Utility to compute thrust value and axis from
 * a collection of candidates.
 *
 * Ported from BaBar implementation.
 *
 * The thrust axis is the vector T which maximises
 * the following expression:
 *
 *       sum_i=1...N ( | T . p_i | )
 *  t = --------------------------------- 
 *      sum_i=1...N ( (p_i . _p_i)^(1/2) )
 *
 * where p_i, i=1...N are the particle momentum vectors.
 * The thrust value is the maximum value of t.
 * 
 * The thrust axis has a two-fold ambiguity due to taking the
 * absolute value of the dot product. This computation returns
 * by convention a thurs axis with a positive component along 
 * the z-direction is positive. 
 *
 * The thrust measure the alignment of a collection of particles
 * along a common axis.  The lower the thrust, the more spherical
 * the event is.  The higher the thrust, the more jet-like the
 *
 */
#include <vector>
#include <TVector3.h>
#include <TLorentzVector.h>

class EventShape  {
public:
  /// spatial vector
  typedef TVector3 Vector;
  /// constructor from first and last iterators
  template<typename const_iterator>
  EventShape(const_iterator begin, const_iterator end) :
    thrust_(0), axis_(0, 0, 0), pSum_(0), 
    n_(end - begin), isTran_(0) {
    if (n_ == 0) return;
    for(const_iterator i = begin; i != end; ++i) {
      cands_.push_back(&*i);
    }
    init(cands_);
  } 
  /// thrust value (in the range [0.5, 1.0])
  double thrust() const { return thrust_; } 
  double sphericity() const;
  double recoil() const; 
  /// thrust axis (with magnitude = 1)
  const Vector& axis() const { return axis_; } 
  void SetTran () { isTran_=1; }

private:

  std::vector<Vector*> cands_;
  double thrust_;
  double sphericity_;
  Vector axis_;
  double pSum_;
  const unsigned int n_;
  bool isTran_;

  struct ThetaPhi {
    ThetaPhi(double t, double p) : theta( t ), phi( p ) { }
    double theta, phi;
  };
  double thrust(const Vector & theAxis) const; 
  ThetaPhi initialAxis() const;
  ThetaPhi finalAxis(ThetaPhi) const;
  Vector axis(double theta, double phi) const;
  Vector axis(const ThetaPhi & tp) const  {
    return axis(tp.theta, tp.phi);
  }
  void parabola(double & a, double & b, double & c, 
		const Vector &, const Vector &, const Vector &) const;
  void init(std::vector<Vector*>&);
};

#endif
