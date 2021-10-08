//
// Created by frongere on 16/11/17.
//

#ifndef MATHUTILS_INTERP1D_H
#define MATHUTILS_INTERP1D_H

#include "Interp.h"
#include "Spline.h"

namespace mathutils {


  template<class XReal, class YScalar>
  class Interp1d {

   protected:
    std::shared_ptr<const std::vector<XReal>> xcoord;
    std::shared_ptr<const std::vector<YScalar>> yval;
    unsigned long ndata = 0;
    XReal xmin;
    XReal xmax;

   public:
    // TODO: voir a separer l'implementation et la mettre en fin de fichier (pas directement dans le corps de la classe)

    virtual void Initialize(std::shared_ptr<const std::vector<XReal>> x,
                            std::shared_ptr<const std::vector<YScalar>> y);

    virtual YScalar Eval(XReal x, bool permissive = false) const = 0;

    virtual std::vector<YScalar> Eval(const std::vector<XReal> &xvector, bool permissive = false) const = 0;

    YScalar operator()(const XReal x, bool permissive = false) const {
      return Eval(x, permissive);
    }

    std::vector<YScalar> operator()(const std::vector<XReal> &xvector, bool permissive = false) const {
      return Eval(xvector, permissive);
    }

    static Interp1d<XReal, YScalar> *MakeInterp1d(INTERP_METHOD method);

    XReal GetXmin() const { return xmin; }

    XReal GetXmax() const { return xmax; }

  };

  template<class XReal, class YScalar>
  void Interp1d<XReal, YScalar>::Initialize(std::shared_ptr<const std::vector<XReal>> x,
                                            std::shared_ptr<const std::vector<YScalar>> y) {

    assert(x->size() == y->size());
    assert (std::is_sorted(x->begin(), x->end()));

    ndata = x->size();
    xmin = x->at(0);
    xmax = x->at(ndata - 1);

    xcoord = x;
    yval = y;

  }


  template<class XReal, class YScalar>
  class Interp1dLinear : public Interp1d<XReal, YScalar> {

   private:
    std::vector<YScalar> a;
    std::vector<YScalar> b;

   public:

    void Initialize(std::shared_ptr<const std::vector<XReal>> x,
                    std::shared_ptr<const std::vector<YScalar>> y) override;

    YScalar Eval(XReal x, bool permissive=false) const override;

    std::vector<YScalar> Eval(const std::vector<XReal> &xvector, bool permissive=false) const override;

  };


  template<class XReal, class YScalar>
  class Interp1dBSpline : public Interp1d<XReal, YScalar> {
   private:
//        Spline m_spline;

   public:
    // TODO: terminer
  };


  template<class XReal, class YScalar>
  void Interp1dLinear<XReal, YScalar>::Initialize(const std::shared_ptr<const std::vector<XReal>> x,
                                                  const std::shared_ptr<const std::vector<YScalar>> y) {

    Interp1d<XReal, YScalar>::Initialize(x, y);

    a.reserve(this->ndata);
    b.reserve(this->ndata);

    XReal xi, xii, xii_m_xi;
    YScalar yi, yii;
    for (unsigned int i = 1; i < this->ndata; ++i) {

      xi = this->xcoord->at(i - 1);
      xii = this->xcoord->at(i);

      yi = this->yval->at(i - 1);
      yii = this->yval->at(i);

      xii_m_xi = xii - xi;

      a.push_back((yii - yi) / xii_m_xi);
      b.push_back((yi * xii - xi * yii) / xii_m_xi);

    }
  }

  template<class XReal, class YScalar>
  YScalar Interp1dLinear<XReal, YScalar>::Eval(const XReal x, bool permissive) const {
    // TODO: il faut que le type de retour soit compatible avec real et complex !!!

    double x_tmp = x;
    if (x < this->xmin or x > this->xmax) {
      if (permissive) {
        x < this->xmin ? x_tmp = this->xmin : x_tmp = this->xmax;
      } else {
        std::cerr << "Interpolation evaluated for value " << x << ", outside of the range : [" << this->xmin << ", "
                  << this->xmax << "]" << std::endl;
      }
    }
    // First, binary search on the x coords
    auto upper = std::lower_bound(this->xcoord->begin(), this->xcoord->end(), x_tmp);
    auto index = std::distance(this->xcoord->begin(), upper);

    if (index == 0) index = 1;  // Bug fix for x == xmin

    YScalar a_ = a.at(index - 1);
    YScalar b_ = b.at(index - 1);

    return a_ * x_tmp + b_;
  }

  template<class XReal, class YScalar>
  std::vector<YScalar> Interp1dLinear<XReal, YScalar>::Eval(const std::vector<XReal> &xvector, bool permissive) const {

    auto n = xvector.size();

    std::vector<YScalar> out;
    out.reserve(n);

    for (int i = 0; i < n; i++) {
      out.push_back(Eval(xvector[i], permissive));
    }
    return out;
  }

  /// Factory method to create 1D interpolation classes
  template<class XReal, class YScalar>
  Interp1d<XReal, YScalar> *Interp1d<XReal, YScalar>::MakeInterp1d(INTERP_METHOD method) {
    switch (method) {
      case LINEAR:
        return new Interp1dLinear<XReal, YScalar>;
        break;
      case BSPLINE:
//        return new Interp1dBSpline<XReal, YScalar>;
        std::cerr << "no B-Slpine interpolator available yet" << std::endl;
        break;
      default:
        std::cerr <<"1D INTERPOLATION METHOD DOES NOT EXIST" << std::endl;
        break;
    }
  }

}  // end namespace mathutils



#endif //MATHUTILS_INTERP1D_H
