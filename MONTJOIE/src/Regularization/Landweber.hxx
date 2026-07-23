#ifndef ITREG_FILE_LANDWEBER_HXX

namespace itreg
{
  
  template<class T>
  class LandweberP : public RegMethodP<T>
  {
  public:
    LandweberP();
    
    T omega; //!< scaling parameter

    void SetParameter(const string& keyword, const string& value);
    
  };


  template<class T>
  class ModLandweberP : public RegMethodP<T>
  {
  public:
    ModLandweberP();
    
    T omega; //!< scaling parameter
    T l0, psi;

    void SetParameter(const string& keyword, const string& value);
  };

  
  //! Scaled Landweber iteration
  /*!
    The step length parameter is chosen by estimating the norm
    of the Frechet derivative using a few steps of the power method
   */
  template<class T, class VectorRhs = Vector<T>, class VectorSol = Vector<T> >
  class ScaledLandweber : public RegularizationMethod<T, VectorRhs, VectorSol>
  {
  private:
    LandweberP<T>& param;
    
  public:
    ScaledLandweber(LandweberP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op);
    
    void Solve(const VectorRhs& ydelta, const T& delta, VectorSol& xn);
    virtual void Scaling(VectorSol& x, VectorRhs& y);
    
  protected:
    void OutputStep(int step, const VectorSol& xn, const T& residualNorm);

  };


  //! Landweber iteration
  template<class T, class VectorRhs = Vector<T>, class VectorSol = Vector<T> >
  class Landweber : public ScaledLandweber<T, VectorRhs, VectorSol>
  {
  private:
    LandweberP<T>& param;
    
  public:
    Landweber(LandweberP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op);
    
    virtual inline void Scaling(VectorSol& x, VectorRhs& y) {}
    
  };


  //! Modified Landweber iteration as suggested by O. Scherzer
  template<class T, class VectorRhs = Vector<T>, class VectorSol = Vector<T> >
  class ModLandweber : public RegularizationMethod<T, VectorRhs, VectorSol>
  {
  private:
    ModLandweberP<T>& param;
    
  public:
    ModLandweber(ModLandweberP<T>& param_, ForwardOperator<T, VectorRhs, VectorSol>& op);
    
    void Solve(const VectorRhs& ydelta, const T& delta, VectorSol& xn);
    T regPar(int k) const;
    
  protected:
    void OutputStep(int step, const VectorSol& xn, const T& residualNorm);

  };

}

#define ITREG_FILE_LANDWEBER_HXX
#endif
