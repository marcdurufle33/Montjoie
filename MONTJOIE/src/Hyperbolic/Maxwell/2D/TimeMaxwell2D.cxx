#ifndef MONTJOIE_FILE_TIME_MAXWELL_2D_CXX

#include "ProdMatVectMaxwell2D.cxx"

namespace Montjoie
{
   
  //! compute Y = G(tn,X) if we put time-maxwell equation to the form dX/dt = G(t,X)
  /*! 
    \param[in] tn input time
    \param[in] nb_deriv derivative number of G,
    if 0 we evaluate G, if 1 dG/dt, if 2 d^2 G /dt^2 ...
    \param[in] X iterate
    \param[out] Y result of G(tn,X)
  */
  template<class TypeEquation>
  void TimeMaxwell_2D<TypeEquation>::
  EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv, const VectReal_wp& X,
			     VectReal_wp& Y, bool invert_mass, bool source)
  {
    int nb_unknowns = this->GetNbScalarUnknowns();
    int nb_unknowns_vec = this->GetNbVectorialUnknowns();
    
    VectReal_wp En, Hn;
    VectReal_wp Prod_En, Prod_Hn;
    
    Real_wp one(1), zero(0);
    
    // X is split in (Uh, Vh, Vh_pml, Vh_diamond, Vh_star)
    // idem for Y
    Real_wp* x_ptr = X.GetData(); Real_wp* y_ptr = Y.GetData();
    En.SetData(nb_unknowns, x_ptr); Prod_En.SetData(nb_unknowns, y_ptr);
    x_ptr += nb_unknowns; y_ptr += nb_unknowns;
    
    Hn.SetData(nb_unknowns_vec, x_ptr); Prod_Hn.SetData(nb_unknowns_vec, y_ptr);
    x_ptr += nb_unknowns_vec; y_ptr += nb_unknowns_vec;
        
    // dE/dt = Dh^{-1} ( Rh^S H + Fh - Sh E)
    this->ApplyOperatorRhScalar(one, tn, Hn, zero, Prod_En);
    if (source)
      this->AddScalarSourceAtTime(one, tn, nb_deriv, Prod_En);
    
    this->ApplyOperatorSh(-one, tn, En, one, Prod_En);
    if (invert_mass)
      {
        this->SolveOperatorDh(Prod_En);
        this->SetDirichletCondition(tn, nb_deriv+1, Prod_En);
      }

    // dH/dt = Bh^{-1} ( Rh^V E - ShV H + FhV)  (for elements outside PML)
    this->ApplyOperatorRhVectorial(one, tn, En, zero, Prod_Hn);
    if (source)
      this->AddVectorialSourceAtTime(one, tn, nb_deriv, Prod_Hn);
    
    this->ApplyOperatorShVectorial(-one, tn, Hn, one, Prod_Hn);
    if (invert_mass)
      this->SolveOperatorBh(Prod_Hn);
    
    // nullify temporary vectors
    En.Nullify(); Hn.Nullify(); Prod_En.Nullify(); Prod_Hn.Nullify();
  }

  template<class TypeEquation>
  void TimeMaxwell_2D<TypeEquation>::
  EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv, const Real_wp& alpha, int level,
                             const VectReal_wp& X, VectReal_wp& Y, bool invert_mass, bool source)
  {
    abort();
  }
  
} // end namespace

#define MONTJOIE_FILE_TIME_MAXWELL_2D_CXX
#endif
  
