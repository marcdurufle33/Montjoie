#ifndef MONTJOIE_FILE_GENERIC_EQUATION_HXX

namespace Montjoie
{
    
  //! base class to specify a new equation
  /*!
    Classes specifying equations are all derived from this base class
   */
  template<class T>
  class GenericEquation_Base
  {
  public :    
    //! equation solved with real or complex numbers
    typedef T Complexe;
    
    static const bool FirstOrderFormulation = false;
    static const bool DiscontinuousDiMatrix = false;
    static const bool TensorStiffnessSymmetric = true;
    
    GenericEquation_Base();
    
    //static void SetInputData(const string&, const Vector<string>&);
    static bool ComputeDFjm1();

    static bool SymmetricGlobalMatrix();
    static bool SymmetricElementaryMatrix();

    enum {type_element = 1 };
    static inline IVect GetOtherElementType() { return IVect(); }

    
    /**************
     * Mass terms *
     **************/
    
    template<class TypeEquation, class Dimension>
    static void ComputeMassMatrix(EllipticProblem<TypeEquation>& var,
                                  int i, const ElementReference_Dim<Dimension>& Fb);
    
    template<class TypeEquation, class T0, class MatMass>
    static void GetTensorMass(const EllipticProblem<TypeEquation>& vars,
			      int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
			      int ref, MatMass& mass);
    
    template<class TypeEquation, class T0, class Vector1>
    static void ApplyTensorMass(const EllipticProblem<TypeEquation>& var,
                                int i, int j, const GlobalGenericMatrix<T0>& nat_mat,
				int ref, Vector1& U, Vector1& V);
    
    /*******************
     * Stiffness terms *
     *******************/
    
    template<class TypeEquation, class NatureMat, class Vector1>
    static void GetNeededDerivative(const EllipticProblem<TypeEquation>& vars,
                                    const NatureMat& nat_mat,
                                    Vector1& unknown_to_derive, Vector1& fct_test_to_derive);
    
    template<class TypeEquation, class NatureMat, class MatStiff>
    static void GetGradGradTensor(const EllipticProblem<TypeEquation>& vars,
                                  int num_elem, int jloc, const NatureMat& nat_mat,
                                  int ref, MatStiff& Cgrad_grad);
    
    template<class TypeEquation, class NatureMat, class MatStiff>
    static void GetGradPhiTensor(const EllipticProblem<TypeEquation>& vars,
				  int num_elem, int jloc, const NatureMat& nat_mat, int ref,
				  MatStiff& Dgrad_phi, MatStiff& Ephi_grad);
    
    template<class TypeEquation, class NatureMat, class Vector1>
    static void ApplyTensorStiffness(const EllipticProblem<TypeEquation>& var,
				     int i, int j, const NatureMat& nat_mat,
                                     int ref, Vector1& dU, Vector1& dV);
    
    template<class TypeEquation, class NatureMat, class Vector1, class Vector2>
    static void ApplyGradientUnknown(const EllipticProblem<TypeEquation>& var,
				     int i, int j, const NatureMat& nat_mat,
                                     int ref, Vector1& dU, Vector2& V);

    template<class TypeEquation, class NatureMat, class Vector1, class Vector2>
    static void ApplyGradientFctTest(const EllipticProblem<TypeEquation>& var,
				     int i, int j, const NatureMat& nat_mat,
				     int ref, Vector1& U, Vector2& dV);
        
    /*************
     * PML terms *
     *************/

    template<class TypeEquation, class NatureMat, class Prop3, class Prop4, int p, int q>
    static void ApplyPmlGradGrad(const EllipticProblem<TypeEquation>& vars,
				  int num_elem, int jloc, const NatureMat& nat_mat, int ref,
				  TinyMatrix<TinyMatrix<Real_wp, Prop3, p, p>,
                                 Prop4, q, q>& Cgrad_grad);
    
    template<class TypeEquation, class NatureMat, class Prop3, class Prop4, int p, int q>
    static void ApplyPmlGradGrad(const EllipticProblem<TypeEquation>& vars,
				 int num_elem, int jloc, const NatureMat& nat_mat, int ref,
				 TinyMatrix<TinyMatrix<Complex_wp, Prop3, p, p>,
                                 Prop4, q, q>& Cgrad_grad);
    
    template<class TypeEquation, class NatureMat, class Prop3, int p, int q>
    static void ApplyPmlGradPhi(const EllipticProblem<TypeEquation>& vars,
				 int num_elem, int jloc, const NatureMat& nat_mat, int ref,
				 TinyMatrix<TinyVector<Real_wp, p>, Prop3, q, q>& Dgrad,
				 TinyMatrix<TinyVector<Real_wp, p>, Prop3, q, q>& Egrad);
    
    template<class TypeEquation, class NatureMat, class Prop3, int p, int q>
    static void ApplyPmlGradPhi(const EllipticProblem<TypeEquation>& vars,
				int num_elem, int jloc, const NatureMat& nat_mat, int ref,
				TinyMatrix<TinyVector<Complex_wp, p>, Prop3, q, q >& Dgrad,
				TinyMatrix<TinyVector<Complex_wp, p>, Prop3, q, q >& Egrad);

    template<class TypeEquation, class NatureMat, class Prop3, int p>
    static void ApplyPmlPhi(const EllipticProblem<TypeEquation>& vars,
                            int num_elem, int jloc, const NatureMat& nat_mat, int ref,
                            TinyMatrix<Real_wp, Prop3, p, p>& A);
    
    template<class TypeEquation, class NatureMat, class Prop3, int p>
    static void ApplyPmlPhi(const EllipticProblem<TypeEquation>& vars,
                            int num_elem, int jloc, const NatureMat& nat_mat, int ref,
                            TinyMatrix<Complex_wp, Prop3, p, p>& A);
    
    /**************
     * DFi impact *
     **************/
    
    template<class TypeEquation, class NatureMat, class MatStiff, class TypeElt>
    static void GetMassPhiDFiTensor(const EllipticProblem<TypeEquation>& vars,
				    int num_elem, int jloc, const NatureMat& nat_mat, int ref,
				    MatStiff& Amass, bool variable,
                                    bool affine, const TypeElt& Fb);
    
    template<class TypeEquation, class NatureMat, class MatStiff, class TypeElt>
    static void GetGradGradDFiTensor(const EllipticProblem<TypeEquation>& vars,
                                     int num_elem, int jloc, const NatureMat& nat_mat,
                                     int ref, MatStiff& Cgrad_grad,
                                     bool variable, bool affine, const TypeElt& Fb);
    
    template<class TypeEquation, class NatureMat,  class MatStiff, class TypeElt>
    static void GetGradPhiDFiTensor(const EllipticProblem<TypeEquation>& vars,
				    int num_elem, int jloc, const NatureMat& nat_mat, int ref,
				    MatStiff& Dgrad_phi, MatStiff& Ephi_grad,
                                    bool variable, bool affine, const TypeElt& Fb);
    
    /***********************
     * Boundary conditions *
     ***********************/
        
    template<class Matrix1, class TypeEquation, class NatureMat, class Dimension>
    static void GetNabc(Matrix1& Nabc, const typename Dimension::R_N& normale, int ref,
                        int iquad, int npoint, const NatureMat& nat_mat,
			int ref_d, const EllipticProblem<TypeEquation>& vars,
                        const ElementReference_Dim<Dimension>& Fb);
    
    template<class Vector1, class TypeEquation, class NatureMat, class Dimension>
    static void MltNabc(const typename Dimension::R_N& normale, int ref,
                        const Vector1& Vn, Vector1& Un,
			int num_elem1, int npoint, const NatureMat& nat_mat, int ref_d,
			const EllipticProblem<TypeEquation>& vars,
			const ElementReference_Dim<Dimension>& Fb);
    
    /************************************
     * Penalization/Stabilization Terms *
     ************************************/
        
    template<class Matrix1, class TypeEquation, class NatureMat, class Dimension>
    static void GetPenalDG(Matrix1& Nabc, const typename Dimension::R_N& normale, int iquad,
			   int npoint, int nf, const NatureMat& nat_mat, int ref, int ref2,
			   const EllipticProblem<TypeEquation>& vars,
                           const ElementReference_Dim<Dimension>& Fb);
    
    template<class Vector1, class TypeEquation, class NatureMat, class Dimension>
    static void MltPenalDG(const typename Dimension::R_N& normale, const Vector1& Vn, Vector1& Un,
			   int i, int npoint, int nf, const NatureMat& nat_mat, int ref, int ref2,
			   const EllipticProblem<TypeEquation>& vars,
                           const ElementReference_Dim<Dimension>& Fb);
    
  };
  
  
  //! class for a generic equation
  template <class T>
  class GenericEquation : public GenericEquation_Base<T>
  {
  };
  
} // namespace Montjoie

#define MONTJOIE_FILE_GENERIC_EQUATION_HXX
#endif
