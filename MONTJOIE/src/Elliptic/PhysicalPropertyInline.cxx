#ifndef MONTJOIE_FILE_PHYSICAL_PROPERTY_INLINE_CXX

#include "OneDimPhysicalIndexInline.cxx"

namespace Montjoie
{
  
  /*************************
   * PhysicalIndices_Sinus *
   *************************/
  
  
  //! default constructor
  template<class Dimension, class T>
  inline PhysicalSinusoidalMedia<Dimension, T>::PhysicalSinusoidalMedia()
    : RectangleCutOff<Dimension, Real_wp>()
  {
    SetComplexOne(mu0);
    SetComplexZero(amplitude);
    inverse = false;
  }
  
  
  //! returns true if some coefficients have an imaginary part
  template<class Dimension, class T>
  inline bool PhysicalSinusoidalMedia<Dimension, T>::IsComplex() const
  {
    if (abs(imagpart(mu0)) > ScalarPhysicalIndice<Dimension, T>::threshold)
      return true;
    
    if (abs(imagpart(mu0)) > ScalarPhysicalIndice<Dimension, T>::threshold)
      return true;
    
    return false;
  }
  
  
  //! user asks to inverse physical index
  template<class Dimension, class T>
  inline void PhysicalSinusoidalMedia<Dimension, T>::GetInverse()
  {
    inverse = !inverse;
  }
  
    
  //! kind of constructor
  template<class Dimension, class T>
  inline void PhysicalSinusoidalMedia<Dimension, T>::
  Init(const R_N& xmin_, const R_N & xmax_, const R_N& coef_,
       const T& offset, const T& ampli, const R_N& kwave_)
  {
    RectangleCutOff<Dimension, Real_wp>::Init(xmin_, xmax_, coef_);
    mu0 = offset;
    amplitude = ampli;
    kwave = kwave_;
    inverse = false;
  }
  
  
  //! multiplication of the physical index by a coefficient
  template<class Dimension, class T>
  inline void PhysicalSinusoidalMedia<Dimension, T>::Mlt(const T& coef)
  {
    amplitude *= coef;
    mu0 *= coef;
  }

  
  /***************************
   * PhysicalRegularMedia *
   ***************************/
  
  
  //! the media will be assumed to be periodic
  /*!
    In that case, we assume that the media is periodic, and values
    are given on the periodic cell [xmin, xmax] x [ymin, ymax] x [zmin, zmax]
   */
  template<class Dimension, class T>
  inline void PhysicalRegularMedia<Dimension, T>::SetPeriodic()
  {
    periodic = PERIODICITY;
  }

  
  //! the media will not be assumed to be periodic
  template<class Dimension, class T>
  inline void PhysicalRegularMedia<Dimension, T>::SetNoPeriodic()
  {
    periodic = USUAL;
  }

  
  //! the media will be assumed to be quasi-periodic
  /*!
    In that case, we assume that the media is periodic, and values
    are given on the periodic cell [xmin, xmax] x [ymin, ymax] x [zmin, zmax]
    Only central cell is non-affected (and index is equal to mu0)
  */  
  template<class Dimension, class T>
  inline void PhysicalRegularMedia<Dimension, T>::SetQuasiPeriodic()
  {
    periodic = QUASI_PERIODICITY;
  }
  
  
  //! returns true if some coefficients have an imaginary part
  template<class Dimension, class T>
  inline bool PhysicalRegularMedia<Dimension, T>::IsComplex() const
  {
    return complex_index;
  }
  
  
  //! user asks to inverse physical index
  template<class Dimension, class T>
  inline void PhysicalRegularMedia<Dimension, T>::GetInverse()
  {
    inverse = !inverse;
  }
  
  
  //! multiplication of the physical index by a coefficient
  template<class Dimension, class T>
  inline void PhysicalRegularMedia<Dimension, T>::Mlt(const T& coef)
  {
    amplitude *= coef;
    mu0 *= coef;
  }
  

  //! finds interval [imin, imax] containing the considered point x
  /*!
    \param[in] i nearest integer to considered point, with i < x
    \param[out] imin first extremity of interval
    \param[out] imax second extremity of interval
    \param[in] n imax has to be lower than n
    \param[in] dx difference x - i
    \return local coordinate in interval [imin, imax] (0 if x = imin and 1 if x = imax)
  */
  template<class Dimension, class T>
  inline Real_wp PhysicalRegularMedia<Dimension, T>::
  GetIndexBounds(int i, int& imin, int& imax,
		 int n, const Real_wp& dx) const
  {
    // for the indices, we want absolutely a continuity of the indice
    // so we don't take "neighboring cells to the point, but rather 
    // we consider that there is an underlying regular mesh, 
    // each element of the mesh containing several cells
    // depending on the order of interpolation
    Real_wp x(0);
    imin = 4*(i/4);
    imax = min(n-1, imin+4);
    if (imax < (imin+4))
      imin = imax-4;
    
    x = 0.25*(dx + i - imin);
    return x;
  }
  
    
  /***************************
   * PhysicalMeshMedia *
   ***************************/
  
  
  //! default constructor
  template<class Dimension, class T>
  inline PhysicalMeshMedia<Dimension, T>::PhysicalMeshMedia()
  {
    inverse = false; complex_index = false;
    mesh_order = 1;
  }
  
  
  //! returns true if the index is complex
  template<class Dimension, class T>
  inline bool PhysicalMeshMedia<Dimension, T>::IsComplex() const
  {
    return complex_index;
  }
  
  
  //! user asks to inverse physical index
  template<class Dimension, class T>
  inline void PhysicalMeshMedia<Dimension, T>::GetInverse()
  {
    inverse = !inverse;
  }
  
  
  //! multiplication of the physical index by a coefficient
  template<class Dimension, class T>
  inline void PhysicalMeshMedia<Dimension, T>::Mlt(const T& coef)
  {
    for (int i = 0; i < nu.GetM(); i++)
      Seldon::Mlt(coef, nu(i));
  }
  
  
  //! returns the file name where the mesh is stored
  template<class Dimension, class T>
  inline const string& PhysicalMeshMedia<Dimension, T>::GetMeshName() const
  {
    return mesh_file_name;
  }
  
  
  //! returns the order of approximation
  template<class Dimension, class T>
  inline int PhysicalMeshMedia<Dimension, T>::GetOrder() const
  {
    return mesh_order;
  }


  //! returns the number of elements of the mesh
  template<class Dimension, class T>
  inline int PhysicalMeshMedia<Dimension, T>::GetNbElt() const
  {
    return nu.GetM();
  }
  
  
  //! returns the number of nodal points of the element i
  template<class Dimension, class T>
  inline int PhysicalMeshMedia<Dimension, T>::GetNbNodes(int i) const
  {
    return nu(i).GetM();
  }
    
  
  //! returns the value of index on element i and nodal point j
  template<class Dimension, class T>
  inline const T& PhysicalMeshMedia<Dimension, T>::GetValue(int i, int j) const
  {
    return nu(i)(j);
  }
  
  
  //! returns the gradient of index on element i and nodal point j
  template<class Dimension, class T>
  inline const TinyVector<T, Dimension::dim_N>& PhysicalMeshMedia<Dimension, T>::
  GetGradient(int i, int j) const
  {
    return grad_nu(i)(j);
  }


  //! returns the row r of hessian of index on element i and nodal point j
  template<class Dimension, class T>
  inline TinyVector<T, Dimension::dim_N> PhysicalMeshMedia<Dimension, T>::
  GetHessian(int i, int j, int r) const
  {
    TinyVector<T, Dimension::dim_N> row;
    GetRow(hessian_nu(i)(j), r, row);
    return row;
  }


  /*************************
   * RadialVaryingMedia *
   *************************/


  //! Default constructor
  template<class Dimension, class T>
  inline RadialVaryingMedia<Dimension, T>
  ::RadialVaryingMedia() 
  {
    index_radial = NULL;
    inverse = false; complex = false; 
  }


  //! Destructor
  template<class Dimension, class T>
  inline RadialVaryingMedia<Dimension, T>::~RadialVaryingMedia()
  {
    Clear();
  }


  //! Copy constructor
  template<class Dimension, class T>
  inline RadialVaryingMedia<Dimension, T>
  ::RadialVaryingMedia(const RadialVaryingMedia<Dimension, T>& v)
  {
    index_radial = NULL;
    *this = v;
  }


  //! Assignment operator
  template<class Dimension, class T>
  inline RadialVaryingMedia<Dimension, T>&
  RadialVaryingMedia<Dimension, T>::operator=(const RadialVaryingMedia<Dimension, T>& v)
  {
    Clear();
    if (v.index_radial != NULL)
      index_radial = v.index_radial->GetDuplicate();

    inverse = v.inverse;
    complex = v.complex;
    return *this;
  }


  /*************************
   * PhysicalVaryingMedia *
   *************************/
  
  
  //! default constructor
  template<class Dimension, class T>
  inline PhysicalVaryingMedia<Dimension, T>::PhysicalVaryingMedia()
    : index_regular(), index_sinus(), index_mesh(), index_radial()
  {
    type = CONSTANT;
    SetComplexZero(cte_coef);
    same_mesh = false;
    SetComplexOne(amplitude_coef);
    SetComplexZero(offset_coef);
  }
  
  
  //! constructor needed when using TinyMatrix<PhysicalVaryingMedia>
  template<class Dimension, class T>
  inline PhysicalVaryingMedia<Dimension, T>::PhysicalVaryingMedia(const T& cte)
    : index_regular(), index_sinus(), index_mesh(), index_radial()
  {
    type = CONSTANT;
    cte_coef = T(cte);
  }

  
  //! constructor needed when using TinyMatrix<PhysicalVaryingMedia>
  template<class Dimension, class T>
  inline PhysicalVaryingMedia<Dimension, T>::PhysicalVaryingMedia(int cte)
    : index_regular(), index_sinus(), index_mesh(), index_radial()
  {
    type = CONSTANT;
    cte_coef = T(cte);
  }
  
  
  //! true if physical index is varying
  template<class Dimension, class T>
  inline bool PhysicalVaryingMedia<Dimension, T>::IsVarying() const
  {
    return (type != CONSTANT);
  }
  
  
  //! true if physical index is equal to 0
  template<class Dimension, class T>
  inline bool PhysicalVaryingMedia<Dimension, T>::IsZero() const
  {
    if (type == CONSTANT)
      if (cte_coef == T(0))
	return true;
    
    return false;
  }
  
  
  //! returns true if the same mesh is used for the discretization
  template<class Dimension, class T>
  inline bool PhysicalVaryingMedia<Dimension, T>::IsSameMesh() const
  {
    return same_mesh;
  }
  
  
  //! setting the type of media (constant or variable)
  template<class Dimension, class T>
  inline void PhysicalVaryingMedia<Dimension, T>::SetMediaType(int type_)
  {
    type = type_;
  }
  

  //! returning the type of media
  template<class Dimension, class T>
  inline int PhysicalVaryingMedia<Dimension, T>::GetMediaType() const
  {
    return type;
  }
  
  
  //! returns constant value of media (when constant)
  template<class Dimension, class T>
  inline T PhysicalVaryingMedia<Dimension, T>::GetConstant() const
  {
    return cte_coef;
  }


  template<class Dimension, class T>
  inline const T& PhysicalVaryingMedia<Dimension, T>::GetOffset() const
  {
    return offset_coef;
  }
  
  
  template<class Dimension, class T>
  inline const T& PhysicalVaryingMedia<Dimension, T>::GetAmplitude() const
  {
    return amplitude_coef; 
  }
  
  
  //! tells that the same mesh is used for the discretization
  template<class Dimension, class T>
  inline void PhysicalVaryingMedia<Dimension, T>::SetSameMesh()
  {
    same_mesh = true;
  }
  
  
  //! returning values of index on quadrature points (if variable index)
  template<class Dimension, class T>
  inline Vector<Vector<T> >& PhysicalVaryingMedia<Dimension, T>::GetValue()
  {
    return eval_coef;
  }
  
  
  //! returning gradient of index on quadrature points (if variable index)
  template<class Dimension, class T>
  inline Vector<Vector<TinyVector<T, Dimension::dim_N> > >&
  PhysicalVaryingMedia<Dimension, T>::GetGradient()
  {
    return grad_coef;
  }
  

  //! returns hessian of index on quadrature points (if variable index)
  template<class Dimension, class T>
  inline Vector<Vector<TinyMatrix<T, Symmetric, Dimension::dim_N, Dimension::dim_N> > >&
  PhysicalVaryingMedia<Dimension, T>::GetHessian()
  {
    return hessian_coef;
  }

  
  //! sets physical index of quadrature point (i,j) to a given value
  template<class Dimension, class T>
  inline void PhysicalVaryingMedia<Dimension, T>::
  SetCoefficient(int i, int j, const T& coef)
  {
    if (type != CONSTANT)
      eval_coef(i)(j) = coef;
  }
  
  
  //! sets a gradient to a given value
  template<class Dimension, class T>
  inline void PhysicalVaryingMedia<Dimension, T>::
  SetCoefGradient(int i, int j, const TinyVector<T, Dimension::dim_N>& grad)
  {
    if (type != CONSTANT)
      grad_coef(i)(j) = grad;
  }


  //! sets a row of the hessian to a given value
  template<class Dimension, class T>
  inline void PhysicalVaryingMedia<Dimension, T>::
  SetCoefHessian(int i, int j, int r, const TinyVector<T, Dimension::dim_N>& row)
  {
    if (type != CONSTANT)
      SetRow(row, r, hessian_coef(i)(j));
  }
  
  
  //! returns value of physical index on quadrature point (i,j)
  template<class Dimension, class T>
  inline const T& PhysicalVaryingMedia<Dimension, T>::
  GetCoefficient(int i, int j) const
  {
    if (type == CONSTANT)
      return cte_coef;
    
    return eval_coef(i)(j);
  }
  
  
  //! returns gradient of physical index on quadrature point (i,j)
  template<class Dimension, class T>
  inline const TinyVector<T, Dimension::dim_N>&
  PhysicalVaryingMedia<Dimension, T>::GetCoefGradient(int i, int j) const
  {
    if (type == CONSTANT)
      return zero_grad;
    
    return grad_coef(i)(j);
  }
  
  
  //! returns hessian of physical index on quadrature point (i,j)
  template<class Dimension, class T>
  inline const TinyMatrix<T, Symmetric, Dimension::dim_N, Dimension::dim_N>&
  PhysicalVaryingMedia<Dimension, T>::GetCoefHessian(int i, int j) const
  {
    if (type == CONSTANT)
      return zero_hess;
    
    return hessian_coef(i)(j);
  }
  
  
  /*************************
   *  ScalarPhysicalIndice *
   *************************/
  
  
  //! default constructor
  template<class Dimension, class T>
  inline ScalarPhysicalIndice<Dimension, T>::ScalarPhysicalIndice()
  {
    type_media = CONSTANT;
    SetComplexOne(cte_rho);
  }
  
  
  //! returns constant index
  template<class Dimension, class T>
  inline ScalarPhysicalIndice<Dimension, T>::operator T() const
  {
    return cte_rho;
  }
  
  
  //! returns true if constant index is complex
  template<class Dimension, class T>
  inline bool ScalarPhysicalIndice<Dimension, T>::IsZero() const
  {
    if (type_media == CONSTANT)
      {
	if (cte_rho == T(0))
	  return true;
        
        return false;
      }
    
    return fct_rho.IsZero();
  }
  
  
  //! returns true if constant index is complex
  template<class Dimension, class T>
  inline bool ScalarPhysicalIndice<Dimension, T>::IsComplex() const
  {
    if (IsVarying())
      return fct_rho.IsComplex();
    
    if (abs(imagpart(cte_rho)) > threshold)
      return true;
    
    return false;
  }
  
  
  //! returns true if physical index is varying
  template<class Dimension, class T>
  inline bool ScalarPhysicalIndice<Dimension, T>::IsVarying() const
  {
    return ((type_media != CONSTANT)&&(fct_rho.IsVarying()));
  }
  
  
  //! sets the type of media
  template<class Dimension, class T>
  inline void ScalarPhysicalIndice<Dimension, T>::SetMediaType(int type_)
  {
    fct_rho.SetMediaType(type_);
    if (type_ > 0)
      type_media = VARYING;
    else
      type_media = CONSTANT;
  }
  
  
  //! returns the type of media
  template<class Dimension, class T>
  inline int ScalarPhysicalIndice<Dimension, T>::GetMediaType() const
  {
    return fct_rho.GetMediaType();
  }
  
  
  //! sets physical index to 1
  template<class Dimension, class T>
  inline void ScalarPhysicalIndice<Dimension, T>::SetIdentity()
  {
    type_media = CONSTANT;
    cte_rho = T(1);
  }
  
  
  //! returns constant (if index is constant)
  template<class Dimension, class T>
  inline const T& ScalarPhysicalIndice<Dimension, T>::GetConstant() const
  {
    return cte_rho;
  }

  
  //! sets physical index to a constant value
  template<class Dimension, class T>
  inline void ScalarPhysicalIndice<Dimension, T>::SetConstant(const T& coef)
  {
    type_media = CONSTANT;
    cte_rho = coef;
  }
  
  
  //! sets physical index to zero
  template<class Dimension, class T>
  inline void ScalarPhysicalIndice<Dimension, T>::Zero()
  {
    type_media = CONSTANT;
    cte_rho = T(0);
  }
  
  
  //! physical index is inversed
  template<class Dimension, class T>
  inline void ScalarPhysicalIndice<Dimension, T>::GetInverse()
  {
    if (type_media == CONSTANT)
      {
	if (cte_rho != T(0))
	  cte_rho = T(1)/cte_rho;
      }
    else
      fct_rho.GetInverse();
  }
  
  
  //! multiplication of the physical index by a coefficient
  template<class Dimension, class T>
  inline void ScalarPhysicalIndice<Dimension, T>::Mlt(const T& coef)
  {
    if (type_media == CONSTANT)
      cte_rho *= coef;
    else
      fct_rho.Mlt(coef);
  }
  
  
  //! returns physical index of quadrature point (i,j)
  template<class Dimension, class T> 
  inline T ScalarPhysicalIndice<Dimension, T>::
  GetCoefficient(const VarPhysicalProblem& var, int i, int j) const
  {
    if (type_media == CONSTANT)
      return cte_rho;
    else
      {
	int i1 = var.ElementRho(i);
	return fct_rho.GetCoefficient(i1,j);
      }
    
    return cte_rho;
  }
  
  
  //! provides coefficient and its gradient
  template<class Dimension, class T> 
  inline void ScalarPhysicalIndice<Dimension, T>::
  GetCoefGradient(const VarPhysicalProblem& var, int i, int j,
                  T& coef, TinyVector<T, Dimension::dim_N>& grad) const
  {
    if (type_media != CONSTANT)
      {
	int i1 = var.ElementRho(i);
	coef = fct_rho.GetCoefficient(i1, j);
	grad = fct_rho.GetCoefGradient(i1, j);
      }
    else
      {
	coef = cte_rho;
	grad.Zero();
      }
  }
  

  template<class Dimension, class T> 
  inline void ScalarPhysicalIndice<Dimension, T>
  ::GetCoefHessian(const VarPhysicalProblem& var, int i, int j,
		   T& coef, TinyVector<T, Dimension::dim_N>& grad,
		   TinyMatrix<T, Symmetric, Dimension::dim_N, Dimension::dim_N>& hess)
  {
    if (type_media != CONSTANT)
      {
	int i1 = var.ElementRho(i);
	coef = fct_rho.GetCoefficient(i1, j);
	grad = fct_rho.GetCoefGradient(i1, j);
	hess = fct_rho.GetCoefHessian(i1, j);
      }
    else
      {
	coef = cte_rho;
	grad.Zero();
	hess.Zero();
      } 
  }
  
  
  //! reads the physical index in a data file
  /*!
    \param[in] var given problem
    \param[in] nb subscript to access array parameters
    \param[in] parameters values retrieved from data file
    \param[in] keyword keyword of the line of data file
   */
  template<class Dimension, class T>
  inline void ScalarPhysicalIndice<Dimension, T>::
  SetInputData(int& nb,
	       const VectString& parameters, const string& keyword)
  {
    this->SetInputVaryingMedia(nb, fct_rho, cte_rho, parameters);
    
    type_media = CONSTANT;
    if (fct_rho.IsVarying())
      type_media = VARYING;
  }
  
  
  //! displays informations about class ScalarPhysicalIndice
  template<class Dimension, class T>
  inline ostream& operator <<(ostream& out, const ScalarPhysicalIndice<Dimension, T>& var)
  {
    if (var.type_media == var.CONSTANT)
      out<<var.cte_rho;
    else
      out<<"varying media";
    
    return out;
  }
  
  
  /*************************
   * VectorPhysical Indice *
   *************************/
  
  
  //! default constructor
  template<class Dimension, int m, class T>
  inline VectorPhysicalIndice<Dimension, m, T>::VectorPhysicalIndice()
  {
    type_media = CONSTANT;
  }
  
  
  //! returns a component of the constant vector
  template<class Dimension, int m, class T>
  inline const T& VectorPhysicalIndice<Dimension, m, T>::operator()(int i) const
  {
    return cte_rho(i);
  }
  
  
  //! returns constant vector
  template<class Dimension, int m, class T>
  inline const TinyVector<T, m>&
  VectorPhysicalIndice<Dimension, m, T>::GetConstant() const
  {
    return cte_rho;
  }
  
  
  //! returns 0 for constant vector, 1 otherwise
  template<class Dimension, int m, class T>
  inline int VectorPhysicalIndice<Dimension, m, T>::GetMediaType() const
  {
    return type_media;
  }
  
  
  //! returns the type of varying media for a component of the vector
  template<class Dimension, int m, class T>
  inline int VectorPhysicalIndice<Dimension, m, T>::GetMediaType(int i) const
  {
    return fct_rho(i).GetMediaType();
  }
  
  
  //! sets all the components of the physical index to 1
  template<class Dimension, int m, class T> 
  inline void VectorPhysicalIndice<Dimension, m, T>::SetIdentity()
  {
    type_media = CONSTANT;
    for (int i = 0; i < m; i++)
      fct_rho(i).SetMediaType(fct_rho(i).CONSTANT);
    
    cte_rho.Fill(1);
  }
  
  
  //! sets physical index to a constant vector
  template<class Dimension, int m, class T> 
  inline void VectorPhysicalIndice<Dimension, m, T>::SetConstant(const TinyVector<T, m>& coef)
  {
    type_media = CONSTANT;
    for (int i = 0; i < m; i++)
      fct_rho(i).SetMediaType(fct_rho(i).CONSTANT);
    
    cte_rho = coef;
  }
  
  
  //! sets physical index to null vector
  template<class Dimension, int m, class T>
  inline void VectorPhysicalIndice<Dimension, m, T>::Zero()
  {
    type_media = CONSTANT;
    for (int i = 0; i < m; i++)
      fct_rho(i).SetMediaType(fct_rho(i).CONSTANT);
    
    cte_rho.Fill(0);
  }
  

  //! multiplication by coef
  template<class Dimension, int m, class T> 
  inline void VectorPhysicalIndice<Dimension, m, T>::Mlt(const T& coef)
  {
    if (type_media == CONSTANT)
      Seldon::Mlt(coef, cte_rho);
    else
      {
	for (int k = 0; k < m; k++)
	  fct_rho(k).Mlt(coef);
      }
  }
  
    
  //! returns physical index of quadrature point (i,j)
  template<class Dimension, int m, class T> 
  inline TinyVector<T, m> VectorPhysicalIndice<Dimension, m, T>::
  GetCoefficient(const VarPhysicalProblem& var, int i, int j) const
  {
    // DISP(type_media); DISP(CONSTANT);
    if (type_media == VARYING)
      {
	int i1 = var.ElementRho(i);
	TinyVector<T, m> rho;
	for (int k = 0; k < m; k++)
	  rho(k) = fct_rho(k).GetCoefficient(i1, j);
	
	return rho;
      }
    
    return cte_rho;
  }
  
    
  //! displays informations about class VectorPhysicalIndice
  template<class Dimension, int m, class T>
  inline ostream& operator <<(ostream& out, const VectorPhysicalIndice<Dimension, m, T>& var)
  {
    if (var.type_media == var.CONSTANT)
      out<<var.cte_rho;
    else
      out<<"varying media";
    
    return out;
  }
  
  
  /************************
   * TensorPhysicalIndice *
   ************************/
  
  
  //! default constructor
  template<class Dimension, int m, class T>
  inline TensorPhysicalIndice<Dimension, m, T>::TensorPhysicalIndice()
  {
    type_anisotropy = ISOTROPE;
    type_media = CONSTANT;
    cte_mu.SetIdentity();
  }
  
  
  //! returns constant index
  template<class Dimension, int m, class T>
  inline const TinyMatrix<T, Symmetric, m, m>&
  TensorPhysicalIndice<Dimension, m, T>::GetConstant() const
  {
    return cte_mu;
  }
  
  
  //! returns component of the constant index
  template<class Dimension, int m, class T>
  inline const T& TensorPhysicalIndice<Dimension, m, T>::operator()(int i, int j) const
  {
    return cte_mu(i,j);
  }
  
  
  //! returns the type of anisotropy (isotrope, orthotrope or anisotrope)
  template<class Dimension, int m, class T>
  inline int TensorPhysicalIndice<Dimension, m, T>::GetAnisotropy() const
  {
    return type_anisotropy;
  }
  
  
  //! sets the type of anisotropy (isotrope, orthotrope or anisotrope)
  template<class Dimension, int m, class T>
  inline void TensorPhysicalIndice<Dimension, m, T>::SetAnisotropy(int type)
  {
    type_anisotropy = type;
  }
  

  //! returns type of media for a component of the tensor
  template<class Dimension, int m, class T>
  inline int TensorPhysicalIndice<Dimension, m, T>::GetMediaType(int i, int j) const
  {
    return fct_mu(i, j).GetMediaType();
  }
  
  
  //! returns 0 if the tensor is constant, 1 otherwise
  template<class Dimension, int m, class T>
  inline int TensorPhysicalIndice<Dimension, m, T>::GetMediaType() const
  {
    return type_media;
  }
  
  
  //! sets the physical index to identity matrix
  template<class Dimension, int m, class T>
  inline void TensorPhysicalIndice<Dimension, m, T>::SetIdentity()
  {
    type_anisotropy = ISOTROPE;
    type_media = CONSTANT;
    for (int i = 0; i < m; i++)
      for (int j = i; j < m; j++)
        fct_mu(i, j).SetMediaType(fct_mu(i, j).CONSTANT);
    
    cte_mu.SetIdentity();
  }


  //! returns number of components of matrix
  template<class Dimension, int m, class T>
  inline int TensorPhysicalIndice<Dimension, m, T>::GetSize() const
  {
    return cte_mu.GetSize();
  }
  
  
  /*************************
   * ElasticPhysicalIndice *
   *************************/
  
  
  //! default constructor
  template<class Dimension, int m, class T>
  inline ElasticPhysicalIndice<Dimension, m, T>::ElasticPhysicalIndice()
  {
    type_media = CONSTANT;
    type_anisotropy = ISOTROPE;
    cte_C.FillIsotrope(T(1), T(1));
  }
  
  
  //! sets the type of media
  template<class Dimension, int m, class T>
  inline void ElasticPhysicalIndice<Dimension, m, T>::SetMediaType(int type_)
  {
    for (int i = 0; i < nb_components; i++)
      for (int j = i; j < nb_components; j++)
	fct_C(i, j).SetMediaType(type_);
    
    if (type_ > 0)
      type_media = VARYING;
    else
      type_media = CONSTANT;
  }
  
  
  //! returns the type of media
  template<class Dimension, int m, class T>
  inline int ElasticPhysicalIndice<Dimension, m, T>::GetMediaType() const
  {
    return fct_C(0, 0).GetMediaType();
  }
  
  
  //! returns the type of anisotropy
  template<class Dimension, int m, class T>
  inline int ElasticPhysicalIndice<Dimension, m, T>::GetAnisotropyType() const
  {
    return type_anisotropy;
  }
  
  
  //! returns component (i,j,k,l) of constant elastic tensor
  template<class Dimension, int m, class T>
  inline const T& ElasticPhysicalIndice<Dimension, m, T>::operator()(int i, int j,
                                                                     int k, int l) const
  {
    return cte_C(i, j, k, l);
  }
  
  
  //! sets physical index to null tensor
  template<class Dimension, int m, class T>
  inline void ElasticPhysicalIndice<Dimension, m, T>::Zero()
  {
    type_anisotropy = ISOTROPE;
    type_media = CONSTANT;
    cte_C.Fill(0);
  }
  
  
  //! sets physical index to isotrope tensor with lambda = mu = 1
  template<class Dimension, int m, class T>
  inline void ElasticPhysicalIndice<Dimension, m, T>::SetIdentity()
  {
    type_media = CONSTANT;
    type_anisotropy = ISOTROPE;
    cte_C.FillIsotrope(T(1), T(1));
  }
  
  
  //! returns the constant elastic tensor
  template<class Dimension, int m, class T>
  inline const TinySymmetricTensor<T, m>& 
  ElasticPhysicalIndice<Dimension, m, T>::GetConstant() const
  {
    return cte_C;
  }
  
  
  //! apply tensor to a "stress vector" (dV = C dU)
  template<class Dimension, int m, class T> template<class Vector1>
  inline void ElasticPhysicalIndice<Dimension, m, T>::
  MltMatrix(const VarPhysicalProblem& var, int num_elem, int jloc,
	    const TinyVector<Vector1, m>& dU, TinyVector<Vector1, m>& dV) const
  {
    if (type_media == VARYING)
      {
        MltMatrixVar(var.ElementRho(num_elem), jloc, indic, dU, dV);
        return;
      }
    
    if (type_anisotropy == ANISOTROPE)
      cte_C.Mlt(dU, dV);
    else
      cte_C.MltOrthotrope(dU, dV);
  }
  
  
  //! apply tensor to a "stress vector" (dV = C dU)
  template<class Dimension, int m, class T> template<class Vector1>
  inline void ElasticPhysicalIndice<Dimension, m, T>::
  MltVector(const VarPhysicalProblem& var, int num_elem, int jloc,
	    const Vector1& dU, Vector1& dV) const
  {
    if (type_media == VARYING)
      {
        MltVectorVar(var.ElementRho(num_elem), jloc, indic, dU, dV);
        return;
      }

    if (type_anisotropy == ANISOTROPE)
      cte_C.Mlt(dU, dV);
    else
      cte_C.MltOrthotrope(dU, dV);
  }
  
  
  //! returns true if the index is varying
  template<class Dimension, int m, class T>
  inline bool ElasticPhysicalIndice<Dimension, m, T>::IsVarying() const
  {
    return (type_media != CONSTANT);
  }
    
}

#define MONTJOIE_FILE_PHYSICAL_PROPERTY_INLINE_CXX
#endif

