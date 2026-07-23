#ifndef MONTJOIE_FILE_MULTI_DIM_PHYSICAL_INDEX_CXX

#include "MultiDimPhysicalIndex.hxx"

namespace Montjoie
{

  //! Default constructor
  template<class T, class Dimension>
  VariableParameterND_Base<T, Dimension>::VariableParameterND_Base()
  {
    nb_points_grid_c = 0;
    order = 1;
  }


  //! Destructor
  template<class T, class Dimension>
  VariableParameterND_Base<T, Dimension>::~VariableParameterND_Base() 
  {
  }


  //! Returns the number of degrees of freedom used to discretize c
  template<class T, class Dimension>
  int VariableParameterND_Base<T, Dimension>::GetNbDof() const
  {
    return nb_points_grid_c; 
  }


  //! Returns the interpolation points 
  template<class T, class Dimension>
  const typename Dimension::VectR_N& VariableParameterND_Base<T, Dimension>
  ::GetInterpolationPoints() const
  {
    return all_points;
  }


  //! returns an object containing the discretization of a multi-dimensional parameter
  template<class T, class Dimension>
  VariableParameterND_Base<T, Dimension>* 
  VariableParameterND_Base<T, Dimension>::GetNewApproximation(DistributedProblem<Dimension>& var,
							      const Vector<string>& param)
  {
    if (param.GetM() <= 0)
      {
	cout << "Enter the discretization method for c and number of points" << endl;
	cout << "Current parameters are " << param << endl;
	abort();
      }

    Real_wp x0 = var.GetXmin(), xN = var.GetXmax();    
    Real_wp y0 = var.GetYmin(), yN = var.GetYmax();    
    //Real_wp z0 = var.GetZmin(), zN = var.GetZmax();    

    VariableParameterND_Base<T, Dimension>* var_index;
    
    if (param(0) == "RectangularBSpline")
      {
	if (param.GetM() <= 1)
	  {
	    cout << "Enter the number of intervals in x y and orders" << endl;
	    cout << "Current parameters are " << param << endl;
	    abort();
	  }

	RectangularBSplineIndex<T, Dimension>* spline_approx;
	spline_approx = new RectangularBSplineIndex<T, Dimension>();
	int nb_intervals_x = to_num<int>(param(1));
	int nb_intervals_y = to_num<int>(param(2));
	int order_x = to_num<int>(param(3));
	int order_y = to_num<int>(param(4));

	VectReal_wp points_x, points_y;
	Linspace(x0, xN, nb_intervals_x+1, points_x);
	Linspace(y0, yN, nb_intervals_y+1, points_y);
	
	spline_approx->SetApproximation(points_x, points_y, order_x, order_y);
	var_index = spline_approx;
      }
    else if (param(0) == "Uniform")
      {
	UniformParameterIndex<T, Dimension>* unif_approx;
	unif_approx = new UniformParameterIndex<T, Dimension>();
	var_index = unif_approx;
      }
    else if (param(0) == "SameMesh")
      {
	SameMeshIndex<T, Dimension>* mesh_approx;
	mesh_approx = new SameMeshIndex<T, Dimension>();
	int r = to_num<int>(param(1));
	bool discont = false;
	if (param.GetM() > 2)
	  if (param(2) == "Discontinuous")
	    discont = true;
	
	mesh_approx->SetApproximation(r, discont, var);
	var_index = mesh_approx;
      }
    else
      {
	cout << "Unknown approximation " << param(0) << endl;
	abort();
      }
    
    return var_index;
  }
  
  
  //! Projects a function to degrees of freedom
  template<class T, class Dimension>
  void VariableParameterND_Base<T, Dimension>::ProjectFunction(T (&f)(const R_N&), Vector<T>& projF)
  {
    // default discretization is assumed to be interpolatory
    projF.Reallocate(nb_points_grid_c);
    for (int i = 0; i < all_points.GetM(); i++)
      projF(i) = f(all_points(i));
  } 


  //! returns weight for each basis function
  template<class T, class Dimension>
  void VariableParameterND_Base<T, Dimension>::GetWeightFunction(VectReal_wp& poids)
  {
    poids.Reallocate(GetNbDof());
    poids.Fill(Real_wp(1));
  }

  
  //! only dofs associated with elements such that num(i) is true are kept
  template<class T, class Dimension>
  void VariableParameterND_Base<T, Dimension>::KeepElements(const Vector<bool>& num)
  {
  }
  
  
  /*****************
   * SameMeshIndex *
   *****************/


  //! Default constructor
  template<class T, class Dimension>
  SameMeshIndex<T, Dimension>::SameMeshIndex()
  {
    discontinuous = true; 
  }


  //! Sets values of parameter on degrees of freedom
  template<class T, class Dimension>
  void SameMeshIndex<T, Dimension>::Init(const Vector<T>& cn)
  {
    val_cn = cn;
  }
  

  //! Initializes a parameter with the same mesh as contained in the given problem
  template<class T, class Dimension>
  void SameMeshIndex<T, Dimension>::SetApproximation(int order, bool discont, DistributedProblem<Dimension>& var)
  {
    mesh = var.mesh;
    mesh.SetGeometryOrder(max(1, order));
    this->order = order;
    discontinuous = discont;

    if (discontinuous)
      {
	offset.Reallocate(mesh.GetNbElt()+1);
	if (order == 0)
	  {
	    offset.Fill();
	    this->nb_points_grid_c = mesh.GetNbElt();
	    this->all_points.Reallocate(this->nb_points_grid_c);
	    R_N center; 	Vector<R_N> s;
	    for (int i = 0; i < mesh.GetNbElt(); i++)
	      {
		mesh.GetVerticesElement(i, s);
		center.Zero();
		for (int j = 0; j < s.GetM(); j++)
		  center += s(j);
		
		center /= s.GetM();
		this->all_points(i) = center;
	      }
	    
	    return;
	  }

	offset(0) = 0;
	for (int i = 0; i < mesh.GetNbElt(); i++)
	  offset(i+1) = offset(i) + mesh.GetNbPointsNodalElt(i);
	
	this->nb_points_grid_c = offset(mesh.GetNbElt());
	this->all_points.Reallocate(this->nb_points_grid_c);
	SetPoints<Dimension> PtsElem;
	Vector<R_N> s;
	int num = 0;
	for (int i = 0; i < mesh.GetNbElt(); i++)
	  {
	    mesh.GetVerticesElement(i, s);
	    mesh.FjElemNodal(s, PtsElem, mesh, i);
	    for (int j = 0; j < PtsElem.GetNbPointsNodal(); j++)
	      this->all_points(num++) = PtsElem.GetPointNodal(j);
	  }
      }
    else
      {
	mesh.GetNodesCurvedMesh(this->all_points, Nodle, order, true);
	this->nb_points_grid_c = this->all_points.GetM();
      }
    
    ElemIndex.Reallocate(mesh.GetNbElt());
    ElemIndex.Fill();
  }


  //! Evaluates the parameter at a given point
  template<class T, class Dimension>
  T SameMeshIndex<T, Dimension>::Evaluate(const R_N& pt_glob, int ne, const R_N& pt_loc)
  {
    int ne_loc = ElemIndex(ne);
    T val; SetComplexZero(val);
    if (ne_loc < 0)
      return val;
    
    if ((discontinuous) && (this->order == 0))
      return val_cn(offset(ne_loc));
    
    VectReal_wp phi;
    mesh.ComputeValuesPhiNodalRef(ne, pt_loc, phi);

    if (discontinuous)
      for (int j = 0; j < phi.GetM(); j++)
	val += phi(j) * val_cn(offset(ne_loc) + j);
    else
      for (int j = 0; j < phi.GetM(); j++)
	val += phi(j) * val_cn(Nodle(ne_loc)(j));      
    
    return val;
  }


  //! Evaluates basis functions at a given point
  template<class T, class Dimension>
  void SameMeshIndex<T, Dimension>
  ::ComputeValuesPhiRef(int i, int jloc, const R_N& point, 
			const R_N& point_loc, Vector<Real_wp, VectSparse>& phi)
  {
    int iloc = ElemIndex(i);
    if (iloc < 0)
      {
	phi.Clear();
	return;
      }
    
    if ((discontinuous) && (this->order == 0))
      {
	phi.Reallocate(1);
	phi.Index(0) = iloc;
	phi.Value(0) = 1.0;
	return;
      }

    VectReal_wp phi_loc;
    mesh.ComputeValuesPhiNodalRef(i, point_loc, phi_loc);
    phi.Reallocate(phi_loc.GetM());
    if (discontinuous)
      for (int j = 0; j < phi_loc.GetM(); j++)
	{
	  phi.Index(j) = offset(iloc) + j;
	  phi.Value(j) = phi_loc(j);
	}
    else
      for (int j = 0; j < phi_loc.GetM(); j++)
	{
	  phi.Index(j) = Nodle(iloc)(j);
	  phi.Value(j) = phi_loc(j);
	}
  }

  
  //! returns the numbers of degrees of freedom located on edge num_loc of element num_elem 
  template<class T, class Dimension>
  IVect SameMeshIndex<T, Dimension>::GetDofsOnBoundary(int ne, int num_loc)
  {
    IVect num;
    int ne_loc = ElemIndex(ne);
    if (ne_loc < 0)
      return num;
    
    num.Reallocate(this->order+1);
    int nv = mesh.Element(ne).GetNbVertices();
    num(0) = Nodle(ne_loc)(num_loc);
    num(this->order) = Nodle(ne_loc)((num_loc+1)%nv);
    for (int i = 1; i < this->order; i++)
      num(i) = Nodle(ne_loc)(nv + num_loc*(this->order-1) + i);
    
    return num;
  }


  //! returns weight for each basis function
  template<class T, class Dimension>
  void SameMeshIndex<T, Dimension>::GetWeightFunction(VectReal_wp& poids)
  {
    poids.Reallocate(this->GetNbDof());
    poids.Zero();
    typename Dimension::VectR_N s;
    SetPoints<Dimension> PtsElem; SetMatrices<Dimension> MatElem;
    for (int i = 0; i < mesh.GetNbElt(); i++)
      if (ElemIndex(i) >= 0)
	{
	  int iloc = ElemIndex(i);
	  mesh.GetVerticesElement(i, s);
	  mesh.FjElemNodal(s, PtsElem, mesh, i);
	  mesh.DFjElemNodal(s, PtsElem, MatElem, mesh, i);
	  const ElementGeomReference<Dimension>& Fb = mesh.GetReferenceElement(i);	  
	  if (discontinuous && (this->order == 0))
	    {
	      poids(iloc) = Det(MatElem.GetPointNodal(0));
	      continue;
	    }
	  
	  for (int j = 0; j < PtsElem.GetNbPointsNodal(); j++)
	    {
	      Real_wp Ji = Det(MatElem.GetPointNodal(j));
	      if (discontinuous)
		poids(offset(iloc) + j) += Ji*Fb.WeightsNodalND(j);
	      else
		poids(Nodle(iloc)(j)) += Ji*Fb.WeightsNodalND(j);
	    }
	}
  }


  //! only dofs associated with elements contained in num are kept
  template<class T, class Dimension>
  void SameMeshIndex<T, Dimension>::KeepElements(const Vector<bool>& num)
  {
    ElemIndex.Reallocate(mesh.GetNbElt());
    ElemIndex.Fill(-1);
    int nb = 0;
    for (int i = 0; i < ElemIndex.GetM(); i++)
      if (num(i))
	ElemIndex(i) = nb++;
    
    if (nb == mesh.GetNbElt())
      return;
    
    if (discontinuous)
      {
	IVect old_offset(offset);
	Vector<R_N> old_points(this->all_points);

	offset.Reallocate(nb+1);
	offset(0) = 0;
	for (int i = 0; i < mesh.GetNbElt(); i++)
	  if (ElemIndex(i) >= 0)
	    {
	      int iloc = ElemIndex(i);
	      offset(iloc+1) = offset(iloc) + mesh.GetNbPointsNodalElt(i);
	    }
	
	this->nb_points_grid_c = offset(nb);
	this->all_points.Reallocate(this->nb_points_grid_c);
	for (int i = 0; i < mesh.GetNbElt(); i++)
	  if (ElemIndex(i) >= 0)
	    {
	      int iloc = ElemIndex(i);
	      int nb_pts = mesh.GetNbPointsNodalElt(i);
	      for (int j = 0; j < nb_pts; j++)
		this->all_points(offset(iloc) + j) = old_points(old_offset(i) + j);
	    }	
      }
    else
      {
	Vector<IVect> OldNodle(Nodle);
	Nodle.Reallocate(nb);
	
	IVect IndexDof(this->nb_points_grid_c);
	IndexDof.Fill(-1);
	this->nb_points_grid_c = 0;
	for (int i = 0; i < mesh.GetNbElt(); i++)
	  if (ElemIndex(i) >= 0)
	    for (int j = 0; j < OldNodle(i).GetM(); j++)
	      {
		int n = OldNodle(i)(j);
		if (IndexDof(n) == -1)
		  {
		    IndexDof(n) = 0;
		    this->nb_points_grid_c++;
		  }
	      }
	
	nb = 0;
	Vector<R_N> old_points(this->all_points);
	this->all_points.Reallocate(this->nb_points_grid_c);
	for (int i = 0; i < IndexDof.GetM(); i++)
	  if (IndexDof(i) >= 0)
	    {
	      this->all_points(nb) = old_points(i);
	      IndexDof(i) = nb++;
	    }
	
	for (int i = 0; i < mesh.GetNbElt(); i++)
	  if (ElemIndex(i) >= 0)
	    {
	      int iloc = ElemIndex(i);
	      Nodle(iloc).Reallocate(OldNodle(i).GetM());
	      for (int j = 0; j < OldNodle(i).GetM(); j++)
		Nodle(iloc)(j) = IndexDof(OldNodle(i)(j));
	    }
      }
    
    //DISP(this->nb_points_grid_c); DISP(this->all_points);
    //DISP(ElemIndex); DISP(Nodle); DISP(offset);
  }
  
  

  /*************************
   * UniformParameterIndex *
   *************************/


  //! Default constructor
  template<class T, class Dimension>
  UniformParameterIndex<T, Dimension>::UniformParameterIndex()
  {
    SetComplexZero(val_cn);
    this->nb_points_grid_c = 1;
    this->all_points.Reallocate(1);
    this->all_points(0).Zero();
    this->order = 0;
  }


  //! Sets values of parameter on degrees of freedom
  template<class T, class Dimension>
  void UniformParameterIndex<T, Dimension>::Init(const Vector<T>& cn)
  {
    val_cn = cn(0);
  }
  

  //! Evaluates the parameter at a given point
  template<class T, class Dimension>
  T UniformParameterIndex<T, Dimension>::Evaluate(const R_N& pt_glob, int ne, const R_N& pt_loc)
  {
    return val_cn;
  }


  //! Evaluates basis functions at a given point
  template<class T, class Dimension>
  void UniformParameterIndex<T, Dimension>
  ::ComputeValuesPhiRef(int i, int jloc, const R_N& point, 
			const R_N& point_loc, Vector<Real_wp, VectSparse>& phi)
  {
    phi.Reallocate(1);
    phi.Index(0) = 0;
    phi.Value(0) = 1.0;
  }

  
  //! returns weight for each basis function
  template<class T, class Dimension>
  void UniformParameterIndex<T, Dimension>::GetWeightFunction(VectReal_wp& poids)
  {
    poids.Reallocate(this->GetNbDof());
    poids.Fill(1);
  }


  /***************************
   * RectangularBSplineIndex *
   ***************************/


  //! Sets values of c on degrees of freedom
  template<class T, class Dimension>
  void RectangularBSplineIndex<T, Dimension>::Init(const Vector<T>& cn)
  {
    val_cn = cn;
  }
  

  //! Sets the knots in x and y and orders
  template<class T, class Dimension>
  void RectangularBSplineIndex<T, Dimension>
  ::SetApproximation(const VectReal_wp& points_x, const VectReal_wp& points_y, int order_x, int order_y)
  {
    this->order = max(order_x, order_y);
    
    VectReal_wp knots_x, knots_y;
    knots_x.Reallocate(points_x.GetM() + 2*order_x);
    for (int i = 0; i < order_x; i++)
      {
	knots_x(i) = points_x(0);
	knots_x(points_x.GetM() + order_x + i) = points_x(points_x.GetM()-1);
      }
    
    for (int i = 0; i < points_x.GetM(); i++)
      knots_x(order_x+i) = points_x(i);

    knots_y.Reallocate(points_y.GetM() + 2*order_y);
    for (int i = 0; i < order_y; i++)
      {
	knots_y(i) = points_y(0);
	knots_y(points_y.GetM() + order_y + i) = points_y(points_y.GetM()-1);
      }
    
    for (int i = 0; i < points_y.GetM(); i++)
      knots_y(order_x+i) = points_y(i);
    
    Vector<T> Ones_x(knots_x.GetM() - order_x - 1); Ones_x.Fill(1);
    spline_x.SetData(knots_x, Ones_x, order_x);

    Vector<T> Ones_y(knots_y.GetM() - order_y - 1); Ones_y.Fill(1);
    spline_y.SetData(knots_y, Ones_y, order_y);

    this->nb_points_grid_c = spline_x.GetNbDof() * spline_y.GetNbDof();
    value_phi_x.Reallocate(spline_x.GetNbDof());
    value_phi_y.Reallocate(spline_y.GetNbDof());
    value_phi_x.Zero(); value_phi_y.Zero();

    // minimization with regular points
    Real_wp x0 = points_x(0);
    Real_wp xN = points_x(points_x.GetM()-1);
    Real_wp y0 = points_y(0);
    Real_wp yN = points_y(points_y.GetM()-1);

    VectReal_wp points_min_x, points_min_y;
    Linspace(x0, xN, points_x.GetM() + 3*order_x+2, points_min_x);
    Linspace(y0, yN, points_y.GetM() + 3*order_y+2, points_min_y);

    this->all_points.Reallocate(points_min_x.GetM()*points_min_y.GetM());
    int num = 0;
    for (int i = 0; i < points_min_x.GetM(); i++)
      for (int j = 0; j < points_min_y.GetM(); j++)
	{
	  this->all_points(num).Init(points_min_x(i), points_min_y(j));
	  num++;
	}
  }

  
  //! Evaluates the index on a given point
  template<class T, class Dimension>  
  T RectangularBSplineIndex<T, Dimension>::Evaluate(const R_N& pt_glob, int ne, const R_N& pt_loc)
  {
    Vector<Real_wp, VectSparse> vec_phi;
    ComputeValuesPhiRef(-1, -1, pt_glob, pt_loc, vec_phi);
    T val; SetComplexZero(val);
    for (int i = 0; i < vec_phi.GetM(); i++)
      val += vec_phi.Value(i)*val_cn(vec_phi.Index(i));
    
    return val;
  }


  //! Projects a function f to degrees of freedom
  template<class T, class Dimension>
  void RectangularBSplineIndex<T, Dimension>::ProjectFunction(T (&f)(const R_N&), Vector<T>& projF)
  {
    // VDM is computed
    Matrix<Real_wp> VDM;
    VDM.Reallocate(this->all_points.GetM(), this->GetNbDof());
    Vector<Real_wp, VectSparse> phi;
    VDM.Zero();
    for (int i = 0; i < this->all_points.GetM(); i++)
      {
	ComputeValuesPhiRef(-1, -1, this->all_points(i), R_N(), phi);
	for (int j = 0; j < phi.GetM(); j++)
	  VDM(i, phi.Index(j)) = phi.Value(j);
      }

    // Least-squares problem is solved with QR factorization
    VectReal_wp tau;
    GetQR(VDM, tau);

    projF.Reallocate(this->all_points.GetM());
    for (int i = 0; i < this->all_points.GetM(); i++)
      projF(i) = f(this->all_points(i));
    
    SolveQR(VDM, tau, projF);
  }


  //! Computes basis functions on a given point
  template<class T, class Dimension>  
  void RectangularBSplineIndex<T, Dimension>
  ::ComputeValuesPhiRef(int i, int j, const R_N& point, 
			const R_N& point_loc, Vector<Real_wp, VectSparse>& phi)
  {
    Real_wp x = point(0), y = point(1);
    int pos_x = spline_x.EvaluateFctBasis(x, value_phi_x, false);
    int pos_y = spline_y.EvaluateFctBasis(y, value_phi_y, false);
    
    int nb_val = (spline_x.GetOrder()+1)*(spline_y.GetOrder()+1);
    phi.Reallocate(nb_val);
    int nb = 0;
    for (int m = pos_x-spline_x.GetOrder(); m <= pos_x; m++)
      for (int n = pos_y-spline_y.GetOrder(); n <= pos_y; n++)
	{
	  int num = m*spline_y.GetNbDof() + n;
	  phi.Index(nb) = num;
	  phi.Value(nb) = value_phi_x(m)*value_phi_y(n);
	  nb++;
	}
  }
  
}

#define MONTJOIE_FILE_MULTI_DIM_PHYSICAL_INDEX_CXX
#endif

