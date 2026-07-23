#ifndef MONTJOIE_FILE_LOCAL_DISCONTINUOUS_GALERKIN_1D_CXX

namespace Montjoie
{
  
  template<class TypeEquation>
  LdgScheme1D<TypeEquation>::LdgScheme1D()
  {
    //mass_lumping = false;
    periodic_condition = true;
    mass_lumping = true;
  }
  
  
  //! modification of parameters with a line of the data file
  template<class TypeEquation>
  void LdgScheme1D<TypeEquation>::SetInputData(const string& keyword, const Vector<string>& param)
  {
    if (keyword == "MassLumping")
      {
	if (param.GetM() <= 0)
	  {
	    cout << "In SetInputData of LdgScheme1D" << endl;
	    cout << "MassLumping needs more parameters, for instance :" << endl;
	    cout << "MassLumping = YES" << endl;
	    cout << "Current parameters are : " << endl << param << endl;
	    abort();
	  }

        if (param(0) == "NO")
          mass_lumping = false;
        else
          mass_lumping = true;
      }
  }
  
  
  // returns quadrature number of element i and quadrature point j
  template<class TypeEquation>
  int LdgScheme1D<TypeEquation>::GetQuadNumber(int i, int j) const
  {
    return mesh.GetNbVertices() + i*gauss.GetNbPointsQuad() + j;
  }
  
  
  // sets the number of intervals of the computational domain
  template<class TypeEquation>
  void LdgScheme1D<TypeEquation>::
  SetInterval(const Real_wp& x0, const Real_wp& xN, int n, int r)
  {
    mesh.CreateRegularMesh(x0, xN, n+1, 1);
  }
  
  
  // initialisation of the interpolation grid
  template<class TypeEquation>
  void LdgScheme1D<TypeEquation>
  ::InitGrid(const Real_wp& xmin, const Real_wp& xmax, int Nd)
  {
    grid.Init(xmin, xmax, Nd);
    grid.LocalizePoints(mesh);
  }
  
  
  // modifies the order of approximation
  template<class TypeEquation>
  void LdgScheme1D<TypeEquation>::SetOrder(int order)
  {
    // we number the mesh
    mesh.SetOrder(order);
    if (!TypeEquation::FormulationDG)
      if (periodic_condition)
        mesh.SetPeriodicExtremity();
    
    mesh.NumberMesh(TypeEquation::FormulationDG);
    
    // we compute h_i and invH
    int nb_elt = mesh.GetNbElt();
    h_subdiv.Reallocate(nb_elt);
    invH_subdiv.Reallocate(nb_elt);
    for (int i = 0; i < nb_elt; i++)
      {
	h_subdiv(i) = mesh.Vertex(i+1) - mesh.Vertex(i);
	invH_subdiv(i) = 1.0/h_subdiv(i);
      }
    
    // computation of interpolation and quadrature points
    lob.ConstructQuadrature(order, lob.QUADRATURE_LOBATTO);
    gauss.ConstructQuadrature(order, lob.QUADRATURE_GAUSS);
    // for mass lumping, computation of invDh
    invDh.Reallocate(order+1);
    Dh.Reallocate(order+1);
    for (int i = 0; i <= order; i++)
      {
	Dh(i) = lob.Weights(i);
	invDh(i) = 1.0/lob.Weights(i);
      }
    
    // computation of phi_i(\xi_j) and dphi_i(\xi_j)
    int N = gauss.GetNbPointsQuad();
    ValPhi.Reallocate(order+1, N);
    GradPhi.Reallocate(order+1, N);
    ValPhiWeight.Reallocate(order+1, N);
    GradPhiWeight.Reallocate(order+1, N);
    G_GL.Reallocate(N, order+1);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j < N; j++)
	{
	  ValPhi(i, j) = lob.EvaluatePhi(i, gauss.Points(j));
	  GradPhi(i, j) = lob.EvaluatePhiGrad(i, gauss.Points(j));
	  ValPhiWeight(i, j) = ValPhi(i, j)*gauss.Weights(j);
	  GradPhiWeight(i, j) = GradPhi(i, j)*gauss.Weights(j);
	  G_GL(j, i) = gauss.EvaluatePhi(j, lob.Points(i));
	}
        
    // computation of mass matrix without mass lumping
    mat_mass.Reallocate(order+1, order+1); mat_mass.Fill(0);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k < N; k++)
	  mat_mass(i, j) += gauss.Weights(k)*ValPhi(i, k)*ValPhi(j, k);
    
    invMass = mat_mass; GetInverse(invMass);
    mat_rigid.Reallocate(order+1, order+1); mat_rigid.Fill(0);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k < N; k++)
	  mat_rigid(i, j) += gauss.Weights(k)*ValPhi(i, k)*GradPhi(j, k);

    // global mass matrix with mass lumping
    if (!TypeEquation::FormulationDG)
      {
	invGlobalMass.Reallocate(mesh.GetNbDof());
	invGlobalMass.Fill(0);
	for (int i = 0; i < nb_elt; i++)
	  for (int j = 0; j <= order; j++)
	    invGlobalMass(mesh.GetNumberDof(i, j)) += h_subdiv(i)*lob.Weights(j);
	
	for (int i = 0; i < mesh.GetNbDof(); i++)
	  invGlobalMass(i) = 1.0 / invGlobalMass(i);
      }
    
    // computing position of degrees of freedom
    Xdof.Reallocate(mesh.GetNbDof());
    for (int i = 0; i < nb_elt; i++)
      for (int j = 0; j <= order; j++)
	Xdof(mesh.GetNumberDof(i, j))
          = mesh.Vertex(i) + (mesh.Vertex(i+1)-mesh.Vertex(i))*lob.Points(j);
    
    var_eq.ComputePhysicalProperty(*this, mesh.GetXmin(), mesh.GetXmax());
  }
  
  
  // solving by M^{-1} on element i
  template<class TypeEquation> template<class T>
  inline void LdgScheme1D<TypeEquation>::SolveMass(int i, Vector<T>& Zn) const
  {
    int order = lob.GetOrder();
    if (mass_lumping)
      {
	for (int j = 0; j <= order; j++)
	  Zn(j) *= invH_subdiv(i)*invDh(j);
      }
    else
      {
	Vector<T> ZnLoc = Zn;
	for (int j = 0; j <= order; j++)
	  {
	    T vloc = 0;
	    for (int k = 0; k <= order; k++)
	      vloc += invMass(j, k)*ZnLoc(k);
	    
	    Zn(j) = vloc*invH_subdiv(i);
	  }
      }
  }
  
  
  // solving by M^{-1} on element i
  template<class TypeEquation> template<class T0, class T1>
  inline void LdgScheme1D<TypeEquation>::
  MltMass(int i, const T0& alpha, const Vector<T1>& Un, Vector<T1>& Zn) const
  {
    int order = lob.GetOrder();
    alpha *= h_subdiv(i);
    if (mass_lumping)
      {
	for (int j = 0; j <= order; j++)
	  Zn(j) += alpha*Dh(j)*Un(j);
      }
    else
      {
	for (int j = 0; j <= order; j++)
	  {
	    T1 vloc = 0;
	    for (int k = 0; k <= order; k++)
	      vloc += mat_mass(j, k)*Un(k);
	    
	    Zn(j) += vloc*alpha;
	  }
      }
  }
  
  
  // computation of Zn_j = Zn(xi_j) on element i
  template<class TypeEquation> template<class T>
  inline void LdgScheme1D<TypeEquation>
  ::GetFunction(int i, const Vector<T>& Zn, Vector<T>& ZnQuad) const
  {
    T vloc; int N = gauss.GetNbPointsQuad();
    int order = lob.GetOrder();
    for (int j = 0; j < N; j++)
      {
	vloc = 0;
	for (int k = 0; k <= order; k++)
	  vloc += ValPhi(k, j)*Zn(k);
	
	ZnQuad(j) = vloc;
      }
  }
    
  
  // adding integral \int Zn phi_j  on element i
  template<class TypeEquation> template<class T0, class T1>
  inline void LdgScheme1D<TypeEquation>::
  IntegrateAgainstFunction(int i, const T0& alpha_, const Vector<T1>& ZnQuad, Vector<T1>& Zn) const
  {
    T1 vloc;
    int order = lob.GetOrder();
    int N = gauss.GetNbPointsQuad();
    T0 alpha = alpha_*h_subdiv(i);
    for (int j = 0; j <= order; j++)
      {
	vloc = 0;
	for (int k = 0; k < N; k++)
	  vloc += ValPhiWeight(j, k)*ZnQuad(k);
	
	Zn(j) += vloc*alpha;
      }
  }
  
  
  // computation of dZn_i = dZn/dx on quadrature points
  template<class TypeEquation> template<class T>
  inline void LdgScheme1D<TypeEquation>
  ::GetDerivative(int i, const Vector<T>& Zn, Vector<T>& dZn) const
  {
    T vloc;
    int N = gauss.GetNbPointsQuad();
    int order = lob.GetOrder();
    for (int k = 0; k < N; k++)
      {
	vloc = 0;
	for (int j = 0; j <= order; j++)
	  vloc += GradPhi(j, k)*Zn(j);
	
	dZn(k) = vloc*invH_subdiv(i);
      }
  }
  
  
  // adding integral \int Zn dphi_j/dx  on element i
  template<class TypeEquation> template<class T0, class T1>
  inline void LdgScheme1D<TypeEquation>::
  IntegrateAgainstDerivative(int i, const T0& alpha, const Vector<T1>& dZn, Vector<T1>& Zn) const
  {
    int N = gauss.GetNbPointsQuad();
    int order = lob.GetOrder();
    for (int k = 0; k <= order; k++)
      {
	T1 vloc = 0;
	for (int j = 0; j < N; j++)
	  vloc += GradPhiWeight(k, j)*dZn(j);
	
	Zn(k) += alpha*vloc;
      }
  }
  
  
  // computation of values of Zn on Gauss-Lobatto points from values on quadrature points
  template<class TypeEquation> template<class T>
  inline void LdgScheme1D<TypeEquation>::ProjectToDof(const Vector<T>& ZnQuad, Vector<T>& Zn) const
  {
    int N = gauss.GetNbPointsQuad();
    int order = lob.GetOrder();
    Zn.Fill(0);
    for (int j = 0; j < N; j++)
      for (int k = 0; k <= order; k++)
	Zn(k) += G_GL(j, k)*ZnQuad(j);
  }
  
  
  // retrieves u- and u+ around vertex i
  template<class TypeEquation> template<class T>
  inline void LdgScheme1D<TypeEquation>::GetUplusMinus(int i, int offset, const Vector<T>& Un,
                                                       const IVect& new_num, T& Um, T& Up) const
  {
    if ((i == 0)||(i == mesh.GetNbElt()))
      {
	Um = Un(new_num(offset+mesh.GetNbDof()-1));
	Up = Un(new_num(offset));
      }
    else
      {
	int j = offset + i*(mesh.GetOrder()+1);
	Um = Un(new_num(j-1));
	Up = Un(new_num(j));
      }
  }
  
  
  // retrieves u- and u+ around vertex i
  template<class TypeEquation> template<int m, class T>
  inline void LdgScheme1D<TypeEquation>
  ::GetUplusMinus(int i, int offset_un, const Vector<T>& Un, const IVect& new_num,
                  TinyVector<T, m>& Um, TinyVector<T, m>& Up) const
  {
    if ((i == 0)||(i == mesh.GetNbElt()))
      {
	int offset = offset_un;
	for (int n = 0; n < m; n++)
	  {
	    Um(n) = Un(new_num(offset+mesh.GetNbDof()-1));
	    Up(n) = Un(new_num(offset));
	    offset += mesh.GetNbDof();
	  }
      }
    else
      {
	int offset = offset_un + i*(mesh.GetOrder()+1);
	for (int n = 0; n < m; n++)
	  {
	    Um(n) = Un(new_num(offset-1));
	    Up(n) = Un(new_num(offset));
	    offset += mesh.GetNbDof();
	  }
      }
  }
  
  
  // we are computing Zn = Zn + alpha*(\int f(u) u_x phi dx 
  // + \int g(u) phi_x + \int h(u) phi dx + flux_plus phi(1) - flux_minus phi(0))
  template<class TypeEquation> template<class OperatorEq, class T0, class T1>
  void LdgScheme1D<TypeEquation>
  ::MltAddStiffness(const T0& alpha, const Vector<T1>& Un,
                    int offset_zn, const OperatorEq& var,
                    const IVect& new_num, Vector<T1>& Zn, bool solve_mass) const
  {
    // loop on elements
    int Nvol = mesh.GetNbDof();
    int order = mesh.GetOrder();
    int N = gauss.GetNbPointsQuad();
    TinyVector<Vector<T1>, TypeEquation::nb_unknowns>
      UnLoc, UnQuad, dUnQuad, feval, feval_grad, ZnLoc;
    
    TinyVector<T1, TypeEquation::nb_unknowns> u, du, Um, Up, fm, fp, flux, ftilde;
    
    // allocating arrays
    for (int m = 0; m < TypeEquation::nb_unknowns; m++)
      {
	UnLoc(m).Reallocate(order+1); UnLoc(m).Fill(0);
	ZnLoc(m).Reallocate(order+1); ZnLoc(m).Fill(0);
	UnQuad(m).Reallocate(N); UnQuad(m).Fill(0);
	dUnQuad(m).Reallocate(N); dUnQuad(m).Fill(0);
	feval(m).Reallocate(N); feval(m).Fill(0);
        feval_grad(m).Reallocate(N); feval_grad(m).Fill(0);
      }
    
    TinyVector<T1, TypeEquation::nb_unknowns> g, h;
    g.Zero(); h.Zero();
    
    // loop over elements
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
	// we retrieve value of u on element
	for (int m = 0; m < TypeEquation::nb_unknowns; m++)
	  for (int j = 0; j <= order; j++)
	    UnLoc(m)(j) = Un(new_num(offset_zn + m*Nvol + mesh.GetNumberDof(i, j)));
		
	// we compute u, du/dx on quadrature points
        for (int m = 0; m < TypeEquation::nb_unknowns; m++)
          {
            GetFunction(i, UnLoc(m), UnQuad(m));
            GetDerivative(i, UnLoc(m), dUnQuad(m));
          }
        
	// computation of f(u), g(u) and h(u) on quadrature points 
	for (int j = 0; j < N; j++)
	  {
	    // we get u, du
	    CopyVector(UnQuad, j, u); 
	    CopyVector(dUnQuad, j, du); 
	    	    	    
	    var.EvaluateFunction(GetQuadNumber(i, j), u, du, g, h);
            	    
	    CopyVector(h, j, feval);
            CopyVector(g, j, feval_grad);
	  }
        
        FillZero(ZnLoc);
        
	// integration against basis functions
	for (int m = 0; m < TypeEquation::nb_unknowns; m++)
	  IntegrateAgainstFunction(i, 1.0, feval(m), ZnLoc(m));
	
        // integration against gradient of basis functions
        for (int m = 0; m < TypeEquation::nb_unknowns; m++)
          IntegrateAgainstDerivative(i, 1.0, feval_grad(m), ZnLoc(m));
        
	// now we evaluate numerical flux
        if (TypeEquation::FormulationDG)
          {
            GetUplusMinus(i, offset_zn, Un, new_num, Um, Up);
	    
            var.EvaluateNumericalFlux(i, 0, Um, Up, flux);
            
            for (int m = 0; m < TypeEquation::nb_unknowns; m++)
              ZnLoc(m)(0) += flux(m);
            
            GetUplusMinus(i+1, offset_zn, Un, new_num, Um, Up);
	    
            var.EvaluateNumericalFlux(i+1, 1, Um, Up, flux);
            
            for (int m = 0; m < TypeEquation::nb_unknowns; m++)
              ZnLoc(m)(order) += flux(m);
          }
	
	// inversion by mass matrix if required
	if (TypeEquation::FormulationDG)
	  {
	    if (solve_mass)
	      for (int m = 0; m < TypeEquation::nb_unknowns; m++)
		SolveMass(i, ZnLoc(m));
	  }
	
	for (int m = 0; m < TypeEquation::nb_unknowns; m++)
	  {
	    for (int j = 0; j <= order; j++)
	      Zn(new_num(offset_zn + m*Nvol + mesh.GetNumberDof(i, j))) += alpha*ZnLoc(m)(j);
	  }
	
      }
    
    if (!TypeEquation::FormulationDG)
      for (int m = 0; m < TypeEquation::nb_unknowns; m++)
	for (int j = 0; j < Nvol; j++)
	  Zn(new_num(offset_zn + m*Nvol + j)) *= invGlobalMass(j);
  }
    
  
  // we are computing Zn = Zn + alpha*(\int f(u) du/dt phi_x + \int h(u) du/dt phi dx 
  // + flux_plus du/dt phi(1) - flux_minus du/dt phi(0))
  template<class TypeEquation> template<class OperatorEq, class T0, class T1>
  void LdgScheme1D<TypeEquation>
  ::MltAddMass(const T0& alpha, const Vector<T1>& Un, const Vector<T1>& dU_dt,
               int offset_zn, const OperatorEq& var, const IVect& new_num, Vector<T1>& Zn) const
  {
    // loop on elements
    int Nvol = mesh.GetNbDof();
    int order = mesh.GetOrder();
    int N = gauss.GetNbPointsQuad();
    TinyVector<Vector<T1>, TypeEquation::nb_unknowns> UnLoc, UnDtLoc, dUnQuad, UnQuad,
      UnDtQuad, dUnDtQuad, feval, feval_grad, ZnLoc;
    TinyVector<T1, TypeEquation::nb_unknowns> u, du, du_dt, du_dxt, Um, Up, dUm, dUp;
    
    // allocating arrays
    for (int m = 0; m < TypeEquation::nb_unknowns; m++)
      {
	UnLoc(m).Reallocate(order+1); UnLoc(m).Fill(0);
	UnDtLoc(m).Reallocate(order+1); UnDtLoc(m).Fill(0);
	ZnLoc(m).Reallocate(order+1); ZnLoc(m).Fill(0);
	UnQuad(m).Reallocate(N); UnQuad(m).Fill(0);
        dUnQuad(m).Reallocate(N); dUnQuad(m).Fill(0);
	UnDtQuad(m).Reallocate(N); UnDtQuad(m).Fill(0);
        dUnDtQuad(m).Reallocate(N); dUnDtQuad(m).Fill(0);
	feval(m).Reallocate(N); feval(m).Fill(0);
        feval_grad(m).Reallocate(N); feval_grad(m).Fill(0);
      }
    
    TinyMatrix<T1, General, TypeEquation::nb_unknowns, TypeEquation::nb_unknowns>
      f, h, fx, hx, flux_Um, flux_Up;
    
    TinyVector<T1, TypeEquation::nb_unknowns> h2, g2;
    
    f.Zero(); h.Zero(); fx.Zero(); hx.Zero(); h2.Zero(); g2.Zero();
    
    // loop over elements
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
	// we retrieve value of u on element
	for (int m = 0; m < TypeEquation::nb_unknowns; m++)
	  for (int j = 0; j <= order; j++)
	    {
              UnLoc(m)(j) = Un(new_num(offset_zn + m*Nvol + mesh.GetNumberDof(i, j)));
              UnDtLoc(m)(j) = dU_dt(new_num(offset_zn + m*Nvol + mesh.GetNumberDof(i, j)));
            }
        
        for (int m = 0; m < TypeEquation::nb_unknowns; m++)
          {
            GetFunction(i, UnLoc(m), UnQuad(m));
            GetDerivative(i, UnLoc(m), dUnQuad(m));
            GetFunction(i, UnDtLoc(m), UnDtQuad(m));
            GetDerivative(i, UnDtLoc(m), dUnDtQuad(m));
          }
        
	// computation of f(u), g(u) and h(u) on quadrature points 
	for (int j = 0; j < N; j++)
	  {
	    // we get u, du
	    CopyVector(UnQuad, j, u); 
            CopyVector(dUnQuad, j, du); 
            CopyVector(UnDtQuad, j, du_dt); 
            CopyVector(dUnDtQuad, j, du_dxt); 
            
	    var.EvaluateFunction(GetQuadNumber(i, j), u, du, f, h, fx, hx);
            h2 = dot(h, du_dt) + dot(hx, du_dxt);
            g2 = dot(f, du_dt) + dot(fx, du_dxt);
            	    
            CopyVector(h2, j, feval);
            CopyVector(g2, j, feval_grad);
	  }
        
        FillZero(ZnLoc);
        
	// integration against basis functions
	for (int m = 0; m < TypeEquation::nb_unknowns; m++)
	  IntegrateAgainstFunction(i, 1.0, feval(m), ZnLoc(m));
	
        // integration against gradient of basis functions
        for (int m = 0; m < TypeEquation::nb_unknowns; m++)
          IntegrateAgainstDerivative(i, 1.0, feval_grad(m), ZnLoc(m));
        
	// now we evaluate numerical flux
        if (TypeEquation::FormulationDG)
          for (int k = 0; k < 2; k++)
            {
              GetUplusMinus(i+k, offset_zn, Un, new_num, Um, Up);
              GetUplusMinus(i+k, offset_zn, dU_dt, new_num, dUm, dUp);
              
              var.EvaluateNumericalFlux(i+k, k, Um, Up, flux_Um, flux_Up);
              g2 = dot(flux_Um, dUm);
              h2 = dot(flux_Up, dUp);
              
              for (int m = 0; m < TypeEquation::nb_unknowns; m++)
                ZnLoc(m)(k*order) += g2(m) + h2(m);
            }
        
	for (int m = 0; m < TypeEquation::nb_unknowns; m++)
	  {
	    for (int j = 0; j <= order; j++)
	      Zn(new_num(offset_zn + m*Nvol + mesh.GetNumberDof(i, j))) += alpha*ZnLoc(m)(j);
	  }
      }
  }
  
  
  //! evaluation of Fn = Fn + F(U)
  template<class TypeEquation> template<class T0, class T1>
  void LdgScheme1D<TypeEquation>
  ::AddSystem(const Vector<T1>& Un_mass, const Vector<T1>& Un, int offset_zn,
              const T0& alpha_mass, const T0& beta_stiff,
              const IVect& new_num, Vector<T1>& Fn) const
  {
    // loop over elements
    if (alpha_mass != T0(0))
      MltAddMass(alpha_mass, Un, Un_mass, offset_zn, var_eq.var_mass, new_num, Fn);
	
    //if (beta_stiff != Real_wp(0))
    MltAddStiffness(beta_stiff, Un, offset_zn, var_eq.var_stiff, new_num, Fn);
  }
  
  
  // computation of differential for Newton algorithm
  template<class TypeEquation> template<class OperatorEq, class MatrixSparse, class T0, class T1>
  void LdgScheme1D<TypeEquation>::
  AddStiffnessDifferential(const Vector<T0>& alpha, const Vector<T1>& Un,
                           int offset_zn, const OperatorEq& var,
                           const IVect& new_num, MatrixSparse& DF) const
  {
    int Nvol = mesh.GetNbDof();
    int nodl = Nvol*TypeEquation::nb_unknowns;
    int order = mesh.GetOrder();
    int N = gauss.GetNbPointsQuad();
    TinyVector<Vector<T1>, TypeEquation::nb_unknowns> UnLoc, UnQuad, dUnQuad, feval, ZnLoc;
    TinyVector<T1, TypeEquation::nb_unknowns> u, du, Um, Up, fm, fp, flux;

    for (int m = 0; m < TypeEquation::nb_unknowns; m++)
      {
	UnLoc(m).Reallocate(order+1); UnLoc(m).Fill(0);
	ZnLoc(m).Reallocate(order+1); ZnLoc(m).Fill(0);
	UnQuad(m).Reallocate(N); UnQuad(m).Fill(0);
	dUnQuad(m).Reallocate(N); dUnQuad(m).Fill(0);
	feval(m).Reallocate(N); feval(m).Fill(0);
      }
        
    TinyMatrix<T1, General, TypeEquation::nb_unknowns, TypeEquation::nb_unknowns>
      dh_du, dg_du, dg_ddu, dh_ddu;
    
    TinyMatrix<T1, General, TypeEquation::nb_unknowns, TypeEquation::nb_unknowns>
      dflux_dum, dflux_dup;
    TinyVector<T1, TypeEquation::nb_unknowns> g, h;
    
    TinyMatrix<Matrix<T1>, General, TypeEquation::nb_unknowns, TypeEquation::nb_unknowns> dF_uu;
    for (int m = 0; m < TypeEquation::nb_unknowns; m++)
      for (int n = 0; n < TypeEquation::nb_unknowns; n++)
	{
	  dF_uu(m, n).Reallocate(order+1, order+1);
	  dF_uu(m, n).Fill(0);
	}
    
    int c_interac; IVect col_interac(TypeEquation::nb_unknowns*(order+1)*alpha.GetM());
    col_interac.Zero();
    Vector<T1> val_interac(TypeEquation::nb_unknowns*(order+1)*alpha.GetM()); val_interac.Fill(0);
    
    int nb_m = 0, nb_p = 0; int nb = 0;
    // main loop over elements
    for (int ne = 0; ne < mesh.GetNbElt(); ne++)
      {
	if (ne == 0)
	  nb_m = Nvol-1;
	else
	  nb_m = nb-1;
	
	if (ne == mesh.GetNbElt()-1)
	  nb_p = 0;
	else
	  nb_p = nb+order+1;
	
        Real_wp he = h_subdiv(ne);
	Real_wp invHe = 1.0/he;
	FillZero(dF_uu);
	
	// we retrieve values of u on element
	for (int m = 0; m < TypeEquation::nb_unknowns; m++)
	  for (int j = 0; j <= order; j++)
	    UnLoc(m)(j) = Un(new_num(offset_zn + m*Nvol + mesh.GetNumberDof(ne, j)));
        
	// we compute u, du/dx
        for (int m = 0; m < TypeEquation::nb_unknowns; m++)
          {
            GetFunction(ne, UnLoc(m), UnQuad(m));
            GetDerivative(ne, UnLoc(m), dUnQuad(m));
          }	
	
	// computation of f(u, p), 
	FillZero(dF_uu); T1 vloc = 0;
        
        // volumic part
        for (int q = 0; q < N; q++)
          {
            // we get u, du, p, dp
            CopyVector(UnQuad, q, u); 
            CopyVector(dUnQuad, q, du); 
            
            var.EvaluateDerivative(GetQuadNumber(ne, q), u, du, g, h,
                                   dg_du, dh_du, dg_ddu, dh_ddu);
            
            for (int i = 0; i <= order; i++)
              for (int k = 0; k <= order; k++)
                for (int m = 0; m < TypeEquation::nb_unknowns; m++)
                  {
                    for (int p = 0; p < TypeEquation::nb_unknowns; p++)
                      {
                        vloc = ValPhiWeight(k, q)*dh_du(m, p)*ValPhi(i, q);
                        vloc += invHe*GradPhiWeight(k, q)*dh_ddu(m, p)*ValPhi(i, q);
                        vloc += invHe*ValPhiWeight(k, q)*dg_du(m, p)*GradPhi(i, q);
                        vloc += invHe*invHe*GradPhiWeight(k, q)*dg_ddu(m, p)*GradPhi(i, q);
                                                
                        dF_uu(m, p)(i, k) += vloc*he;
                      }
                  }
          }
        
        if (TypeEquation::FormulationDG)
          {
            // now we differentiate numerical flux	
            GetUplusMinus(ne, offset_zn, Un, new_num, Um, Up);
            
            var.EvaluateDerivativeNumericalFlux(ne, 0, Um, Up, flux, dflux_dum, dflux_dup);
            
            for (int m = 0; m < TypeEquation::nb_unknowns; m++)
              {
                for (int p = 0; p < TypeEquation::nb_unknowns; p++)
                  {
                    dF_uu(m, p)(0,0) += dflux_dup(m, p);
                    c_interac = alpha.GetM();
                    for (int nc = 0; nc < alpha.GetM(); nc++)
                      {
                        col_interac(nc) = new_num(nb_m + p*Nvol + nc*nodl);
                        val_interac(nc) = alpha(nc)*dflux_dum(m, p);
                      }
                    
                    DF.AddInteractionRow(new_num(offset_zn + nb + m*Nvol),
                                         c_interac, col_interac, val_interac);
                  }
              }
            
            
            // extremity r
            GetUplusMinus(ne+1, offset_zn, Un, new_num, Um, Up);
            
            var.EvaluateDerivativeNumericalFlux(ne+1, 1, Um, Up, flux, dflux_dum, dflux_dup);
            
            for (int m = 0; m < TypeEquation::nb_unknowns; m++)
              {
                for (int p = 0; p < TypeEquation::nb_unknowns; p++)
                  {
                    dF_uu(m, p)(order,order) += dflux_dum(m, p);
                    c_interac = alpha.GetM();
                    for (int nc = 0; nc < alpha.GetM(); nc++)
                      {
                        col_interac(nc) = new_num(nb_p + p*Nvol + nc*nodl);
                        val_interac(nc) = alpha(nc)*dflux_dup(m, p);
                      }
                    
                    DF.AddInteractionRow(new_num(offset_zn + nb + order + m*Nvol),
                                         c_interac, col_interac, val_interac);
                  }
              }
          }
        
        // we add contributions to sparse matrix
        for (int m = 0; m < TypeEquation::nb_unknowns; m++)
          for (int i = 0; i <= order; i++)
            {
              c_interac = 0;
              int irow = new_num(m*Nvol + mesh.GetNumberDof(ne, i) + offset_zn);
              for (int p = 0; p < TypeEquation::nb_unknowns; p++)
                for (int k = 0; k <= order; k++)
                  {
		    if (abs(dF_uu(m, p)(i, k)) > epsilon_machine)
		      {
			for (int nc = 0; nc < alpha.GetM(); nc++)
			  {
			    col_interac(c_interac) = new_num(mesh.GetNumberDof(ne, k)
                                                             + p*Nvol + nc*nodl);
			    val_interac(c_interac) = alpha(nc)*dF_uu(m, p)(i, k);
			    c_interac++;
			  }
		      }
		  }
	      
	      DF.AddInteractionRow(irow, c_interac, col_interac, val_interac);
              
	    }	
        
        if (TypeEquation::FormulationDG)
          nb += order+1;
        else
          nb += order;
      }
  }
  
  
  // computation of differential for Newton algorithm
  template<class TypeEquation> template<class OperatorEq, class MatrixSparse, class T0, class T1>
  void LdgScheme1D<TypeEquation>
  ::AddMassDifferential(const Vector<T0>& alpha, const Vector<T1>& Un,
                        const Vector<T1>& dU_dt, int offset_zn,
                        const OperatorEq& var, const IVect& new_num, MatrixSparse& DF) const
  {
    int Nvol = mesh.GetNbDof(); int nb = 0;
    int nodl = Nvol*TypeEquation::nb_unknowns;
    int order = mesh.GetOrder();
    int N = gauss.GetNbPointsQuad();
    TinyVector<Vector<T1>, TypeEquation::nb_unknowns> UnLoc, UnDtLoc,
      UnQuad, dUnQuad, UnDtQuad, dUnDtQuad, feval, ZnLoc;
    TinyVector<T1, TypeEquation::nb_unknowns> u, du, du_dt, du_dxt, Um, Up, dUm, dUp, fm, fp, flux;

    for (int m = 0; m < TypeEquation::nb_unknowns; m++)
      {
	UnLoc(m).Reallocate(order+1); UnLoc(m).Fill(0);
        UnDtLoc(m).Reallocate(order+1); UnDtLoc(m).Fill(0);
	ZnLoc(m).Reallocate(order+1); ZnLoc(m).Fill(0);
	UnQuad(m).Reallocate(N); UnQuad(m).Fill(0);
        dUnQuad(m).Reallocate(N); dUnQuad(m).Fill(0);
	UnDtQuad(m).Reallocate(N); UnDtQuad(m).Fill(0);
        dUnDtQuad(m).Reallocate(N); dUnDtQuad(m).Fill(0);
	feval(m).Reallocate(N); feval(m).Fill(0);
      }
    
    TinyArray3D<T1, TypeEquation::nb_unknowns,
      TypeEquation::nb_unknowns, TypeEquation::nb_unknowns> df_du, dh_du,
      dfluxUm_dup, dfluxUm_dum, dfluxUp_dum, dfluxUp_dup, dh_ddu, df_ddu;

    TinyArray3D<T1, TypeEquation::nb_unknowns,
      TypeEquation::nb_unknowns, TypeEquation::nb_unknowns> dfx_du, dhx_du, dhx_ddu, dfx_ddu;
        
    TinyMatrix<T1, General, TypeEquation::nb_unknowns, TypeEquation::nb_unknowns>
      f, h, fx, hx, flux_Um, flux_Up;
    
    TinyMatrix<Matrix<T1>, General, TypeEquation::nb_unknowns, TypeEquation::nb_unknowns>
      dF_uu, dF_du;
    
    for (int m = 0; m < TypeEquation::nb_unknowns; m++)
      for (int n = 0; n < TypeEquation::nb_unknowns; n++)
	{
	  dF_uu(m, n).Reallocate(order+1, order+1);
	  dF_uu(m, n).Fill(0);
	  dF_du(m, n).Reallocate(order+1, order+1);
	  dF_du(m, n).Fill(0);
	}

    int c_interac; IVect col_interac(TypeEquation::nb_unknowns*(order+1)*(alpha.GetM()+1));
    col_interac.Zero();
    Vector<T1> val_interac(TypeEquation::nb_unknowns*(order+1)*(alpha.GetM()+1));
    val_interac.Fill(0);
    
    int nb_m = 0, nb_p = 0; 
    // main loop over elements
    for (int ne = 0; ne < mesh.GetNbElt(); ne++)
      {
	if (ne == 0)
	  nb_m = Nvol-1;
	else
	  nb_m = nb-1;
	
	if (ne == mesh.GetNbElt()-1)
	  nb_p = 0;
	else
	  nb_p = nb+order+1;
	
        Real_wp he = h_subdiv(ne);
	Real_wp invHe = 1.0/he;
	FillZero(dF_uu);
	
	// we retrieve values of u on element
	for (int m = 0; m < TypeEquation::nb_unknowns; m++)
	  for (int j = 0; j <= order; j++)
	    {
              UnLoc(m)(j) = Un(new_num(offset_zn + m*Nvol + mesh.GetNumberDof(ne, j)));
              UnDtLoc(m)(j) = dU_dt(new_num(offset_zn + m*Nvol + mesh.GetNumberDof(ne, j)));
            }
        
	// we compute u, du/dx
        for (int m = 0; m < TypeEquation::nb_unknowns; m++)
          {
            GetFunction(ne, UnLoc(m), UnQuad(m));
            GetDerivative(ne, UnLoc(m), dUnQuad(m));
            GetFunction(ne, UnDtLoc(m), UnDtQuad(m));
            GetDerivative(ne, UnDtLoc(m), dUnDtQuad(m));
          }	
	
	// computation of f(u, p), 
	FillZero(dF_uu); FillZero(dF_du);
        T1 vloc = 0;
        
        // volumic part
        for (int q = 0; q < N; q++)
          {
            // we get u, du, p, dp
            CopyVector(UnQuad, q, u); 
            CopyVector(dUnQuad, q, du); 
            CopyVector(UnDtQuad, q, du_dt); 
            CopyVector(dUnDtQuad, q, du_dxt); 
            
            var.EvaluateDerivative(GetQuadNumber(ne, q), u, du, f, h, fx, hx,
                                   df_du, dh_du, df_ddu, dh_ddu, dfx_du, dhx_du, dfx_ddu, dhx_ddu);
            
            for (int i = 0; i <= order; i++)
              for (int k = 0; k <= order; k++)
                for (int m = 0; m < TypeEquation::nb_unknowns; m++)
                  {
                    for (int p = 0; p < TypeEquation::nb_unknowns; p++)
                      {
                        vloc = invHe*GradPhiWeight(i, q)*f(m, p)*ValPhi(k, q);
                        vloc += ValPhiWeight(i, q)*h(m, p)*ValPhi(k, q);
                        vloc += invHe*invHe*GradPhiWeight(i, q)*fx(m, p)*GradPhi(k, q);
                        vloc += invHe*ValPhiWeight(i, q)*hx(m, p)*GradPhi(k, q);
                        
                        dF_du(m, p)(i, k) += vloc*he;
                        
                        vloc = 0;
                        for (int n = 0; n < TypeEquation::nb_unknowns; n++)
                          {
                            vloc += invHe*ValPhiWeight(k, q)*df_du(m, n, p)*du_dt(n)*GradPhi(i, q);
                            vloc += ValPhiWeight(k, q)*dh_du(m, n, p)*du_dt(n)*ValPhi(i, q);
                            vloc += invHe*invHe*GradPhiWeight(k, q)
                              *df_ddu(m, n, p)*du_dt(n)*GradPhi(i, q);
                            
                            vloc += invHe*GradPhiWeight(k, q)*dh_ddu(m, n, p)
                              *du_dt(n)*ValPhi(i, q);
                            
                            vloc += invHe*ValPhiWeight(k, q)*dfx_du(m, n, p)*du_dxt(n)
                              *GradPhi(i, q);
                            
                            vloc += ValPhiWeight(k, q)*dhx_du(m, n, p)*du_dxt(n)*ValPhi(i, q);
                            vloc += invHe*invHe*GradPhiWeight(k, q)*dfx_ddu(m, n, p)
                              *du_dxt(n)*GradPhi(i, q);
                            
                            vloc += invHe*GradPhiWeight(k, q)*dhx_ddu(m, n, p)
                              *du_dxt(n)*ValPhi(i, q);
                          }
                        
                        dF_uu(m, p)(i, k) += vloc*he;
                      }
                  }
          }
        
        // now we differentiate numerical flux	
        if (TypeEquation::FormulationDG)
          for (int ext = 0; ext < 2; ext++)
            {
              GetUplusMinus(ne+ext, offset_zn, Un, new_num, Um, Up);
              GetUplusMinus(ne+ext, offset_zn, dU_dt, new_num, dUm, dUp);
              
              var.EvaluateDerivativeNumericalFlux(ne+ext, ext, Um, Up, flux_Um, flux_Up,
                                                  dfluxUm_dum, dfluxUm_dup,
                                                  dfluxUp_dum, dfluxUp_dup);
              
              for (int m = 0; m < TypeEquation::nb_unknowns; m++)
                {
                  for (int p = 0; p < TypeEquation::nb_unknowns; p++)
                    {
                      if (ext == 0)
                        {
                          vloc = 0;
                          for (int n = 0; n < TypeEquation::nb_unknowns; n++)
                            {
                              vloc += dfluxUp_dup(m,n,p)*dUp(n);
                              vloc += dfluxUm_dup(m,n,p)*dUm(n);
                            }
                          dF_uu(m, p)(0, 0) += vloc;
                          
                          dF_du(m, p)(0, 0) += flux_Up(m, p);
                          
                          c_interac = 0;
                          col_interac(c_interac) = new_num(offset_zn + nb_m + p*Nvol);
                          val_interac(c_interac) = flux_Um(m, p);
                          c_interac++;
                          
                          for (int nc = 0; nc < alpha.GetM(); nc++)
                            {
                              vloc = 0;
                              for (int n = 0; n < TypeEquation::nb_unknowns; n++)
                                {
                                  vloc += dfluxUp_dum(m,n,p)*dUp(n);
                                  vloc += dfluxUm_dum(m,n,p)*dUm(n);
                                }
                              
                              col_interac(c_interac) = new_num(nb_m + p*Nvol + nc*nodl);
                              val_interac(c_interac) = alpha(nc)*vloc;
                              c_interac++;
                            }
                        }
                      else
                        {
                          vloc = 0;
                          for (int n = 0; n < TypeEquation::nb_unknowns; n++)
                            {
                              vloc += dfluxUp_dum(m,n,p)*dUp(n);
                              vloc += dfluxUm_dum(m,n,p)*dUm(n);
                            }
                          dF_uu(m, p)(order, order) += vloc;
                          
                          dF_du(m, p)(order, order) += flux_Um(m, p);
                          
                          c_interac = 0;
                          col_interac(c_interac) = new_num(offset_zn + nb_p + p*Nvol);
                          val_interac(c_interac) = flux_Up(m, p);
                          c_interac++;
                          
                          for (int nc = 0; nc < alpha.GetM(); nc++)
                            {
                              vloc = 0;
                              for (int n = 0; n < TypeEquation::nb_unknowns; n++)
                                {
                                  vloc += dfluxUp_dup(m,n,p)*dUp(n);
                                vloc += dfluxUm_dup(m,n,p)*dUm(n);
                                }
                              
                              col_interac(c_interac) = new_num(nb_p + p*Nvol + nc*nodl);
                              val_interac(c_interac) = alpha(nc)*vloc;
                              c_interac++;
                            }
                        }
                      
                      DF.AddInteractionRow(new_num(offset_zn + nb + order*ext + m*Nvol),
                                           c_interac, col_interac, val_interac);
                    }
                }
            }
        
	// we add contributions to sparse matrix
	for (int m = 0; m < TypeEquation::nb_unknowns; m++)
	  for (int i = 0; i <= order; i++)
	    {
	      c_interac = 0;
	      int irow = new_num(m*Nvol + mesh.GetNumberDof(ne, i) + offset_zn);
	      for (int p = 0; p < TypeEquation::nb_unknowns; p++)
		for (int k = 0; k <= order; k++)
		  {
		    if (abs(dF_uu(m, p)(i, k)) > epsilon_machine)
		      {
			for (int nc = 0; nc < alpha.GetM(); nc++)
			  {
			    col_interac(c_interac) = new_num(p*Nvol + mesh.GetNumberDof(ne, k)
                                                             + nc*nodl);
			    val_interac(c_interac) = alpha(nc)*dF_uu(m, p)(i, k);
			    c_interac++;
			  }
		      }

                    if (abs(dF_du(m, p)(i, k)) > epsilon_machine)
		      {
                        col_interac(c_interac) = new_num(offset_zn + p*Nvol
                                                         + mesh.GetNumberDof(ne, k));
                        val_interac(c_interac) = dF_du(m, p)(i, k);
                        c_interac++;
                      }
		  }
	      
	      DF.AddInteractionRow(irow, c_interac, col_interac, val_interac);
              
	    }	
	
        if (TypeEquation::FormulationDG)
          nb += order+1;
        else
          nb += order;
      }
  }
  
  
  // computation of differential for Newton algorithm
  template<class TypeEquation> template<class MatrixSparse, class T0, class T1>
  void LdgScheme1D<TypeEquation>
  ::GetDifferential(const Vector<T1>& Un_mass, const Vector<T1>& Un, int offset_zn,
                    const Vector<T0>& beta_stiff, const IVect& new_num, MatrixSparse& DF) const
  {    
    AddMassDifferential(beta_stiff, Un, Un_mass, offset_zn, var_eq.var_mass, new_num, DF);
    
    AddStiffnessDifferential(beta_stiff, Un, offset_zn, var_eq.var_stiff, new_num, DF);
    
  }
  
  
  //! returns energy of the scheme
  template<class TypeEquation>
  Real_wp LdgScheme1D<TypeEquation>::GetEnergy(const VectReal_wp& Zn, const IVect& new_num) const
  {
    VectReal_wp ProdZn(mesh.GetNbDof()*TypeEquation::nb_unknowns);
    ProdZn.Fill(0);
    MltAddMass(1.0, Zn, Zn, 0, this->var_eq.var_mass, new_num, ProdZn);
    Real_wp energy = DotProd(Zn, ProdZn);
    energy = sqrt(energy);
    return energy;
  }
  
  
  template<class TypeEquation>
  Real_wp LdgScheme1D<TypeEquation>
  ::GetEnergy(const VectComplex_wp& Zn, const IVect& new_num) const
  {
    VectComplex_wp ProdZn(mesh.GetNbDof()*TypeEquation::nb_unknowns);
    ProdZn.Fill(0);
    MltAddMass(1.0, Zn, Zn, 0, this->var_eq.var_mass, new_num, ProdZn);
    Real_wp energy = abs(DotProdConj(Zn, ProdZn));
    energy = sqrt(energy);
    return energy;
  }
  
  
  //! interpolation on a regular grid
  template<class TypeEquation> template<class T>
  void LdgScheme1D<TypeEquation>
  ::GetInterpolateUn(int m, const Vector<T>& Zn, const IVect& new_num, Vector<T>& ZnInterp) const
  {
    ZnInterp.Reallocate(grid.GetNbPointsGrid());
    int offset = m*mesh.GetNbDof();
    int order = lob.GetOrder();
    for (int i = 0; i < ZnInterp.GetM(); i++)
      {
	T vloc = 0;
	int ne = grid.GetElementNumber(i);
	if (ne >= 0)
	  {
	    Real_wp x = grid.GetLocalCoordinate(i);
	    for (int j = 0; j <= order; j++)
	      vloc += lob.EvaluatePhi(j, x)*Zn(new_num(offset+mesh.GetNumberDof(ne, j)));
	  }
	
	ZnInterp(i) = vloc;
      }
  }
  
  
  //! default constructor
  template<class TypeEquation>
  TimeLdgScheme1D<TypeEquation>::TimeLdgScheme1D() : LdgScheme1D<TypeEquation>()
  {
    type_scheme = TimeSchemeEnum::GAUSS_RUNGE_KUTTA;
    order_time_scheme = 2;
    stopping_criterion = 20.0*epsilon_machine;
    compute_differential = true;
    nb_eval_diff = 0;
  }
  

  //! computation of coefficients for different time schemes
  template<class TypeEquation>
  void TimeLdgScheme1D<TypeEquation>
  ::InitTimeScheme(int type_scheme_, int order_, const Real_wp& dt, const Real_wp& Tf)
  {
    nb_eval_diff = 0;
    type_scheme = type_scheme_;
    order_time_scheme = order_;
    compute_differential = true;
    int Nvol = this->mesh.GetNbDof();
    int nodl = Nvol*TypeEquation::nb_unknowns;
    VectReal_wp Y0(Nvol*TypeEquation::nb_time_unknowns);
    Y0.Fill(0);
    this->var_eq.ComputeInitialCondition(this->mesh.GetXmin(), this->mesh.GetXmax(),
                                         this->Xdof, Y0);
    
    deltat = dt;
        
    NewRowNumber.Reallocate(nodl);
    NewRowNumber.Fill();
    
    switch (type_scheme)
      {
      case TimeSchemeEnum::RUNGE_KUTTA :
	{          
	  RK_scheme.SetOrder(4);
	  RK_scheme.SetInitialCondition(0.0, dt, Y0, *this);
	}
	break;
      case TimeSchemeEnum::GAUSS_RUNGE_KUTTA :
        {
          // computing coefficients of Gauss Runge-Kutta method	  
          int order = order_time_scheme/2-1;	  
          AgaussRK.Reallocate(order+1, order+1);
	  Globatto<Real_wp> gauss;
	  gauss.ConstructQuadrature(order, gauss.QUADRATURE_GAUSS);
	  for (int i = 0; i <= order; i++)
	    for (int j = 0; j <= order; j++)
	      {
		Real_wp vloc = 0;
		for (int k = 0; k <= order; k++)
		  vloc += gauss.Points(j)*gauss.Weights(k)
                    *gauss.EvaluatePhi(i, gauss.Points(j)*gauss.Points(k));
		
		AgaussRK(j, i) = vloc;
	      }
	  
	  BgaussRK = gauss.Weights();
          
	  // computing intermediary derivatives p from initial condition	  
	  Un_time.Reallocate(nodl); 
          Un_time.Fill(0);
	  for (int i = 0; i < Y0.GetM(); i++)
	    Un_time(i) = Y0(i);
                             
          Matrix<Real_wp, General, ArrayRowSparse> matDF;
          matDF.Reallocate(Un_time.GetM(), Un_time.GetM());
          VectReal_wp beta(1);
          beta.Fill(1);
          
          this->GetDifferential(Un_time, Un_time, 0, beta, NewRowNumber, matDF);

          // extracting part relative to intermediary unknowns p
          Matrix<Real_wp, General, ArrayRowSparse> DF;
          int offset = TypeEquation::nb_time_unknowns*Nvol;
          int N = Un_time.GetM() - offset;
          if (N > 0)
            {
              GetSubMatrix(matDF, offset, Un_time.GetM(), DF);
              MatrixLU facto;
              facto.HideMessages();
              GetLU(DF, facto);
              
              // newton algorithm to find p
              Real_wp test = 1e200, test_prec = 1e300; int nb_iter = 0;
              Real_wp norme_rhs = Norm2(Un_time);
              VectReal_wp Pn_sol(N), rhs(N); Pn_sol.Fill(0); rhs.Fill(0);
              VectReal_wp Fn(Un_time.GetM()); Fn.Fill(0);
              while ((nb_iter < 20)&&(test > stopping_criterion*norme_rhs)&&(test < test_prec))
                {
                  for (int i = offset; i < Un_time.GetM(); i++)
                    Un_time(i) = Pn_sol(i-offset);
                  
                  Fn.Fill(0);
                  this->AddSystem(Un_time, Un_time, 0, 1.0, 1.0, NewRowNumber, Fn);
                  
                  for (int i = offset; i < Un_time.GetM(); i++)
                    rhs(i - offset) = Fn(i);
                  
                  test_prec = test;
                  test = Norm2(rhs);
                  SolveLU(facto, rhs);
                  
                  Add(-1.0, rhs, Pn_sol);
                  nb_iter++;
                }
              
            }
          
          GaussRowNumber.Reallocate(nodl*(order+1));
          GaussRowNumber.Fill();
          /*int nb = 0;
          for (int i = 0; i < this->mesh.GetNbElt(); i++)
            for (int j = 0; j <= this->lob.GetOrder(); j++)
              for (int k = 0; k < TypeEquation::nb_unknowns; k++)
                for (int m = 0; m <= order; m++)
                  GaussRowNumber(i*(this->lob.GetOrder()+1)+j+k*Nvol+m*nodl) = nb++;
          */
          // allocating mat diff for further computations
          mat_diff.Clear();
          mat_diff.Reallocate(Un_time.GetM()*(order+1), Un_time.GetM()*(order+1));          
          
	}
	break;
      }
  }
  
  
  template<class TypeEquation>
  void TimeLdgScheme1D<TypeEquation>
  ::EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv,
                               const VectReal_wp& Xn, VectReal_wp& ProdXn,
                               bool invert_mass, bool source)
  {
    // for Runge-Kutta algorithm
    ProdXn.Fill(0);
    this->MltAddStiffness(-1.0, Xn, 0, this->var_eq.var_stiff, NewRowNumber, ProdXn, true);
  }
  

  template<class TypeEquation>
  void TimeLdgScheme1D<TypeEquation>
  ::EvaluateFunction(const Real_wp& tn, const VectReal_wp& Xn, VectReal_wp& ProdXn,
                     bool invert_mass, bool source)
  {
    EvaluateDerivativeFunction(tn, 0, Xn, ProdXn);
  }


  template<class TypeEquation>
  void TimeLdgScheme1D<TypeEquation>::
  CheckDifferential(VectReal_wp& Un_sol, VectReal_wp& Un_half)
  {
    int order = order_time_scheme/2;    
    int nodl = this->mesh.GetNbDof()*TypeEquation::nb_unknowns;
    Matrix<Real_wp> A = AgaussRK; Mlt(deltat, A);
    Matrix<Real_wp> DF_num(Un_sol.GetM(), Un_sol.GetM()); DF_num.Fill(0);
    VectReal_wp ZnTmp = Un_sol;
    VectReal_wp ProdZnp = Un_sol, ProdZnm = Un_sol;
    ProdZnp.Fill(0); ProdZnm.Fill(0); DF_num.Fill(0);
    Real_wp h = 1e-7;
    for (int j = 0; j < Un_sol.GetM(); j++)
      {
        ZnTmp(j) += h;
        for (int m = 0; m < order; m++)
          for (int n = 0; n < nodl; n++)
            {
              Un_half(GaussRowNumber(m*nodl + n)) = Un_time(n);
              for (int k = 0; k < order; k++)
                Un_half(GaussRowNumber(m*nodl + n)) += A(m, k)*ZnTmp(GaussRowNumber(k*nodl + n));
            }
        
        ProdZnp.Fill(0);
        for (int m = 0; m < order; m++)
          this->AddSystem(ZnTmp, Un_half, m*nodl, 1.0, 1.0, GaussRowNumber, ProdZnp);
        
        ZnTmp(j) -= 2.0*h;
        for (int m = 0; m < order; m++)
          for (int n = 0; n < nodl; n++)
            {
              Un_half(GaussRowNumber(m*nodl + n)) = Un_time(n);
              for (int k = 0; k < order; k++)
                Un_half(GaussRowNumber(m*nodl + n)) += A(m, k)*ZnTmp(GaussRowNumber(k*nodl + n));
            }
        
        ProdZnm.Fill(0);
        for (int m = 0; m < order; m++)
          this->AddSystem(ZnTmp, Un_half, m*nodl, 1.0, 1.0, GaussRowNumber, ProdZnm);
        
        for (int i = 0; i < Un_sol.GetM(); i++)
          DF_num(i, j) = (ProdZnp(i) - ProdZnm(i))/(2.0*h);
        
        ZnTmp(j) += h;
      }
    
    DF_num.WriteText("df_num.dat");
    mat_diff.WriteText("df_exact.dat");
    
    for (int j = 0; j < Un_sol.GetM(); j++)
      for (int i = 0; i < Un_sol.GetM(); i++)
        {
          Real_wp val_num = DF_num(i, j);
          Real_wp val_exact = mat_diff(i, j);
          if (abs(val_num-val_exact) > 1e-5)
            {
              DISP(i); DISP(j);
              DISP(val_num); DISP(val_exact);
              // DF(i,j) = val_num;
              abort();
            }
        }                  
  }
   
  
  template<class TypeEquation>
  void TimeLdgScheme1D<TypeEquation>::Advance(const Real_wp& t, int n)
  {
    glob_chrono.Start(VirtualTimer::ALL);
    
    switch (type_scheme)
      {
      case TimeSchemeEnum::RUNGE_KUTTA :
	RK_scheme.Advance(t, n, *this);
	break;
      case TimeSchemeEnum::GAUSS_RUNGE_KUTTA :
	{          
	  int order = order_time_scheme/2;
	  Matrix<Real_wp> A = AgaussRK; Mlt(deltat, A);
	  
          // initial guess = 0
	  VectReal_wp Un_sol(order*Un_time.GetM()); Un_sol.Fill(0);
	  VectReal_wp Un_half = Un_sol, Fn = Un_sol;	  
	  int nodl = this->mesh.GetNbDof()*TypeEquation::nb_unknowns;
          
	  // Newton iteration to solve non-linear system
	  Real_wp test = 1e200;// test_prec = 1e300; 
	  int nb_iter = 0;
          //Real_wp norme_rhs = Norm2(Un_time);
          // we may have already computed f when checking Newton algorithm
          bool rhs_computed = false; Real_wp norm_fn = 1e200; 
          VectReal_wp Un_tmp = Un_sol, FnComputed = Un_sol;
          
          bool test_convergence = false;
	  while (!test_convergence)
            {
	      // computing Un_half = y_n + deltat \sum a_ij z_j
	      for (int i = 0; i < order; i++)
		for (int n = 0; n < nodl; n++)
		  {
		    Un_half(GaussRowNumber(i*nodl + n)) = Un_time(n);
		    for (int j = 0; j < order; j++)
		      Un_half(GaussRowNumber(i*nodl + n))
                        += A(i, j)*Un_sol(GaussRowNumber(j*nodl + n));
		  }
              
              // computation of differential
	      glob_chrono.Start(VirtualTimer::MASS);
              
	      if (compute_differential)
                {
                  mat_diff.Zero();
                  for (int i = 0; i < order; i++)
                    {
                      VectReal_wp beta(order);
                      for (int j = 0; j < order; j++)
                        beta(j) = A(i, j);
                      
                      this->GetDifferential(Un_sol, Un_half, i*nodl, beta,
                                            GaussRowNumber, mat_diff);
                    }                  
                  
                  //CheckDifferential(Un_sol, Un_half);
                }
              
	      glob_chrono.Stop(VirtualTimer::MASS);
              
              // inversion of differential
	      glob_chrono.Start(VirtualTimer::STIFFNESS);
	      
	      if (compute_differential)
                {
                  mat_lu.HideMessages();
                  GetLU(mat_diff, mat_lu, true);
                  compute_differential = false;
                  nb_eval_diff++;
                }
	      
              glob_chrono.Stop(VirtualTimer::STIFFNESS);
              
	      glob_chrono.Start(VirtualTimer::FLUX);
	      
	      //test_prec = test;
              
              if (rhs_computed)
                {
                  Copy(FnComputed, Fn);
                  test = norm_fn;
                }
              else
                {
                  // evaluation of f(Un_half)
                  Fn.Fill(0);
                  for (int i = 0; i < order; i++)
                    this->AddSystem(Un_sol, Un_half, i*nodl, 1.0, 1.0, GaussRowNumber, Fn);
                  
                  test = Norm2(Fn);
                }
              
	      glob_chrono.Stop(VirtualTimer::FLUX);
                            
              glob_chrono.Start(VirtualTimer::STIFFNESS);

              SolveLU(mat_lu, Fn);
	      
              glob_chrono.Stop(VirtualTimer::STIFFNESS);
                            
              // new iterate Xn+1 = Xn - DF^{-1} F
              // or Xn+1 = Xn - gamma*DF^{-1} F in order to have convergent algorithm
              bool test_loop = true;
              Real_wp gamma = 1.0;
              glob_chrono.Start(VirtualTimer::PROD);
              while (test_loop)
                {
                  Copy(Un_sol, Un_tmp);
                  Add(-gamma, Fn, Un_tmp);
                  
                  // computing the residual with this choice of gamma
                  for (int i = 0; i < order; i++)
                    for (int n = 0; n < nodl; n++)
                      {
                        Un_half(GaussRowNumber(i*nodl + n)) = Un_time(n);
                        for (int j = 0; j < order; j++)
                          Un_half(GaussRowNumber(i*nodl + n))
                            += A(i, j)*Un_tmp(GaussRowNumber(j*nodl + n));
                      }
                  
                  FnComputed.Fill(0);
                  for (int i = 0; i < order; i++)
                    this->AddSystem(Un_tmp, Un_half, i*nodl, 1.0, 1.0, GaussRowNumber, FnComputed);
                  
                  norm_fn = Norm2(FnComputed);
                  if (norm_fn < test)
                    {
                      // okay, we have found a gamma for which
                      // next residual is lower
                      test_loop = false;
                      Copy(Un_tmp, Un_sol);
                      rhs_computed = true;
                    }
                  else
                    {
                      if (test < 1e-10)
                        {
                          norm_fn = test;
                          test_convergence = true;
                          test_loop = false;
                        }
                      else
                        {
                          // we decrease gamma
                          gamma *= 0.5;
                                                    
                          // if gamma is too small, we stop the processus
                          if (gamma < 1e-5)
                            {
                              test_loop = false;
                              cout << "Energy " << this->GetEnergy() << endl;
                              cout << "Newton did not converge at all " << endl;
                              abort();
                            }
                        }
                    }
                }
              glob_chrono.Stop(VirtualTimer::PROD);
              
              if (nb_iter >= 50)
                test_convergence = true;
              
              if (norm_fn <= stopping_criterion)
                test_convergence = true;
            
	      nb_iter++;
              
              if (nb_iter%11 == 0)            
                compute_differential = true;
              
	    }
          
          if (nb_iter > 10)            
            compute_differential = true;
          
	  if (norm_fn > 100.0*stopping_criterion)
	    {
	      cout << "Newton iteration stopped at " << test << endl;
	      cout << "Number of iterations " << nb_iter << endl;
	      // abort();
	    }
	  
	  for (int i = 0; i < order; i++)
	    for (int n = 0; n < Un_time.GetM(); n++)
	      Un_time(n) += deltat*BgaussRK(i)*Un_sol(GaussRowNumber(n+i*nodl));
	}
	break;
      }

    glob_chrono.Stop(VirtualTimer::ALL);
  }
  

  template<class TypeEquation>
  Real_wp TimeLdgScheme1D<TypeEquation>::GetEnergy() const
  {
    switch (type_scheme)
      {
      case TimeSchemeEnum::RUNGE_KUTTA :
	return LdgScheme1D<TypeEquation>::GetEnergy(RK_scheme.GetIterate(), NewRowNumber);
	break;
      case TimeSchemeEnum::GAUSS_RUNGE_KUTTA :
	{
	  return LdgScheme1D<TypeEquation>::GetEnergy(Un_time, NewRowNumber);
	}
	break;
      }
    
    return Real_wp(0);
  }
  
  
  template<class TypeEquation>
  const VectReal_wp& TimeLdgScheme1D<TypeEquation>::GetIterate() const
  {
    return GetIterateReal();
  }
  
  
  template<class TypeEquation>
  const VectReal_wp& TimeLdgScheme1D<TypeEquation>::GetIterateReal() const
  {
    switch (type_scheme)
      {
      case TimeSchemeEnum::RUNGE_KUTTA :
	return RK_scheme.GetIterate();
      case TimeSchemeEnum::GAUSS_RUNGE_KUTTA :
	return Un_time;
      }	
    
    return Un_time;
  }
  
  
  template<class TypeEquation>
  void TimeLdgScheme1D<TypeEquation>::GiveIterate(int n, const Real_wp& t, VectReal_wp& Y)
  {
  }
  
  
  template<class TypeEquation>
  void TimeLdgScheme1D<TypeEquation>::GetInterpolateUn(int m, VectReal_wp& ZnInterp) const
  {
    LdgScheme1D<TypeEquation>::GetInterpolateUn(m, this->GetIterate(), NewRowNumber, ZnInterp);
  }

}

#define MONTJOIE_FILE_LOCAL_DISCONTINUOUS_GALERKIN_1D_CXX
#endif
