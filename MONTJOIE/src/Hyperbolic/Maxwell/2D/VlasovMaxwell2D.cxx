#ifndef MONTJOIE_FILE_VLASOV_MAXWELL2D_CXX

#include "VlasovMaxwell2D.hxx"

namespace Montjoie
{
  ////////////////////
  // INITIALIZATION //
  
  //! constructor for 2-D Vlasov-Maxwell problem
  template<class TypeEquation>
  VlasovMaxwell<Dimension2>::VlasovMaxwell(HyperbolicProblem<TypeEquation>& var)
    : VlasovMaxwell_Base<Dimension2>(var)
  {
    uniform_magnetic_field = false;
    value_uniform_magnetic = 0; 
    mirror_right_side = false;
  }
  
  
  //! Reading of an input file line
  void VlasovMaxwell<Dimension2>
  ::SetInputData(const string& description_field, const VectString& parameters)
  {
    VlasovMaxwell_Base<Dimension2>::SetInputData(description_field, parameters);
    
    if (!description_field.compare("UniformMagneticField"))
      {
	this->value_uniform_magnetic = to_num<Real_wp>(parameters(0));
	this->uniform_magnetic_field = true;
      }
    else if (!description_field.compare("MirrorRightSide"))
      {
	if (!parameters(0).compare("YES"))
	  this->mirror_right_side = true;
      }
  }
  
  
  void VlasovMaxwell<Dimension2>
  ComputeProjectors(VectR2& PtsInterpol2D, VectR2& PtsQuad2D, VectReal_wp& Weights2D)
  {
    //const Mesh<Dimension2>& mesh = this->var_harmonic.mesh;
    const MeshNumbering<Dimension2>& mesh_num = this->var_harmonic.mesh_num;
    int r = mesh_num.GetOrder();
    const ElementReference<Dimension2, 1>& Fb = this->var_harmonic.GetReferenceElement(0);
    
    // constructing quadrature points
    QuadrangleQuadrature::
      ConstructQuadrature(2*this->order_integration_current, PtsQuad2D, Weights2D);
    
    Globatto<Real_wp> gauss;
    gauss.ConstructQuadrature(this->order_integration_current);
    VectReal_wp PtsQuad1D = gauss.Points();
    
    // we assume that dofs are defined on Gauss-Lobatto points
    // (edge finite element or discontinuous Galerkin)
    const Matrix<int>& NumNodes2D = Fb.GetNumNodes2D();
    Globatto<Real_wp>& phi_dof = lob_boundary;
    phi_dof.ConstructQuadrature(r, phi_dof.QUADRATURE_LOBATTO);    
    this->proj_nodal_to_quadrature.Init(phi_dof, NumNodes2D, PtsQuad1D, PtsQuad2D);
    
    // basis functions on boundary
    WeightsBoundary = phi_dof.Weights();    
    
    if (this->order_interpolation_current > 0)
      {
	// interpolation points are defined on regular points (or lobatto)
	// with a given number of subdivisions
	SubdivGlobatto phi_interp;
	phi_interp.Init(this->regular_interpolation_current,
                        this->nb_subdiv_interpolation, this->order_interpolation_current);
	int ri = phi_interp.GetOrder();
	PtsInterpol2D.Reallocate((ri+1)*(ri+1));
	Matrix<int> NumNodes2D_int(ri+1, ri+1);
	for (int i = 0; i <= ri; i++)
	  for (int j = 0; j <= ri; j++)
	    {
	      NumNodes2D_int(i, j) = j + i*(ri+1);
	      PtsInterpol2D(NumNodes2D_int(i, j)).Init(phi_interp.Points(i), phi_interp.Points(j));
	    }
	
	this->proj_interpolate_to_quadrature.Init(phi_interp, NumNodes2D_int,
                                                  PtsQuad1D, PtsQuad2D);
      }
  }

  
  //! Initialization of a localization grid to have a fast computation for coupling terms
  /*!
    \param[in] pts interpolation points where the current J need to be evaluated
    \param[in] ElemPoint number of the element to which each interpolation point belongs
    \param[in] LocalPosPoint local coordinates inside the element, for each interpolation point
    \param[in] radius influence radius R for each particle
    \param[out] Cloud_elem list of interpolation points
    (number of element) to scan for each subdivision of the grid
    \param[out] Cloud_pt list of interpolation points
    (local number in the element) to scan for each subdivision of the grid
    \param[out] grid interpolation grid
   */
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_2D<TypeElement,TypeEquation>::
  PreComputeCloudsGrid(VectR2& pts, const Real_wp& radius,
                       Vector<IVect>& Cloud_elem, GridInterpolation<Dimension2>& grid)
  {
    Mesh<Dimension2>& mesh = this->var_harmonic.mesh;
    MeshNumbering<Dimension2>& mesh_num = this->var_harmonic.mesh_num;
    GeneratePeriodicDirections();

    // is there a periodicity condition for x or y axis ?
    bool periodic_x = false, periodic_y = false;
    for (int num = 0; num < mesh.GetNbPeriodicReferences(); num++)
      {
	if (abs(mesh_num.GetTranslationPeriodicBoundary(num)(1)) < 1e-5)
	  periodic_x = true;
	
	if (abs(mesh_num.GetTranslationPeriodicBoundary(num)(0)) < 1e-5)
	  periodic_y = true;
      }
    
    // first step : determination of the grid
    // space step of the grid = R/3
    Real_wp delta = radius/3;
    if (!periodic_x)
      {
	grid.SetXmin(mesh.GetXmin() - 1.5*radius);
	grid.SetXmax(mesh.GetXmax() + 1.5*radius);
      }
    else
      {
	grid.SetXmin(mesh.GetXmin());
	grid.SetXmax(mesh.GetXmax());
      }
    
    if (!periodic_y)
      {
	grid.SetYmin(mesh.GetYmin() - 1.5*radius);
	grid.SetYmax(mesh.GetYmax() + 1.5*radius);
      }
    else
      {
	grid.SetYmin(mesh.GetYmin());
	grid.SetYmax(mesh.GetYmax());
      }
    
    int nbx = max(toInteger(ceil( (grid.GetXmax()-grid.GetXmin())/delta )), 1);
    int nby = max(toInteger(ceil( (grid.GetYmax()-grid.GetYmin())/delta )), 1);
    Real_wp step_x(0), step_y(0);
    step_x = (grid.GetXmax()-grid.GetXmin())/nbx;
    step_y = (grid.GetYmax()-grid.GetYmin())/nby;
    grid.SetNbSubdivisions(nbx, nby);
    grid.SetSubdivisionStep(step_x, step_y);
    
    // second step : all the quadrature points are localized on this grid
    grid.GlobalCoord = pts;
    
    // NumBoxGrid_Point(i) : number of the small box of the grid where the point i is
    // ListPoints_Grid(i) : list of points contained by the small box i
    // NbPoints_Grid(i) : amount of points contained by the small box i
    IVect NumBoxGrid_Point; Vector<IVect> ListPoints_Grid; IVect NbPoints_Grid;
    grid.PreLocalizePoints(NumBoxGrid_Point, ListPoints_Grid, NbPoints_Grid, true);
    
    // third step  for each small box of the regular grid,
    // we look for neighbouring quadrature points
    VectR2 pts_box(8); Real_wp dist; R2 point, pt_glob;
    Cloud_elem.Reallocate(nbx*nby);
    int imin(0), imax(0), jmin(0), jmax(0);
    // DISP(step_x); DISP(step_y); DISP(nbx); DISP(nby);
    // DISP(NbPoints_Grid);
    for (int ix = 0; ix < nbx; ix++)
      for (int iy = 0; iy < nby; iy++)
	{
	  int num_box = iy*nbx + ix;
	  int size_res = 200;
	  Cloud_elem(num_box).Reallocate(size_res);
	  
	  // four vertices of the box
	  pts_box(0)(0) = grid.GetXmin() + ix*step_x; pts_box(0)(1) = grid.GetYmin() + iy*step_y; 
	  pts_box(1)(0) = pts_box(0)(0) + step_x; pts_box(1)(1) = pts_box(0)(1); 
	  pts_box(2)(0) = pts_box(0)(0) + step_x; pts_box(2)(1) = pts_box(0)(1)+step_y; 
	  pts_box(3)(0) = pts_box(0)(0);          pts_box(3)(1) = pts_box(0)(1)+step_y;
	  // four middles of the edges of the box
	  pts_box(4)(0) = pts_box(0)(0) + 0.5*step_x; pts_box(4)(1) = pts_box(0)(1); 
	  pts_box(5)(0) = pts_box(0)(0) + step_x; pts_box(5)(1) = pts_box(0)(1)+0.5*step_y; 
	  pts_box(6)(0) = pts_box(0)(0) + 0.5*step_x; pts_box(6)(1) = pts_box(0)(1)+step_y; 
	  pts_box(7)(0) = pts_box(0)(0);          pts_box(7)(1) = pts_box(0)(1)+0.5*step_y;
	  
	  // DISP(num_box); DISP(pts_box);
	  // loop on surrounding boxes
	  int nb = 0;
	  if (periodic_x)
	    {
	      imin = ix-3;
	      imax = ix+3;
	    }
	  else
	    {
	      imin = max(0,ix-3);
	      imax = min(ix+3,nbx-1);
	    }
	  
	  if (periodic_y)
	    {
	      jmin = iy-3;
	      jmax = iy+3;
	    }
	  else
	    {
	      jmin = max(0,iy-3);
	      jmax = min(iy+3,nby-1);
	    }
	  
	  // DISP(imin); DISP(imax); DISP(jmin); DISP(jmax);
	  for (int i1 = imin; i1 <= imax; i1++)
	    for (int i2 = jmin; i2 <= jmax; i2++)
	      {
		int i1b = i1%nbx;
		int i2b = i2%nby;
		if (i1b < 0)
		  i1b += nbx;
		
		if (i2b < 0)
		  i2b += nby;
		
		int nbox = i2b*nbx + i1b;
		// DISP(NbPoints_Grid(nbox));
		for (int k1 = 0; k1 < NbPoints_Grid(nbox); k1++)
		  {
		    int k = ListPoints_Grid(nbox)(k1);
		    point = grid.GlobalCoord(k);
		    dist = 2.0*radius;
		    for (int m = 0; m < this->periodic_vector.GetM(); m++)
		      {
			Add(point, this->periodic_vector(m), pt_glob);
			for (int p = 0; p < pts_box.GetM(); p++)
			  dist = min(dist, pts_box(p).Distance(pt_glob));
		      }
		    
		    // DISP(point); DISP(dist);
		    if (dist < 1.005*radius)
		      {
			Cloud_elem(num_box)(nb) = k;
			nb++;
			
			if (nb >= size_res)
			  {
			    size_res *= 2;
			    Cloud_elem(num_box).Resize(size_res);
			  }
		      }
		  }
	      }
	  
	  // resizing arrays to their true size
	  Cloud_elem(num_box).Resize(nb);
	  // DISP(num_box); DISP(Cloud_elem(num_box));
	}
    // exit(0);
  }

  
  //! initialization of periodic_vector
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_2D<TypeElement,TypeEquation>::GeneratePeriodicDirections()
  {
    Mesh<Dimension2>& mesh = this->var_harmonic.mesh;
    MeshNumbering<Dimension2>& mesh_num = this->var_harmonic.mesh_num;
    
    // periodicity vectors
    R2 vec, point, vec2;
    this->periodic_vector.Reallocate(1); this->periodic_vector(0).Zero();
    if (mesh.GetNbPeriodicReferences() == 1)
      {
	this->periodic_vector.Reallocate(3);
	this->periodic_vector(0).Zero();
	vec = mesh_num.GetTranslationPeriodicBoundary(0);
	this->periodic_vector(1) = vec;
	Mlt(-1.0, vec); this->periodic_vector(2) = vec;
      }
    else if (mesh.GetNbPeriodicReferences() == 2)
      {
	this->periodic_vector.Reallocate(9);
	this->periodic_vector(0).Zero();
	vec = mesh_num.GetTranslationPeriodicBoundary(0);
	this->periodic_vector(1) = vec;
	Mlt(-1.0, vec); this->periodic_vector(2) = vec;
	vec2 = mesh_num.GetTranslationPeriodicBoundary(1);
	this->periodic_vector(3) = vec2;
	Mlt(-1.0, vec2); this->periodic_vector(4) = vec2;
	
	Add(vec, vec2, this->periodic_vector(5));
	Subtract(vec, vec2, this->periodic_vector(6));
	Mlt(-1.0, vec);
	Add(vec, vec2, this->periodic_vector(7));
	Subtract(vec, vec2, this->periodic_vector(8));
      }
  }
  
  
  //! returns offset in the vector of unknowns to access to the parameters of the first particle
  template<class TypeElement, class TypeEquation>
  int VlasovMaxwell_2D<TypeElement,TypeEquation>::GetOffset_PosParticle() const
  {
    GhostIf<TypeEquationStationary::FormulationDG> dg_form;
    int offset = this->GetNbScalarUnknowns(dg_form) + this->GetNbVectorialUnknowns(dg_form);
    if (this->type_correction_divergence == this->CORRECTION_HYPERBOLIC)
      offset += this->var_laplace.mesh_num.GetNbDof();
    
    return offset;
  }
  
  
  //! computation of E and H on quadrature points
  /*!
    \param[in] En components of E  (dofs)
    \param[in] Hn components of H (dofs)
    \param[out] Equad values of E on quadrature points
    \param[out] Hquad values of H on quadrature points
    \param[in] normal_quadrature integration done on normal quadrature points ?
   */
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_2D<TypeElement,TypeEquation>::
  EvaluateE_H_Quadrature(const VectReal_wp& En, const VectReal_wp& Hn,
			 VectR2& En_quad, VectReal_wp& Hn_quad,
			 bool normal_quadrature, GhostIf<false>& dg_form)
  {
    const Mesh<Dimension2>& mesh = this->var_harmonic.mesh;
    const MeshNumbering<Dimension2>& mesh_num = this->var_harmonic.mesh_num;
    int nb_elt = mesh.GetNbElt();
    int r = mesh_num.GetOrder();
    const TypeElement& Fb = this->var_harmonic.GetFaceBasis(r);
    int N = Fb.GetNbDof();
    int Nquad = this->nb_points_quadrature_current;
    int Nh = Fb.GetNbPointsQuadratureInside();
    if (this->order_integration_current <= 0)
      normal_quadrature = true;
    
    if (normal_quadrature)
      Nquad = Nh;
    
    //EllipticProblem<TypeElement,TypeEquationStationary>& var = this->var_harmonic;
    Real_wp jacob; int Nall = this->PointsQuadrature.GetM();
    if (normal_quadrature)
      Nall = nb_elt*Nquad;
    
    En_quad.Reallocate(Nall); Hn_quad.Reallocate(Nall);
    VectReal_wp Eloc(N), Hloc(Nh);
    VectR2 Enode(Nh), Equad(Nquad);
    R2 vec_u, vec_v; Matrix2_2 mat_dfj; VectReal_wp Hquad(Nquad); 
    
    // loop over elements
    for (int i = 0; i < nb_elt; i++)
      {
	// Eloc -> values of E on local dofs
        abort();
	/*for (int j = 0; j < N; j++)
	  {
	    int num_dof = mesh_num.Element(i).GetNumberDof(j);
	    if (var.IsNegativeSign(i, j))
	      Eloc(j) = -En(num_dof);
	    else
	      Eloc(j) = En(num_dof);
              }*/
	
	// we get values of H on the element i
	for (int j = 0; j < Nh; j++)
	  Hloc(j) = Hn(i*Nh + j);
	
	// computation of E on nodal points
	//Fb.ComputeNodalValuesRef(Eloc, Enode);
        abort();
        
	// projection on quadrature points
	if (!normal_quadrature)
	  {
	    this->proj_nodal_to_quadrature.Project(Enode, Equad);
	    this->proj_nodal_to_quadrature.Project(Hloc, Hquad);
	  }
	else
	  {
	    Equad = Enode;
	    Hquad = Hloc;
	  }
	
	// multiplying by DF_i^{*-1}
	if (normal_quadrature)
	  {
	    if (this->var_harmonic.Glob_DFjm1(i).GetM() == 1)
	      {
		mat_dfj = this->var_harmonic.Glob_DFjm1(i)(0);
		jacob = Det(mat_dfj); jacob = 1.0/jacob;
		for (int j = 0; j < Nquad; j++)
		  {
		    int num = i*Nquad+j;
		    vec_u(0) = Equad(j)(0); vec_u(1) = Equad(j)(1);
		    MltTrans(mat_dfj, vec_u, vec_v);
		    En_quad(num)(0) = jacob*vec_v(0);
		    En_quad(num)(1) = jacob*vec_v(1);
		    Hn_quad(num) = Hquad(j);
		  }
	      }
	    else
	      for (int j = 0; j < Nquad; j++)
		{
		  int num = i*Nquad+j;
		  mat_dfj = this->var_harmonic.Glob_DFjm1(i)(j);
		  jacob = Det(mat_dfj); jacob = 1.0/jacob;
		  vec_u(0) = Equad(j)(0); vec_u(1) = Equad(j)(1);
		  MltTrans(mat_dfj, vec_u, vec_v);
		  En_quad(num)(0) = jacob*vec_v(0);
		  En_quad(num)(1) = jacob*vec_v(1);
		  Hn_quad(num) = Hquad(j);
		}
	  }
	else
	  for (int j = 0; j < Nquad; j++)
	    {
	      int num = i*Nquad+j;
	      mat_dfj = this->Glob_DFjm1(i)(j);
	      jacob = Det(mat_dfj); jacob = 1.0/jacob;
	      vec_u(0) = Equad(j)(0); vec_u(1) = Equad(j)(1);
	      MltTrans(mat_dfj, vec_u, vec_v);
	      En_quad(num)(0) = jacob*vec_v(0);
	      En_quad(num)(1) = jacob*vec_v(1);
	      Hn_quad(num) = Hquad(j);
	    }
      }
  }
  
  
  //! computation of E and H on quadrature points
  /*!
    \param[in] En components of E  (dofs)
    \param[in] Hn components of H (dofs)
    \param[out] Equad values of E on quadrature points
    \param[out] Hquad values of H on quadrature points
    \param[in] normal_quadrature integration done on normal quadrature points ?
   */
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_2D<TypeElement,TypeEquation>::
  EvaluateE_H_Quadrature(const VectReal_wp& En, const VectReal_wp& Hn,
			 VectR2& En_quad, VectReal_wp& Hn_quad,
			 bool normal_quadrature, GhostIf<true>& dg_form)
  {
    const Mesh<Dimension2>& mesh = this->var_harmonic.mesh;
    const MeshNumbering<Dimension2>& mesh_num = this->var_harmonic.mesh_num;
    int nb_elt = mesh.GetNbElt();
    int r = mesh_num.GetOrder();
    const TypeElement& Fb = this->var_harmonic.GetFaceBasis(r);
    int Nquad = this->nb_points_quadrature_current;
    int Nh = Fb.GetNbPointsQuadratureInside();
    int Nvol = mesh_num.GetNbDof();
    if (this->order_integration_current <= 0)
      normal_quadrature = true;
    
    if (normal_quadrature)
      Nquad = Nh;
    
    // EllipticProblem<TypeElement,TypeEquationStationary>& var = this->var_harmonic;
    //Real_wp jacob; 
    int Nall = this->PointsQuadrature.GetM();
    En_quad.Reallocate(Nall); Hn_quad.Reallocate(Nall);
    VectReal_wp Hloc(Nh);
    VectR2 Enode(Nh), Equad(Nquad);
    R2 vec_u, vec_v; Matrix2_2 mat_dfj; VectReal_wp Hquad(Nquad); 
    
    // loop over elements
    for (int i = 0; i < nb_elt; i++)
      {
	// Eloc -> values of E on local dofs
	for (int j = 0; j < Nh; j++)
	  {
	    int num = i*Nh + j;
	    Enode(j)(0) = En(num);
	    Enode(j)(1) = En(num + Nvol);
	    Hloc(j) = Hn(num);
	  }
	
	// projection on quadrature points
	if (!normal_quadrature)
	  {
	    this->proj_nodal_to_quadrature.Project(Enode, Equad);
	    this->proj_nodal_to_quadrature.Project(Hloc, Hquad);
	  }
	else
	  {
	    Equad = Enode;
	    Hquad = Hloc;
	  }
	
	// filling En_quad, Hn_quad
	for (int j = 0; j < Nquad; j++)
	  {
	    int num = i*Nquad + j;
	    En_quad(num) = Equad(j);
	    Hn_quad(num) = Hquad(j);
	  }
      }
  }
  
  
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_2D<TypeElement,TypeEquation>::
  EvaluateE_H_QuadratureBoundary(const VectReal_wp& En, const VectReal_wp& Hn,
				 VectR2& Equad, VectReal_wp& Hquad, GhostIf<true>& dg_form)
  {
    const Mesh<Dimension2>& mesh = this->var_harmonic.mesh;
    //int num_edge, num_elem, nb_pts_quad_edge, nb_dof_edge, num_loc;
    //int r = mesh_num.GetOrder();
    //const TypeElement& Fb = this->var_harmonic.GetFaceBasis(r);
    //Equad.Reallocate((r+1)*mesh.GetNbBoundaryRef());
    //Hquad.Reallocate((r+1)*mesh.GetNbBoundaryRef());
    //int Nh = Fb.GetNbPointsQuadratureInside();
    //int Nvol = mesh_num.GetNbDof();
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      {
        /*
	this->var_harmonic.GetBoundaryValues(mesh, i, num_edge, num_elem,
					     nb_pts_quad_edge, nb_dof_edge, num_loc);
	
	for (int k = 0; k < nb_dof_edge; k++)
	  {
	    int num = Fb.GetQuadNumber(num_loc, k) + num_elem*Nh;
	    int num_quad = i*nb_dof_edge + k;
	    Equad(num_quad)(0) = En(num);
	    Equad(num_quad)(1) = En(num+Nvol);
	    Hquad(num_quad) = Hn(num);
	  }
        */
      }
  }
  
  
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_2D<TypeElement,TypeEquation>::
  EvaluateE_H_QuadratureBoundary(const VectReal_wp& En, const VectReal_wp& Hn,
				 VectR2& Equad, VectReal_wp& Hquad, GhostIf<false>& dg_form)
  {
    const Mesh<Dimension2>& mesh = this->var_harmonic.mesh;
    /*int num_edge, num_elem, nb_pts_quad_edge, nb_dof_edge, num_loc;
    int r = mesh_num.GetOrder();
    const TypeElement& Fb = this->var_harmonic.GetFaceBasis(r);
    VectReal_wp Uloc(Fb.GetNbDof());
    VectR2 Unode(Fb.GetNbPointsNodalElt());
    Equad.Reallocate((r+1)*mesh.GetNbBoundaryRef());
    Hquad.Reallocate((r+1)*mesh.GetNbBoundaryRef());
    int Nh = Fb.GetNbPointsQuadratureInside();
    Real_wp jacob; R2 vec_u; */
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      {
        abort();
	/*this->var_harmonic.GetBoundaryValues(mesh, i, num_edge, num_elem,
					     nb_pts_quad_edge, nb_dof_edge, num_loc);
	
	for (int j = 0; j < Fb.GetNbDof(); j++)
	  {
	    int num_dof = mesh_num.Element(num_elem).GetNumberDof(j);
	    if (this->var_harmonic.IsNegativeSign(num_elem, j))
	      Uloc(j) = -En(num_dof);
	    else
	      Uloc(j) = En(num_dof);
	  }
	
	Fb.ComputeNodalValuesRef(Uloc, Unode);
	
	for (int k = 0; k < nb_dof_edge; k++)
	  {
	    int num = Fb.GetQuadNumber(num_loc, k);
	    int num_quad = i*nb_dof_edge + k;
	    int num2 = 0;
	    if (this->var_harmonic.Glob_DFjm1(num_elem).GetM() > 1)
	      num2 = num;
	    
	    jacob = Det(this->var_harmonic.Glob_DFjm1(num_elem)(num2));
	    MltTrans(this->var_harmonic.Glob_DFjm1(num_elem)(num2), Unode(num), vec_u);
	    Mlt(1.0/jacob, vec_u);
	    Equad(num_quad) = vec_u;
	    Hquad(num_quad) = Hn(num + num_elem*Nh);
	  }
        */
      }
  }
  
  
  //! computation of integral of E and H against interpolation functions basis
  /*!
    \param[in] En components of E  (dofs)
    \param[in] Hn components of H (dofs)
    \param[out] Equad values of E on interpolation points
    \param[out] Hquad values of H on interpolation points
  */
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_2D<TypeElement,TypeEquation>::
  EvaluateE_H_Interpolate(const VectR2& En_quad, const VectReal_wp& Hn_quad,
			  VectR2& En_interpol, VectReal_wp& Hn_interpol)
  {
    Vector<VectReal_wp>& Jacobian = this->Glob_jacobian;
    int nb_pts_interpol = this->nb_points_interpolation_current;
    int Nall = this->PointsInterpolation.GetM();
    En_interpol.Reallocate(Nall); Hn_interpol.Reallocate(Nall);
    int N = this->nb_points_quadrature_current;
    VectR2 Einterpol(nb_pts_interpol), Eloc(N);
    VectReal_wp Hinterpol(nb_pts_interpol), Hloc(N);
    
    for (int i = 0; i < this->var_harmonic.mesh.GetNbElt(); i++)
      {
	for (int j = 0; j < N; j++)
	  {
	    int nquad = i*N + j;
	    Eloc(j) = En_quad(nquad);
	    Mlt(Jacobian(i)(j), Eloc(j));
	    Hloc(j) = Jacobian(i)(j)*Hn_quad(nquad);
	  }
	
	// projection on interpolation points
	this->proj_interpolate_to_quadrature.TransposeProject(Eloc, Einterpol);
	this->proj_interpolate_to_quadrature.TransposeProject(Hloc, Hinterpol);
	
	for (int j = 0; j < nb_pts_interpol; j++)
	  {
	    En_interpol(i*nb_pts_interpol+j) = Einterpol(j);
	    Hn_interpol(i*nb_pts_interpol+j) = Hinterpol(j);
	  }
      }
  }
  
  
  //! Direct integration of current J (evaluated on quadrature points)
  /*!
    \param[in] Jsrc evaluation of J
    \param[in,out] Prod_En we add the integral \f$ -\int J \varphi \f$ to this vector
   */
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_2D<TypeElement,TypeEquation>::
  AddIntegralCurrentJ(const Real_wp& alpha, const VectR2& Jn_quad,
                      VectReal_wp& Prod_En, GhostIf<false>& dg_form)
  {
    const Mesh<Dimension2>& mesh = this->var_harmonic.mesh;
    const MeshNumbering<Dimension2>& mesh_num = this->var_harmonic.mesh_num;
    int nb_elt = mesh.GetNbElt();
    const TypeElement& Fb = this->var_harmonic.GetFaceBasis(mesh_num.GetOrder());
    int N = Fb.GetNbDof();
    int Nquad = this->nb_points_quadrature_current;
    int Nh = Fb.GetNbPointsQuadratureInside();
    
    VectR2 Jquad(Nquad), Jloc(Nh); VectReal_wp Edof(N);
    Matrix2_2 mat_dfj; R2 vec_u, vec_v;
    for (int i = 0; i < nb_elt; i++)
      {
	  for (int j = 0; j < Nquad; j++)
	    {
	      int num = i*Nquad + j;
	      vec_u(0) = Jn_quad(num)(0);
	      vec_u(1) = Jn_quad(num)(1);
	      Mlt(this->Glob_DFjm1(i)(j), vec_u, vec_v);
	      Jquad(j)(0) = vec_v(0)*this->WeightsQuadrature(j);
	      Jquad(j)(1) = vec_v(1)*this->WeightsQuadrature(j);
	    }
	  
	  this->proj_nodal_to_quadrature.TransposeProject(Jquad, Jloc);
	
	  // integration 
	  //Fb.ComputeIntegralRef(Jloc, Edof);
	
	  // summation
          abort();
	  /*for (int j = 0; j < N; j++)
	    {
	      int num_dof = this->var_harmonic.mesh_num.Element(i).GetNumberDof(j);
	      if (this->var_harmonic.IsNegativeSign(i,j))
		Prod_En(num_dof) -= alpha*Edof(j);
	      else
		Prod_En(num_dof) += alpha*Edof(j);
                } */
      }
  }
  
  
  //! Direct integration of current J (evaluated on quadrature points)
  /*!
    \param[in] Jsrc evaluation of J
    \param[in,out] Prod_En we add the integral \f$ -\int J \varphi \f$ to this vector
   */
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_2D<TypeElement,TypeEquation>::
  AddIntegralCurrentJ(const Real_wp& alpha, const VectR2& Jn_quad,
                      VectReal_wp& Prod_En, GhostIf<true>& dg_form)
  {
    const Mesh<Dimension2>& mesh = this->var_harmonic.mesh;
    const MeshNumbering<Dimension2>& mesh_num = this->var_harmonic.mesh_num;
    int nb_elt = mesh.GetNbElt();
    const TypeElement& Fb = this->var_harmonic.GetFaceBasis(mesh_num.GetOrder());
    int Nquad = this->nb_points_quadrature_current;
    int Nh = Fb.GetNbPointsQuadratureInside();
    int Nvol = mesh_num.GetNbDof();
    
    VectR2 Jquad(Nquad), Jloc(Nh);
    for (int i = 0; i < nb_elt; i++)
      {
	  for (int j = 0; j < Nquad; j++)
	    {
	      int num = i*Nquad + j;
	      Jquad(j)= Jn_quad(num);
	      Mlt(alpha*this->Glob_jacobian(i)(j), Jquad(j));
	    }
	  
	  this->proj_nodal_to_quadrature.TransposeProject(Jquad, Jloc);
	
	  // summation
	  for (int j = 0; j < Nh; j++)
	    {
	      int num_dof = this->var_harmonic.mesh_num.Element(i).GetNumberDof(j);
	      Prod_En(num_dof) += Jloc(j)(0);
	      Prod_En(num_dof+Nvol) += Jloc(j)(1);
	    }
      }
  }
  
  
  //! Interpolation of J on quadrature points
  /*!
    \param[in] Jsrc evaluation of J on interpolation points
    \param[in,out] Jquad evaluation of J on quadrature points
   */
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_2D<TypeElement,TypeEquation>::
  EvaluateCurrentQuadrature(const VectR2& Jn_interpol, VectR2& Jn_quad)
  {
    const Mesh<Dimension2>& mesh = this->var_harmonic.mesh;
    int nb_elt = mesh.GetNbElt();
    int Nquad = this->nb_points_quadrature_current;
    int nb_points_interpol = this->nb_points_interpolation_current;
    
    VectR2 Jquad(Nquad), Jloc(nb_points_interpol);
    Jn_quad.Reallocate(nb_elt*Nquad);
    for (int i = 0; i < nb_elt; i++)
      {
	for (int j = 0; j < nb_points_interpol; j++)
	  Jloc(j) = Jn_interpol(i*nb_points_interpol + j);
	
	this->proj_interpolate_to_quadrature.Project(Jloc, Jquad);
	
	for (int j = 0; j < Nquad; j++)
	  Jn_quad(i*Nquad + j) = Jquad(j);
      }
  }
  
  
  //! Computation of points for which P(x) = cte
  /*!
    \param[in] lob definition of basis functions \f$ \varphi_i \f$ 
    \param[in] value polynom \f$ P(x) = \sum value(i) \varphi_i(x) \f$
    \param[in] cte right hand side of equation to solve P(x) = cte
    \param[out] pts found solutions
   */
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_2D<TypeElement,TypeEquation>::
  ComputePoints_ValEqualConstant(const Globatto<Real_wp>& lob, const VectReal_wp& value,
				 const Real_wp& cte, VectReal_wp& pts)
  {
    // newton algorithm on regular points
    Real_wp t0, err, err_prec, f0, df0; int nb_iter; pts.Clear();
    Real_wp threshold = 100*epsilon_machine;
    for (int i = 0; i <= lob.GetOrder(); i++)
      {
	t0 = Real_wp(i)/lob.GetOrder();
	err_prec = 1e30; nb_iter = 0;
	f0 = -cte;
	for (int j = 0; j <= lob.GetOrder(); j++)
	  f0 += lob.EvaluatePhi(j, t0)*value(j);
	
	err = abs(f0);
	err_prec = err+1e5;
	while ((err < (1.0-epsilon_machine)*err_prec)&&(nb_iter <= 20))
	  {
	    // evaluating f'
	    df0 = 0;
	    for (int j = 0; j <= lob.GetOrder(); j++)
	      df0 += lob.EvaluatePhiGrad(j, t0)*value(j);
	    
	    // updating t0
	    t0 -= f0/df0;
	    
	    // we force t0 to belong to [0,1]
	    t0 = min(1.0,t0); t0 = max(0.0,t0);
	    
	    err_prec = err;
	    
	    f0 = -cte;
	    for (int j = 0; j <= lob.GetOrder(); j++)
	      f0 += lob.EvaluatePhi(j, t0)*value(j);
	    
	    err = abs(f0);
	    
	    nb_iter++;
	  }
	
	if (abs(err) < abs(cte)*threshold)
	  pts.PushBack(t0);
      }
    
    // now assembling pts
    if (pts.GetM() >= 2)
      {
	Sort(pts.GetM(), pts);
	// eliminating doublons
	
	Real_wp tprec = pts(0);
	int nb = 1;
	for (int i = 1; i < pts.GetM(); i++)
	  {
	    if (abs(pts(i)-tprec) > threshold)
	      {
		pts(nb) = pts(i);
		tprec = pts(nb);
		nb++;
	      }
	  }
	pts.Resize(nb);
      }
    
  }
  
  //! Computation of \f$ \int rho \varphi \f$ (for Boris correction)
  /*!
    \param[in,out] b_rhs vector to which \f$ \int rho \varphi \f$ is added
    \param[in] evalRho evaluation of rho on quadrature/interpolation points
   */
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_2D<TypeElement,TypeEquation>
  ::EvaluateSourceBoris(VectReal_wp& b_rhs, VectReal_wp& evalRho)
  {
    int nb_pts_interp = this->nb_points_interpolation_current;
    int nb_pts_quad = this->nb_points_quadrature_current;
    int N = this->var_laplace.mesh_num.GetNbLocalDof(0);
    VectReal_wp rho_interp(nb_pts_interp), rho_quad(nb_pts_quad), rho_dof(N);
    b_rhs.Reallocate(this->var_laplace.mesh_num.GetNbDof()); b_rhs.Zero();
    for (int i = 0; i < this->var_harmonic.mesh.GetNbElt(); i++)
      {
	if (this->order_interpolation_current > 0)
	  {
	    for (int j = 0; j < nb_pts_interp; j++)
	      rho_interp(j) = evalRho(i*nb_pts_interp + j);
	    
	    this->proj_interpolate_to_quadrature.Project(rho_interp, rho_quad);
	    
	    for (int j = 0; j < nb_pts_quad; j++)
	      rho_quad(j) *= this->Glob_jacobian(i)(j);
	  }
	else
	  for (int j = 0; j < nb_pts_quad; j++)
	    rho_quad(j) = evalRho(i*nb_pts_quad + j)*this->Glob_jacobian(i)(j);
	
	this->proj_nodal_to_quadrature.TransposeProject(rho_quad, rho_dof);
	
	for (int j = 0; j < N; j++)
	  b_rhs(this->var_laplace.mesh_num.Element(i).GetNumberDof(j)) += rho_dof(j);
      }
    
    this->var_laplace.ImposeNullDirichletCondition(b_rhs);
  }
      
  
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_2D<TypeElement,TypeEquation>::
  EvaluateSourceBoris(VectReal_wp& b_rhs, VectReal_wp& evalRho, VectR2& evalEn,
		      VectR2& evalEnBoundary)
  {
    const Mesh<Dimension2>& mesh = this->var_harmonic.mesh;
    //const MeshNumbering<Dimension2>& mesh_num = this->var_harmonic.mesh_num;
    
    // part \int_\Omega \rho \phi
    EvaluateSourceBoris(b_rhs, evalRho);
    
    // part \int_\Omega E \nabla \phi
    int nb_pts_quad = this->nb_points_quadrature_current;
    int N = this->var_laplace.mesh_num.GetNbLocalDof(0);
    Real_wp coef; VectR2 feval(nb_pts_quad), Edof(N); R2 vec_u, vec_v;
    //int r = mesh_num.GetOrder();
    VectReal_wp contrib(N);
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
	int ref = mesh.Element(i).GetReference();
	Real_wp epsilon = this->var_harmonic.ref_epsilon(ref)(0,0);
	
	for (int j = 0; j < nb_pts_quad; j++)
	  {
	    int num = i*nb_pts_quad + j;
	    vec_v = evalEn(num);
	    coef = this->Glob_jacobian(i)(j)*epsilon;
	    Mlt(coef, vec_v);
	    feval(j) = vec_v;
	  }
	
	this->proj_nodal_to_quadrature.TransposeProject(feval, Edof);
	// this->var_laplace.GetFaceBasis(r).ComputeIntegralGradientRef(Edof, contrib);
	abort();
        
	for (int j = 0; j < N; j++)
	  b_rhs(this->var_laplace.mesh_num.Element(i).GetNumberDof(j)) += contrib(j);
      }
    
    // part -\int_\Gamma E.n \varphi
    //int num_edge, num_elem, nb_pts_quad_edge, nb_dof_edge, num_loc;
    for (int num_edge_ref = 0; num_edge_ref < mesh.GetNbBoundaryRef(); num_edge_ref++)
      {
        /*
	this->var_laplace.GetBoundaryValues(mesh, num_edge_ref, num_edge,
					    num_elem, nb_pts_quad_edge, nb_dof_edge, num_loc);
	
	int ref = mesh.Element(num_elem).GetReference();
	Real_wp epsilon = this->var_harmonic.ref_epsilon(ref)(0,0);
	R2 normale;
	
	if (mesh.Boundary(num_edge).GetNbElements() < 2)
	  {
	    for (int k = 0; k < nb_dof_edge; k++)
	      {
		int num = num_edge_ref*nb_dof_edge + k;
		vec_u = evalEnBoundary(num);
		normale = this->Glob_normale(num_edge_ref)(k);
		
		coef = -DotProd(normale, vec_u);
		coef *= epsilon*this->WeightsBoundary(k);
		
		int num_dof = this->var_laplace.mesh_num.Element(num_elem)
                .GetNumberDof(this->var_laplace.GetLocalNumber(num_elem, num_loc, k));
		
		b_rhs(num_dof) += coef;
	      }
	  }
        */
      }
    
    this->var_laplace.ImposeNullDirichletCondition(b_rhs);
  }
  
  
  //! adding integral \f$ \alpha * \int_\Gamma \phi \varphi \cdot n \f$
  /*!
    \param[in] alpha coefficient
    \param[in] phi values of phi (phi lies on nodal space H^1)
    \param[in,out] Prod_En vector to which integral is added
   */
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_2D<TypeElement,TypeEquation>::
  AddSurfacicIntegral_PhiNodal(const Real_wp& alpha, const VectReal_wp& phi,
                               VectReal_wp& Prod_En, GhostIf<false>& dg_form)
  {
    //int nb_pts_quad_edge, nb_dof_edge, num_loc, num_edge, num_elem, num_dof;
    const Mesh<Dimension2>& mesh = this->var_laplace.mesh;
    //nb_dof_edge = mesh_num.GetOrder()+1;
    //int nb_dof_loc = this->var_laplace.mesh_num.GetNbLocalDof(0);
    //VectR2 feval(nb_dof_edge); VectReal_wp res(nb_dof_loc);
    //const TypeElement& Fb = this->var_harmonic.GetFaceBasis(mesh_num.GetOrder());
    //R2 vec_u, vec_v; Real_wp val;
    for (int num_edge_ref = 0; num_edge_ref < mesh.GetNbBoundaryRef(); num_edge_ref++)
      {
        abort();
	/* this->var_laplace.GetBoundaryValues(mesh, num_edge_ref, num_edge,
					    num_elem, nb_pts_quad_edge, nb_dof_edge, num_loc);
	
	if (mesh.Boundary(num_edge).GetNbElements() < 2)
	  {
	    feval.Fill(R2(0,0));
	    
	    for (int k = 0; k < nb_dof_edge; k++)
	      {
		int num_dof_loc = this->var_laplace.GetLocalNumber(num_elem, num_loc, k);
		num_dof = mesh_num.Element(num_elem).GetNumberDof(num_dof_loc);
		val = phi(num_dof);
		
		int num = num_edge*nb_dof_edge + k;
		vec_u = this->Glob_normale(num_edge)(k);
		Mlt(alpha*WeightsBoundary(k)*val, vec_u);
		
		// multiplying by DFjm1
		Mlt(this->Glob_DFjm1_Boundary(num_edge)(k), vec_u, feval(k));
	      }
	    
	    Fb.ComputeIntegralAllSurfaceRef(feval, res, num_loc);
	    for (int j = 0; j < nb_dof_loc; j++)
	      {
		num_dof = this->var_harmonic.mesh_num.Element(num_elem).GetNumberDof(j);
		if (this->var_harmonic.IsNegativeSign(num_elem,j))
		  Prod_En(num_dof) -= res(j);
		else
		  Prod_En(num_dof) += res(j);
	      }
	  }
        */
      }
  }
  
  
  //! adding integral \f$ \alpha * \int_\Gamma \phi \varphi \cdot n \f$
  /*!
    \param[in] alpha coefficient
    \param[in] phi values of phi (phi lies on nodal space H^1)
    \param[in,out] Prod_En vector to which integral is added
   */
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_2D<TypeElement,TypeEquation>::
  AddSurfacicIntegral_PhiNodal(const Real_wp& alpha, const VectReal_wp& phi,
                               VectReal_wp& Prod_En, GhostIf<true>& dg_form)
  {
    //int nb_pts_quad_edge, nb_dof_edge, num_loc, num_edge, num_elem, num_dof;
    const Mesh<Dimension2>& mesh = this->var_laplace.mesh;
    //nb_dof_edge = mesh_num.GetOrder()+1;
    //int nb_dof_loc = this->var_laplace.mesh_num.GetNbLocalDof(0);
    //VectReal_wp fx(nb_dof_edge), fy(nb_dof_edge);
    //VectReal_wp res_x(nb_dof_loc), res_y(nb_dof_loc);
    //const TypeElement& Fb = this->var_harmonic.GetFaceBasis(mesh_num.GetOrder());
    //R2 vec_u;
    //int Nvol = mesh_num.GetNbDof();
    for (int num_edge_ref = 0; num_edge_ref < mesh.GetNbBoundaryRef(); num_edge_ref++)
      {
        abort();
	/*this->var_laplace.GetBoundaryValues(mesh, num_edge_ref, num_edge,
					    num_elem, nb_pts_quad_edge, nb_dof_edge, num_loc);
	
	if (mesh.Boundary(num_edge).GetNbElements() < 2)
	  {
	    fx.Fill(0); fy.Fill(0);
	    
	    for (int k = 0; k < nb_dof_edge; k++)
	      {
		int num_dof_loc = this->var_laplace.GetLocalNumber(num_elem, num_loc, k);
		num_dof = mesh_num.Element(num_elem).GetNumberDof(num_dof_loc);
		Real_wp val = phi(num_dof);
		
		int num = num_edge*nb_dof_edge + k;
		vec_u = this->Glob_normale(num_edge)(k);
		Mlt(alpha*WeightsBoundary(k)*val, vec_u);
		
		fx(k) = vec_u(0);
		fy(k) = vec_u(1);
	      }
	    
	    Fb.ComputeIntegralSurface_Ref(fx, res_x, num_loc);
	    Fb.ComputeIntegralSurface_Ref(fy, res_y, num_loc);
	    
	    for (int j = 0; j < nb_dof_loc; j++)
	      {
		num_dof = num_elem*nb_dof_loc + j;
		Prod_En(num_dof) += res_x(j); 
		Prod_En(Nvol+num_dof) += res_y(j); 
	      }
	  }
        */
      }
  }
  
  
  //! adding \f$ \int \nabla \phi \cdot \varphi \f$
  /*!
    \param[in] alpha coefficient
    \param[in] phi vector phi which lies on the nodal space H^1
    \param[in,out] Prod_En vector to which integral is added
   */
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_2D<TypeElement,TypeEquation>::
  AddVolumetricIntegral_PhiNodal(const Real_wp& alpha, const VectReal_wp& phi,
                                 VectReal_wp& Prod_En, GhostIf<false>& dg_form)
  {
    // number of elements and quadrature points on each element
    const Mesh<Dimension2>& mesh = this->var_harmonic.mesh;
    const MeshNumbering<Dimension2>& mesh_num = this->var_harmonic.mesh_num;
    const TypeElement& Fb = this->var_harmonic.GetFaceBasis(mesh_num.GetOrder());
    const QuadrangleLobatto& Fb_laplace = this->var_laplace.GetFaceBasis(mesh_num.GetOrder());
    int nb_elt = mesh.GetNbElt();
    
    // variables for interpolation
    int r = mesh_num.GetOrder(), N = r+1;
    int nb_dof = N*N;
    VectReal_wp Uloc(nb_dof);
    Vector<R2> Enode(nb_dof), grad_phi(nb_dof);
    VectReal_wp contrib(Fb.GetNbDof());
    R2 vec_u, vec_v, vec_w; Matrix2_2 dfjm1;
    
    for (int i = 0; i < nb_elt; i++)
      {
	// we get values of phi on the element i
	for (int j = 0; j < nb_dof; j++)
	  Uloc(j) = phi(this->var_laplace.mesh_num.Element(i).GetNumberDof(j));

	// derivation	
	Fb_laplace.ComputeNodalGradientRef(Uloc, grad_phi);
	
	// multiplying by DF_i^{-1} DF_i^{*-1}
	for (int j = 0; j < nb_dof; j++)
	  {
	    dfjm1 = this->var_harmonic.Glob_DFjm1(i)(j);
	    Real_wp jacob = Det(dfjm1);
	    vec_u = grad_phi(j);
	    MltTrans(dfjm1, vec_u, vec_w);
	    Mlt(dfjm1, vec_w, vec_v);
	    jacob = alpha*this->var_harmonic.Glob_jacobian(i)(j)/(jacob*jacob);
	    Enode(j)(0) = jacob*vec_v(0); Enode(j)(1) = jacob*vec_v(1);
	  }
	
	//Fb.ComputeIntegralRef(Enode, contrib);
	
        abort();
	/* for (int j = 0; j < Fb.GetNbDof(); j++)
	  {
	    int num_dof = this->var_harmonic.mesh_num.Element(i).GetNumberDof(j);
	    if (this->var_harmonic.IsNegativeSign(i,j))
	      Prod_En(num_dof) -= contrib(j);
	    else
	      Prod_En(num_dof) += contrib(j);
              }*/
	
      }
    
  }
  
  
  //! adding \f$ \int \nabla \phi \cdot \varphi \f$
  /*!
    \param[in] alpha coefficient
    \param[in] phi vector phi which lies on the nodal space H^1
    \param[in,out] Prod_En vector to which integral is added
   */
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_2D<TypeElement,TypeEquation>::
  AddVolumetricIntegral_PhiNodal(const Real_wp& alpha, const VectReal_wp& phi,
                                 VectReal_wp& Prod_En, GhostIf<true>& dg_form)
  {
    // number of elements and quadrature points on each element
    const Mesh<Dimension2>& mesh = this->var_harmonic.mesh;
    const MeshNumbering<Dimension2>& mesh_num = this->var_harmonic.mesh_num;
    const TypeElement& Fb = this->var_harmonic.GetFaceBasis(mesh_num.GetOrder());
    const QuadrangleLobatto& Fb_laplace = this->var_laplace.GetFaceBasis(mesh_num.GetOrder());
    int nb_elt = mesh.GetNbElt();
    
    // variables for interpolation
    int r = mesh_num.GetOrder(), N = r+1;
    int nb_dof = N*N;
    VectReal_wp Uloc(nb_dof);
    Vector<R2> Enode(nb_dof), grad_phi(nb_dof);
    R2 vec_u, vec_v; Matrix2_2 dfjm1;
    int Nvol = mesh_num.GetNbDof();
    
    for (int i = 0; i < nb_elt; i++)
      {
	// we get values of phi on the element i
        abort();
	//for (int j = 0; j < nb_dof; j++)
        //Uloc(j) = phi(this->var_laplace.mesh_num.GetNumberDof(i,j));

	// derivation	
	Fb_laplace.ComputeNodalGradientRef(Uloc, grad_phi);
	
	// multiplying by DF_i^{*-1}
	for (int j = 0; j < nb_dof; j++)
	  {
	    dfjm1 = this->var_harmonic.Glob_DFjm1(i)(j);
	    Real_wp jacob = Det(dfjm1);
	    vec_u = grad_phi(j);
	    MltTrans(dfjm1, vec_u, vec_v);
	    jacob = alpha*Fb.WeightsND(j);
	    Enode(j)(0) = jacob*vec_v(0); Enode(j)(1) = jacob*vec_v(1);
	  }
	
	for (int j = 0; j < Fb.GetNbDof(); j++)
	  {
	    int num_dof = this->var_harmonic.mesh_num.Element(i).GetNumberDof(j);
	    Prod_En(num_dof) += Enode(j)(0);
	    Prod_En(Nvol+num_dof) += Enode(j)(1);
	  }
	
      }
    
  }
  
  
  //! adding \f$ \nabla \phi \f$
  /*!
    \param[in] alpha coefficient
    \param[in] phi components of phi on nodal points (H^1)
    \param[in,out] Prod_En vector to which \f$ \alpha \nabla \phi \f$ is added
   */
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_2D<TypeElement,TypeEquation>::
  AddVolumetricProjection_PhiNodal(const Real_wp& alpha, const VectReal_wp& phi,
				VectReal_wp& Prod_En, GhostIf<false>& dg_form)
  {
    // number of elements and quadrature points on each element
    Mesh<Dimension2>& mesh = this->var_harmonic.mesh;
    MeshNumbering<Dimension2>& mesh_num = this->var_harmonic.mesh_num;
    const TypeElement& Fb = this->var_harmonic.GetFaceBasis(mesh_num.GetOrder());
    const QuadrangleLobatto& Fb_laplace = this->var_laplace.GetFaceBasis(mesh_num.GetOrder());
    int nb_elt = mesh.GetNbElt();
    
    VectBool DofUsed(this->var_harmonic.mesh_num.GetNbDof()); DofUsed.Fill(false);
    VectReal_wp contrib(Fb.GetNbDof()); contrib.Fill(0);
    VectR2 grad_phi(Fb_laplace.GetNbDof()), Edof(Fb_laplace.GetNbDof());
    VectReal_wp Uloc(Fb_laplace.GetNbDof());
    
    for (int i = 0; i < nb_elt; i++)
      {
	// we get values of phi on the element i
	for (int j = 0; j < Fb_laplace.GetNbDof(); j++)
	  Uloc(j) = phi(this->var_laplace.mesh_num.Element(i).GetNumberDof(j));
	
	// computation of gradient
	//Fb_laplace.ComputeNodalGradientRef(Uloc, grad_phi);
	
	// projection on dofs
	//Fb.ProjectQuadratureToDofRef(Fb, grad_phi, Edof);
	//Fb.ComputeProjectionDofRef(Edof, contrib);
	
        abort();
        /*
	for (int j = 0; j < Fb.GetNbDof(); j++)
	  if (!DofUsed(mesh_num.Element(i).GetNumberDof(j)))
	    {
	      DofUsed(mesh_num.Element(i).GetNumberDof(j)) = true;
	      int num_dof = mesh_num.Element(i).GetNumberDof(j);
	      if (this->var_harmonic.IsNegativeSign(i,j))
		Prod_En(num_dof) -= alpha*contrib(j);
	      else
		Prod_En(num_dof) += alpha*contrib(j);
	    }
        */
      }
  }

  
  //! adding \f$ \nabla \phi \f$
  /*!
    \param[in] alpha coefficient
    \param[in] phi components of phi on nodal points (H^1)
    \param[in,out] Prod_En vector to which \f$ \alpha \nabla \phi \f$ is added
  */
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_2D<TypeElement,TypeEquation>::
  AddVolumetricProjection_PhiNodal(const Real_wp& alpha, const VectReal_wp& phi,
				VectReal_wp& Prod_En, GhostIf<true>& dg_form)
  {
    // number of elements and quadrature points on each element
    Mesh<Dimension2>& mesh = this->var_harmonic.mesh;
    MeshNumbering<Dimension2>& mesh_num = this->var_harmonic.mesh_num;
    const TypeElement& Fb = this->var_harmonic.GetFaceBasis(mesh_num.GetOrder());
    const QuadrangleLobatto& Fb_laplace = this->var_laplace.GetFaceBasis(mesh_num.GetOrder());
    int nb_elt = mesh.GetNbElt();
    
    VectReal_wp Uloc(Fb_laplace.GetNbDof());
    VectR2 grad_phi(Fb_laplace.GetNbDof());
    R2 vec_u; int Nvol = mesh_num.GetNbDof();
    for (int i = 0; i < nb_elt; i++)
      {
	// we get values of phi on the element i
	for (int j = 0; j < Fb_laplace.GetNbDof(); j++)
	  Uloc(j) = phi(this->var_laplace.mesh_num.Element(i).GetNumberDof(j));
	
	// computation of gradient
	Fb_laplace.ComputeNodalGradientRef(Uloc, grad_phi);
		
	for (int j = 0; j < Fb.GetNbDof(); j++)
	  {
	    //Real_wp jacob = Det(this->var_harmonic.Glob_DFjm1(i)(j));
	    MltTrans(this->var_harmonic.Glob_DFjm1(i)(j), grad_phi(j), vec_u);
	    int num_dof = mesh_num.Element(i).GetNumberDof(j);
	    Prod_En(num_dof) += alpha*vec_u(0);
	    Prod_En(num_dof+Nvol) += alpha*vec_u(1);
	  }
      }
  }
  
  
  //! projection of rho on nodal dofs
  /*!
    \param[out] phi result of projection
    \param[in] evalRho evaluation of rho on quadrature/interpolate points
  */
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_2D<TypeElement,TypeEquation>
  ::ComputeValue_PhiNodal(VectReal_wp& phi, const VectReal_wp& evalRho)
  {
    int nb_elt = this->var_harmonic.mesh.GetNbElt();
    int N = this->var_laplace.mesh_num.GetNbLocalDof(0);
    int Nquad = this->nb_points_quadrature_current;
    int Ninterp = this->nb_points_interpolation_current;
    VectReal_wp Udof(N), Uquad(Nquad), Uinterp(Ninterp);
    
    phi.Reallocate(this->var_laplace.GetNbDof()); phi.Fill(0);
    for (int i = 0; i < nb_elt; i++)
      {
	if (this->order_interpolation_current > 0)
	  {
	    for (int j = 0; j < Ninterp; j++)
	      Uinterp(j) = evalRho(i*Ninterp + j);
	    
	    this->proj_interpolate_to_quadrature.Project(Uinterp, Uquad);
	  }
	else
	  for (int j = 0; j < Nquad; j++)
	    Uquad(j) = evalRho(i*Nquad + j);
	
	for (int j = 0; j < Nquad; j++)
	  Uquad(j) *= this->Glob_jacobian(i)(j);
	
	this->proj_nodal_to_quadrature.TransposeProject(Uquad, Udof);
	
	for (int j = 0; j < N; j++)
	  {
	    int num_dof = this->var_laplace.mesh_num.Element(i).GetNumberDof(j);
	    phi(num_dof) += Udof(j);
	  }
      }
    
    for (int i = 0; i < this->var_laplace.mesh_num.GetNbDof(); i++)
      phi(i) *= this->inv_mass_matrix(i);
    
  }
  
  
  //! setting the magnetic field to an uniform value
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_2D<TypeElement,TypeEquation>::
  SetUniform_MagneticField(const Real_wp& valH, VectReal_wp& Y, GhostIf<false>& dg_form)
  {
    int N = this->var_harmonic.GetNbPointsQuadratureInside(0);
    // uniform magnetic field
    for (int i = 0; i < this->var_harmonic.mesh.GetNbElt(); i++)
      for (int j = 0; j < N; j++)
	Y(this->var_harmonic.mesh_num.GetNbDof() + i*N + j) = valH;
  }
  
  
  //! setting the magnetic field to an uniform value
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_2D<TypeElement,TypeEquation>::
  SetUniform_MagneticField(const Real_wp& valH, VectReal_wp& Y, GhostIf<true>& dg_form)
  {
    int N = this->var_harmonic.GetNbPointsQuadratureInside(0);
    // uniform magnetic field
    for (int i = 0; i < this->var_harmonic.mesh.GetNbElt(); i++)
      for (int j = 0; j < N; j++)
	Y(2*this->var_harmonic.mesh_num.GetNbDof() + i*N + j) = valH;
  }
  
  
  //! not used
  template<class TypeElement,class TypeEquation>
  Real_wp VlasovMaxwell_2D<TypeElement,TypeEquation>
  ::GetNormeSolution(const VectReal_wp& Uh) const
  {
    return 0;
  }
  
  
  // INITIALIZATION //
  ////////////////////
  
  ////////////////////////
  // RUNNING TIME STEPS //
  
  //! time iterations for leap frog scheme
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_2D<TypeElement,TypeEquation>::RunTimeIterations()
  {
    HyperbolicProblem<TypeElement, TypeEquation>& var_leaf = this->GetLeafClass();
    CurrentSource_Vlasov<TypeElement,TypeEquation> Jsrc(var_leaf);
    GhostIf<TypeEquationStationary::FormulationDG> dg_form;
    IVect NumBoxGrid_Point;
    int size_E = this->GetNbScalarUnknowns(dg_form);
    int size_H = this->GetNbVectorialUnknowns(dg_form);
    int offset = this->GetOffset_PosParticle();
    
    int Npart, num_elem, nb_iter = 0;
    R2 pt_loc, pk, vk, pa, va, vb;
    Real_wp poids_acc, poids_J, t = this->initial_time;
    VectR2 evalEn, evalEnBoundary; VectReal_wp evalHn, evalHnBoundary;
    Real_wp val_Ex, val_Ey, val_H, dt = this->deltat;
    int nb; Real_wp q2m_dt, gamma, gamma2;
    Real_wp coef1, coef2, delta, coef_diag, coef_extra, one(1), zero(0);
    
    GridInterpolation<Dimension2>& grid = this->loc_particle;
    VectR2 Jn_quad, En_quad; VectReal_wp Hn_quad;
    if (this->GetTimeSchemeType() == TimeSchemeEnum::LEAP_FROG_SYSTEM)
      {
	// allocation of iteration vectors
	VectReal_wp X, Y, Prod_En, Prod_Hn, En, Hn; VectReal_wp Phin, Prod_Phi;
	VectReal_wp Enodal, Hnodal, Jnodal;
	
	X.Reallocate(this->nodl); X.Fill(0);
	this->SetInitialVector(this->initial_time, X);
	Y.Reallocate(X.GetM()); Y.Fill(0);
	
	while (t <= this->final_time+this->epsilon_time)
	  {
	    this->GiveIterate(nb_iter, t, X);
	    // if new particles, increasing the size of arrays
	    if (X.GetM() > Y.GetM())
	      Y.Reallocate(X.GetM());
	    
	    En.SetData(size_E, &X(0)); Hn.SetData(size_H, &X(size_E));
	    Prod_En.SetData(size_E, &Y(0)); Prod_Hn.SetData(size_H, &Y(size_E));
	    
	    // first step : advancing electric field E
	    // (E^{n+1} - E^n)/ dt = R_h B^{n+1/2} - J^{n+1/2}
	    
	    // stiffness matrix for E
	    var_leaf.ApplyOperatorRhScalar(one, t, Hn, Real_wp(0), Prod_En);
	    // source term
	    var_leaf.AddScalarSourceAtTime(one, t, 0, Prod_En);

	    // adding the term - J^{n-1/2}
	    this->LocalizeParticles(X, NumBoxGrid_Point);
	    
	    Npart = 0; Jsrc.evalJ.Fill(R2()); Jsrc.evalRho.Fill(0);
	    nb = offset;
	    for (int i = 0; i < this->nb_particles.GetM(); i++)
	      {
		poids_acc = this->ratio_qm_species(i);
		for (int j = 0; j < this->nb_particles(i); j++)
		  {
		    num_elem = NumBoxGrid_Point(Npart);
		    if (num_elem >= 0)
		      {
			// particle inside the domain
			
			// we get the velocity of the particle
			// quantity of movment
			pk(0) = X(nb+4*j+2); pk(1) = X(nb+4*j+3);
			// deducing the velocity
			this->GetVelocity_FromMomentum(pk, vk); 
			
			// updating value of J
			poids_J = this->weight_particle(i)(j)*this->charge_species(i);
			Jsrc.UpdateCurrent(poids_J, vk, grid.GlobalCoord(Npart),
					   num_elem, pt_loc);
		      }
		    Npart++;
		  }
		nb += 4*this->nb_max_particles(i);
	      }
	    
	    // term of hyperbolic correction S_h \phi^{n+1/2}
	    if (this->type_correction_divergence == this->CORRECTION_HYPERBOLIC)
	      {
		Phin.SetData(this->var_laplace.mesh_num.GetNbDof(), &X(size_E+size_H));
		
		var_leaf.AddVolumetricIntegral_PhiNodal(-one, Phin, Prod_En, dg_form);
		var_leaf.AddSurfacicIntegral_PhiNodal(one, Phin, Prod_En, dg_form);
		
		Phin.Nullify();
	      }
	    
	    if (this->nb_points_interpolation_current > 0)
	      {
		var_leaf.EvaluateCurrentQuadrature(Jsrc.evalJ, Jn_quad);
		var_leaf.AddIntegralCurrentJ(-one, Jn_quad, Prod_En, dg_form);
	      }
	    else
	      var_leaf.AddIntegralCurrentJ(-one, Jsrc.evalJ, Prod_En, dg_form);
	    
	    // inverting by the mass matrix
	    var_leaf.ApplyOperatorDhMinusdtSh(one, t+0.5*dt, En, dt, Prod_En);
	    Copy(Prod_En, En); var_leaf.SolveOperatorDhPlusdtSh(En);
	    var_leaf.SetDirichletCondition(t+dt, 0, En);
	    
	    // advancing scheme in
	    // (H^{n+3/2} - H^{n+1/2})/dt + ShV (H^{n+3/2} + H^{n+1/2})/2 = RhV E^{n+1} 
	    var_leaf.ApplyOperatorRhVectorial(dt, t+dt, En, zero, Prod_Hn);
	    var_leaf.AddVectorialSourceAtTime(dt, t+dt, 0, Prod_Hn);
	    
	    // Prod_Hn = H^{n+3/2}
	    var_leaf.ApplyOperatorBhMinusdtSh(one, t+dt, Hn, one, Prod_Hn);
	    var_leaf.SolveOperatorBhPlusdtSh(Prod_Hn);
	    // Hn = H^{n+1}
	    Add(one, Prod_Hn, Hn); Mlt(0.5, Hn);
	    	    
	    // evaluating E^{n+1} and H^{n+1}
	    if (this->nb_points_interpolation_current > 0)
	      {
		var_leaf.EvaluateE_H_Quadrature(En, Hn, En_quad, Hn_quad, false, dg_form);
		var_leaf.EvaluateE_H_Interpolate(En_quad, Hn_quad, evalEn, evalHn);
	      }
	    else
	      var_leaf.EvaluateE_H_Quadrature(En, Hn, evalEn, evalHn, false, dg_form);
	    
	    Copy(Prod_Hn, Hn);
	    
	    // now advancing the particles
	    Npart = 0; nb = offset;
	    Jsrc.evalJ.Fill(0); Jsrc.evalRho.Fill(0);
	    for (int i = 0; i < this->nb_particles.GetM(); i++)
	      {
		poids_acc = this->ratio_qm_species(i);
		q2m_dt = 0.5*poids_acc*this->deltat;
		for (int j = 0; j < this->nb_particles(i); j++)
		  {
		    // int num_elem = grid.ElementInterp(Npart);
		    // pt_loc = grid.CoorInterp(Npart);
		    num_elem = NumBoxGrid_Point(Npart);
		    if (num_elem >= 0)
		      {
			// particle inside the domain
			
			// we get the velocity of the particle
			// quantity of movment
			pk(0) = X(nb+4*j+2); pk(1) = X(nb+4*j+3);
			// deducing the velocity
			this->GetVelocity_FromMomentum(pk, vk); 
			
			// updating value of J
			poids_J = this->weight_particle(i)(j)*this->charge_species(i);
			Jsrc.UpdateCurrent(poids_J, vk, grid.GlobalCoord(Npart),
					   num_elem, pt_loc);
			
			val_Ex = 0; val_Ey = 0; val_H = 0;
			for (int num_point = 0; num_point <
                               this->CloudElem_QuadraturePoint(num_elem).GetM(); num_point++)
			  {
			    Real_wp coefS = Jsrc.CoefficientShape(num_point);
			    if (coefS != Real_wp(0))
			      {
				int num = this->CloudElem_QuadraturePoint(num_elem)(num_point);
				val_Ex += coefS*evalEn(num)(0);
				val_Ey += coefS*evalEn(num)(1);
				val_H += coefS*evalHn(num);
			      }
			    
			  }
			
			val_H *= this->var_harmonic.mu0;
			
			// advancing momentum p
			// (p^{n+3/2} - p^{n+1/2}) = q/m (E^{n+1} + v \times mu_0 H^{n+1})
			
			// p_a = p^{n+1/2} + q/2m dt E^{n+1} (advancing of E of dt/2)
			pa(0) = pk(0) + q2m_dt*val_Ex;
			pa(1) = pk(1) + q2m_dt*val_Ey;
			
			// now solving dp/dt = q/m v \times B (rotation)
			// (p_b - p_a)/dt = q/m (v_a + v_b)/2 \times B
			this->GetVelocity_FromMomentum(pa, va, gamma, gamma2);
			coef1 = q2m_dt*val_H; coef2 = coef1*coef1;
			delta = 1.0/(gamma2 + coef2);
			coef_diag = delta*(gamma2 - coef2); coef_extra = delta*2.0*gamma*coef1;
			vb(0) = coef_diag*va(0) + coef_extra*va(1);
			vb(1) = coef_diag*va(1) - coef_extra*va(0);
						
			// p^{n+3/2} = p_b + q/2m dt E^{n+1} (advancing of E of dt/2)
			pk(0) = gamma*vb(0) + q2m_dt*val_Ex;
			pk(1) = gamma*vb(1) + q2m_dt*val_Ey; // DISP(pk);
			X(nb+4*j+2) = pk(0); X(nb+4*j+3) = pk(1);
			
			// (x^{n+2}-x^{n+1})/dt = v^{n+3/2}
			this->GetVelocity_FromMomentum(pk, vk); //DISP(vk);
			X(nb+4*j) += this->deltat*vk(0); 
			X(nb+4*j+1) += this->deltat*vk(1);
			
		      }
		    Npart++;
		  }
		nb += 4*this->nb_max_particles(i);
	      }
	    
	    // term of hyperbolic correction (\phi^{n+3/2} - \phi^{n+1/2})/dt
            //   = S_h^* E^{n+1} + \rho^{n+1}
	    if (this->type_correction_divergence == this->CORRECTION_HYPERBOLIC)
	      {
		Phin.SetData(this->var_laplace.mesh_num.GetNbDof(), &X(size_E+size_H));
		Prod_Phi.SetData(this->var_laplace.mesh_num.GetNbDof(), &Y(size_E+size_H));
		Prod_Phi.Fill(0);
		var_leaf.EvaluateE_H_Quadrature(En, Hn, evalEn, evalHn, false, dg_form);
		var_leaf.EvaluateE_H_QuadratureBoundary(En, Hn, evalEnBoundary,
                                                        evalHnBoundary, dg_form);
		var_leaf.EvaluateSourceBoris(Prod_Phi, Jsrc.evalRho, evalEn, evalEnBoundary);
		
		//Real_wp coef(0);
		//if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_NO)
                //coef = square(this->xsi_hyperbolic_correction*PhysicalConstant::speed_light);
		//else
                //coef = square(this->xsi_hyperbolic_correction);
		
		abort();
		// for (int i = 0; i < this->var_laplace.mesh_num.nodl; i++)
		//Prod_Phi(i) = this->var_harmonic.weight_function(i)
                // *(coef*Prod_Phi(i)-this->mass_sigma(i)*Phin(i));
		
		Add(this->deltat, Prod_Phi, Phin);
		Phin.Nullify(); Prod_Phi.Nullify();
	      }
		
	    nb_iter++;
	    t = this->initial_time + this->deltat*nb_iter;
	    
	    En.Nullify(); Hn.Nullify();
	    Prod_En.Nullify(); Prod_Hn.Nullify();
	  }

      }
    else
      TimeMaxwell_2D<TypeElement,TypeEquation>::RunTimeIterations();
    
  }
  
  //! evaluation of operator G, so that time-evolution system reas as dY/dt = G(t,Y)
  /*!
    \param[in] tn time t
    \param[in] nb_deriv n-th derivative of G with respect to t (nb_deriv = n)
    \param[in] X input vector
    \param[out] Y output vector containing the result of G(t,X)
  */
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_2D<TypeElement,TypeEquation>::
  EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv, const VectReal_wp& X,
                             VectReal_wp& Y, bool invert_mass, bool source)
  {
    HyperbolicProblem<TypeElement, TypeEquation>& var_leaf = this->GetLeafClass();
    GhostIf<TypeEquationStationary::FormulationDG> dg_form;
    Real_wp one(1), zero(0);
    // number of dofs for unknown E and H
    int size_E = this->GetNbScalarUnknowns(dg_form);    
    int size_H = this->GetNbVectorialUnknowns(dg_form);
    int offset = this->GetOffset_PosParticle();
    
    // we split the vector X in the electric field part (En) and magnetic field (Hn)
    // same splitting for Y -> Prod_En, Prod_Hn
    VectReal_wp En, Hn, Prod_En, Prod_Hn;
    Real_wp* x_ptr = X.GetData(), *y_ptr = Y.GetData();
    En.SetData(size_E, x_ptr); Prod_En.SetData(size_E, y_ptr);
    x_ptr += size_E; y_ptr += size_E;
    
    Hn.SetData(size_H, x_ptr); Prod_Hn.SetData(size_H, y_ptr);
    x_ptr += size_H; y_ptr += size_H;
    
    // Maxwell equation in E
    this->ApplyOperatorRhScalar(one, tn, Hn, zero, Prod_En);
    if (source)
      this->AddScalarSourceAtTime(one, tn, 0, Prod_En);    
    
    this->ApplyOperatorSh(-one, tn, En, one, Prod_En);
    
    // Maxwell equation in H
    this->ApplyOperatorRhVectorial(one, tn, En, zero, Prod_Hn);
    if (source)
      this->AddVectorialSourceAtTime(one, tn, nb_deriv, Prod_Hn);
    
    this->ApplyOperatorShVectorial(-one, tn, Hn, one, Prod_Hn);
    this->SolveOperatorBh(Prod_Hn);
    
    // treatment of particles
    // first part : localization
    IVect NumBoxGrid_Point;
    GridInterpolation<Dimension2> & grid = this->loc_particle;
    this->LocalizeParticles(X, NumBoxGrid_Point);
        
    // evaluation of E and H on all the quadrature points
    VectR2 evalEn; VectReal_wp evalHn;
    if (this->order_interpolation_current > 0)
      {
	VectR2 En_quad; VectReal_wp Hn_quad;
	this->EvaluateE_H_Quadrature(En, Hn, En_quad, Hn_quad, false, dg_form);
	this->EvaluateE_H_Interpolate(En_quad, Hn_quad, evalEn, evalHn);
      }
    else
      this->EvaluateE_H_Quadrature(En, Hn, evalEn, evalHn, false, dg_form);

    // now we compute J at the quadrature points
    CurrentSource_Vlasov<TypeElement,TypeEquation> Jsrc(var_leaf);
    int Npart = 0; R2 vk, pk, pt_loc; int nb = offset;
    Real_wp poids_J, poids_acc;
    Real_wp val_Ex(0), val_Ey(0), val_H(0);
    VectR2 Jn_quad;
    for (int i = 0; i < this->nb_particles.GetM(); i++)
      {
	poids_acc = this->ratio_qm_species(i); // DISP(i); DISP(poids_acc);
	for (int j = 0; j < this->nb_particles(i); j++)
	  {
	    // int num_elem = grid.ElementInterp(Npart);
	    // pt_loc = grid.CoorInterp(Npart);
	    int num_elem = NumBoxGrid_Point(Npart); // DISP(j); DISP(num_elem);
	    if (num_elem >= 0)
	      {
		// particle inside the domain
		
		// we get the velocity of the particle
		// quantity of movment
		pk(0) = X(nb+4*j+2); pk(1) = X(nb+4*j+3);
		// deducing the velocity
		this->GetVelocity_FromMomentum(pk, vk); 
		// DISP(pk); DISP(vk);
		
		// updating value of J
		poids_J = this->weight_particle(i)(j)*this->charge_species(i);
		Jsrc.UpdateCurrent(poids_J, vk, grid.GlobalCoord(Npart), num_elem, pt_loc);
		
		// computation of the mean value of E and H
		val_Ex = 0; val_Ey = 0; val_H = 0;
		for (int num_point = 0; num_point <
                       this->CloudElem_QuadraturePoint(num_elem).GetM(); num_point++)
		  {
		    Real_wp coefS = Jsrc.CoefficientShape(num_point);
		    if (coefS != Real_wp(0))
		      {
			int num = this->CloudElem_QuadraturePoint(num_elem)(num_point);
			val_Ex += coefS*evalEn(num)(0);
			val_Ey += coefS*evalEn(num)(1);
			val_H += coefS*evalHn(num);
		      }
		    // DISP(evalHn(i_elt,jloc));
		  }
		
		// multypling by mu0 to obtain B and not H
		val_H *= this->var_harmonic.mu0;
		// DISP(Npart); DISP(vk); DISP(X(nb+4*j)); DISP(val_H); DISP(val_Ex); DISP(val_Ey);
		
		// val_H = 10;
		// DISP(val_H); DISP(val_Ex); DISP(val_Ey); DISP(vk); DISP(pk);
		// moving the particle
		// speed part dx_k/dt = v_k
		Y(nb+4*j) = vk(0); Y(nb+4*j+1) = vk(1);
		
		// acceleration part dp_k/dt = q ( E + v \times B)
		Y(nb+4*j+2) = poids_acc*(val_Ex + vk(1)*val_H);
		Y(nb+4*j+3) = poids_acc*(val_Ey - vk(0)*val_H);
		// DISP( poids_acc*(val_Ex + vk(1)*val_H));
	      }
	    
	    Npart++;
	  }
	nb += 4*this->nb_max_particles(i);
      }

    
    // DISP(type_correction_divergence); DISP(CORRECTION_HYPERBOLIC);
    if (this->type_correction_divergence == this->CORRECTION_HYPERBOLIC)
      {
	VectReal_wp Phin, Prod_Phi;
	Phin.SetData(this->var_laplace.mesh_num.GetNbDof(), x_ptr);
        Prod_Phi.SetData(this->var_laplace.mesh_num.GetNbDof(), y_ptr);
	x_ptr += this->var_laplace.mesh_num.GetNbDof();
        y_ptr += this->var_laplace.mesh_num.GetNbDof();
	
	// computation of the source for phi
	Prod_Phi.Fill(0);
	this->EvaluateE_H_Quadrature(En, Hn, evalEn, evalHn, false, dg_form);
	VectReal_wp evalHnBoundary; VectR2 evalEnBoundary;
	this->EvaluateE_H_QuadratureBoundary(En, Hn, evalEnBoundary, evalHnBoundary, dg_form);
	this->EvaluateSourceBoris(Prod_Phi, Jsrc.evalRho, evalEn, evalEnBoundary);
	
	// absorbing boundary condition on phi, and solving of the mass matrix
	//Real_wp coef(0);
	//if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_NO)
        //coef = square(this->xsi_hyperbolic_correction*PhysicalConstant::speed_light);
	//else
        //coef = square(this->xsi_hyperbolic_correction);
	
	abort();
	// for (int i = 0; i < this->var_laplace.mesh_num.nodl; i++)
	//Prod_Phi(i) = this->var_harmonic.weight_function(i)
        // *(coef*Prod_Phi(i)-this->mass_sigma(i)*Phin(i));
	
	// DISP(Norm2(Prod_Phi));
	// computation of grad phi
	this->AddVolumetricIntegral_PhiNodal(-1.0, Phin, Prod_En, dg_form);
	this->AddSurfacicIntegral_PhiNodal(1.0, Phin, Prod_En, dg_form);
	
	Phin.Nullify(); Prod_Phi.Nullify();
      }
    
    // we add current source -J to Prod_En
    if (this->order_interpolation_current > 0)
      {
	var_leaf.EvaluateCurrentQuadrature(Jsrc.evalJ, Jn_quad);
	var_leaf.AddIntegralCurrentJ(one, Jn_quad, Prod_En, dg_form);
      }
    else
      var_leaf.AddIntegralCurrentJ(one, Jsrc.evalJ, Prod_En, dg_form);
    
    // solving mass matrix for E and H
    var_leaf.SolveOperatorDhPlusdtSh(Prod_En);
    
    // for inhomogeneous Dirichlet condition
    this->SetDirichletCondition(tn, 0, Prod_En);
    
    En.Nullify(); Prod_En.Nullify(); Hn.Nullify(); Prod_Hn.Nullify();
  }
  
  
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_2D<TypeElement,TypeEquation>::
  EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv, const Real_wp& alpha, int level,
                             const VectReal_wp& X, VectReal_wp& Y, bool invert_mass, bool source)
  {
    abort();
  }
  
  
  //! localization of particles on the grid
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_2D<TypeElement,TypeEquation>
  ::LocalizeParticles(const VectReal_wp& X, IVect& NumBoxGrid_Point)
  {
    int Npart = this->GetNumberParticles();
    
    // we get the coordinates of all the particles to be localized
    GridInterpolation<Dimension2> & grid = this->loc_particle;
    grid.GlobalCoord.Reallocate(Npart);
    int offset = this->GetOffset_PosParticle(); Npart = 0; int nb = offset;
    for (int i = 0; i < this->nb_particles.GetM(); i++)
      {
	for (int j = 0; j < this->nb_particles(i); j++)
	  {
	    grid.GlobalCoord(Npart)(0) = X(nb+4*j);
	    grid.GlobalCoord(Npart)(1) = X(nb+4*j+1);
	    Npart++;
	  }
	nb += 4*this->nb_max_particles(i);
      }
    
    // if periodicity, points can be shifted
    this->MovePoints_Periodicity(grid.GlobalCoord); // DISP(grid.GlobalCoord);
    Vector<IVect> ListPoints_Grid; IVect NbPoints_Grid;
    grid.PreLocalizePoints(NumBoxGrid_Point, ListPoints_Grid, NbPoints_Grid, false);
    
  }

  
  //! Evaluation of rho
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_2D<TypeElement,TypeEquation>
  ::EvaluateRho(const VectReal_wp& Y, VectReal_wp& bsrc,
                CurrentSource_Vlasov<TypeElement,TypeEquation>& Jsrc)
  {
    Jsrc.evalRho.Fill(0);
    // DISP(*this->var_laplace.FaceBasis);
    // allocation of source
    bsrc.Reallocate(this->var_laplace.GetNbDof()); bsrc.Fill(0);
    
    // evaluation of rho
    IVect NumBoxGrid_Point;
    this->LocalizeParticles(Y, NumBoxGrid_Point);
    
    // computation of rho
    GridInterpolation<Dimension2> & grid = this->loc_particle;
    
    //Real_wp charge;
    int Npart = 0; int nb = 0; R2 vk;
    R2 pt_loc;
    for (int i = 0; i < this->nb_particles.GetM(); i++)
      {
	//charge = this->charge_species(i);
	for (int j = 0; j < this->nb_particles(i); j++)
	  {
	    // int num_elem = grid.ElementInterp(Npart);
	    // pt_loc = grid.CoorInterp(Npart);
	    int num_elem = NumBoxGrid_Point(Npart);
	    if (num_elem >= 0)
	      {
		// particle inside the domain
		// updating value of rho
		Real_wp poids = this->weight_particle(i)(j)*this->charge_species(i);
		Jsrc.UpdateCurrent(poids, vk, grid.GlobalCoord(Npart),
				   num_elem, pt_loc);
	      }
	    Npart++;
	  }
	
	nb += 4*this->nb_max_particles(i);
      }
    
    this->ComputeValue_PhiNodal(bsrc, Jsrc.evalRho);
  }
  
  
  // RUNNING TIME STEPS //
  ////////////////////////
  
  ///////////////////////////
  // DIVERGENCE CORRECTION //
  
  //! writing snapshots if necessary, and position/velocity of the first particle
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_2D<TypeElement, TypeEquation>
  ::WriteSnapshot(int nb_iter, const Real_wp& t, const VectReal_wp& Uh)
  {
    abort();
    //TimeMaxwell_2D<TypeElement, TypeEquation>::WriteSnapshot(nb_iter, t, Uh);
    int offset = this->GetOffset_PosParticle();
    if (nb_iter%10 == 0)
      cout<<"Number of particles "<<this->GetNumberParticles()<<endl;
    
    // if there is at least one particle, we write position and velocity of this particle
    if (this->GetNumberParticles() > 0)
      {
	// DISP(offset); DISP(Uh(offset));
	ofstream file_out(this->name_file_position.data(),ios::app);
	file_out<<Uh(offset)<<" "<<Uh(offset+1)<<"  "<<Uh(offset+2)<<"  "<<Uh(offset+3)<<endl;
	file_out.close();
	
      }
    
    // int test_input; cout<<"we wait"<<endl; cin>>test_input;
  }
  

  template<class TypeElement,class TypeEquation>
  void VlasovMaxwell_2D<TypeElement,TypeEquation>::
  CreateParticlesBreakdownField(const VectReal_wp& En, const VectReal_wp& Hn, int nb_iter,
				VectReal_wp& Y, int offset,
                                VectReal_wp& WeightParticle, VectR2& CoorParticle)
  {
    const Mesh<Dimension2>& mesh = this->var_harmonic.mesh;
    //int order = mesh_num.GetOrder();
    GhostIf<TypeEquationStationary::FormulationDG> dg_form;
    // evaluation of E and H on quadrature points
    VectReal_wp evalHn; VectR2 evalEn;
    this->EvaluateE_H_QuadratureBoundary(En, Hn, evalEn, evalHn, dg_form);

    //int num_elem, num_edge, num_loc, nb_pts_quad_edge, nb_dof_edge;
    //Globatto<Real_wp>& lob1D = this->lob_boundary;
    //list<Real_wp> list_all_weight; list<R2> list_all_point;
    for (int num_edge_ref = 0; num_edge_ref < mesh.GetNbBoundaryRef(); num_edge_ref++)
      for (int num = 0; num < this->breakdown_currents.GetM(); num++)
	if (nb_iter%this->breakdown_currents(num).delta_iter == 0)
	  if (this->breakdown_currents(num).ref == mesh.BoundaryRef(num_edge_ref).GetReference())
	    {
              /*
	      this->var_laplace.GetBoundaryValues(this->var_laplace.mesh, num_edge_ref, num_edge,
						  num_elem, nb_pts_quad_edge,
                                                  nb_dof_edge, num_loc);
	      
	      // charge and mass of particles
	      Real_wp q0 = this->breakdown_currents(num).charge;
	      Real_wp coef_qm = this->breakdown_currents(num).ratio_q0_over_m0;
	      
	      Real_wp En_breakdown = this->breakdown_currents(num).breakdown_field;
	      if (q0 < 0)
		En_breakdown = -En_breakdown;
	      
	      // emitting boundary, checking values of E
	      bool emit = false;
	      VectReal_wp value_En(nb_pts_quad_edge), value_DsEn(nb_pts_quad_edge),
              value_Ds(nb_pts_quad_edge);
	      
              Real_wp int_En(0), dsj(1), val;
	      Real_wp length_edge = 0; R2 valE;
	      VectR2 normale(nb_pts_quad_edge), tangente(nb_pts_quad_edge);
	      VectR2 pointEdge(nb_pts_quad_edge);
	      for (int k = 0; k < nb_dof_edge; k++)
		{
		  valE = evalEn(num_edge_ref*nb_dof_edge + k);
		  normale(k) = this->Glob_normale(num_edge)(k);
		  pointEdge(k) = this->Glob_point_boundary(num_edge)(k);
		  value_DsEn(k) = -DotProd(valE, normale(k));
		  dsj = Norm2(normale(k));
		  value_Ds(k) = dsj;
		  length_edge += dsj*this->WeightsBoundary(k);
		  value_En(k) = value_DsEn(k)/dsj;
		  if (abs(value_En(k)) > abs(En_breakdown))
		    if (q0*value_En(k) > 0)
		      emit = true;
		  
		  Mlt(1.0/dsj, normale(k));
		  tangente(k)(0) = -normale(k)(1);
		  tangente(k)(1) = normale(k)(0);
		}
	      
	      // DISP(value_En); DISP(q0); DISP(emit);
	      
	      if (emit)
		{
		  // creation of particles allowed
		  list<Real_wp> list_weight; list<R2> list_point; R2 pt;
		  VectR2 NewCoor; VectReal_wp NewWeight; Real_wp coef(1);
		  if (this->type_space_charge_limited_emission == this->FIND_INTERVAL)
		    {
		      
		      // computation of the points for which E = En_breakdown
		      // sign_first = 0 if |E| > |En_breakdown|
                      // at the first interval [0,pts(0)], otherwise 1
		      VectReal_wp pts, pts_subdiv; int sign_first = 1;
		      this->ComputePoints_ValEqualConstant(lob1D, value_En, En_breakdown, pts);
		      int_En = 0;
		      for (int j = 0; j < nb_dof_edge; j++)
			int_En += lob1D.EvaluatePhi(j, 0.0)*value_En(j);
		      
		      if (abs(int_En) > abs(En_breakdown))
			sign_first = 0;
		      
		      pts_subdiv.Reallocate(pts.GetM()+2); pts_subdiv(0) = 0;
                      pts_subdiv(pts.GetM()+1) = 1;
		      for (int j = 0; j < pts.GetM(); j++)
			pts_subdiv(j+1) = pts(j);
		      
		      // DISP(pts); DISP(pts_subdiv);
		      // loop over each subdivision for which |E| > breakdown
		      for (int num_subdiv = sign_first;
                      num_subdiv < (pts_subdiv.GetM()-1); num_subdiv += 2)
			{
			  Real_wp t0 = pts_subdiv(num_subdiv), t1 = pts_subdiv(num_subdiv+1);
			  
			  // number of particles to be created
			  Real_wp real_part = this->breakdown_currents(num).nb_part_per_unit
                          *length_edge*(t1-t0);
			  
                          // true number
			  int nb_part = toInteger(round(real_part));
			  // coefficient for the weight
			  coef = 1.0/(nb_part*q0)*PhysicalConstant::epsilon0_permittivity;
			  
			  if (nb_part > 0)
			    {
			      VectReal_wp gamma(nb_part);
			      this->var_random.GenerateRandomNumbers(nb_part, gamma);
			      if (this->breakdown_currents(num).uniform_spatial)
				for (int n = 0; n < nb_part; n++)
				  gamma(n) = (2.0*Real_wp(n)+1)/(2.0*nb_part);
			      
			      for (int n = 0; n < nb_part; n++)
				{
				  // position of the new particle
				  Real_wp lambda = gamma(n)*(t1-t0) + t0;

				  // weight
				  valE(0) = 0; pt.Zero();
				  for (int k = 0; k < nb_dof_edge; k++)
				    {
				      val = lob1D.EvaluatePhi(k, lambda);
				      valE(0) += val*value_DsEn(k);
				      Add(val, pointEdge(k), pt);
				    }
				  
				  list_weight.push_back(coef*(t1-t0)*valE(0));
				  list_point.push_back(pt);
				}
			    }
			}
		    }
		  else if (this->type_space_charge_limited_emission == this->SUBDIVISION_INTERVAL)
		    {
		      // computations of En on regular subdivisions
		      VectReal_wp val_reg_En(order+1); Real_wp mean_value(0);
		      for (int i = 0; i <= order; i++)
			{
			  val_reg_En(i) = 0;
			  Real_wp x = Real_wp(i)/order;
			  for (int j = 0; j < nb_dof_edge; j++)
			    val_reg_En(i) += lob1D.EvaluatePhi(j, x)*value_En(j);
			  
			  mean_value += val_reg_En(i);
			}
		      
		      // DISP(value_En); DISP(val_reg_En);
		      // mean value
		      mean_value /= (order+1); // val_reg_En.Fill(mean_value);
		      // loop on each subdivision
		      for (int i = 0; i < order; i++)
			if (abs(val_reg_En(i) + val_reg_En(i+1)) > 2.0*abs(En_breakdown))
			  if (val_reg_En(i)*En_breakdown > 0)
			    if (val_reg_En(i+1)*En_breakdown > 0)
			      {
				
				// two extremities of the subdivision
				Real_wp t0 = Real_wp(i)/order, t1 = Real_wp(i+1)/order;
				
				// number of particles to be created
				Real_wp real_part
                                = this->breakdown_currents(num).nb_part_per_unit
                                *length_edge*(t1-t0);
				
                                // true number
				int nb_part = toInteger(round(real_part));
				// coefficient for the weight
				coef = 1.0/(nb_part*q0)*PhysicalConstant::epsilon0_permittivity;
				
				if (nb_part > 0)
				  {
				    VectReal_wp gamma(nb_part);
				    this->var_random.GenerateRandomNumbers(nb_part, gamma);
				    if (this->breakdown_currents(num).uniform_spatial)
				      for (int n = 0; n < nb_part; n++)
					gamma(n) = (2.0*Real_wp(n)+1)/(2.0*nb_part);
				    
				    for (int n = 0; n < nb_part; n++)
				      {
					// position of the new particle
					Real_wp lambda = gamma(n)*(t1-t0) + t0;
										
					// weight
					valE(0) = 0; pt.Zero();
					for (int k = 0; k < nb_dof_edge; k++)
					  {
					    val = lob1D.EvaluatePhi(k, lambda);
					    valE(0) += val*value_DsEn(k);
					    Add(val, pointEdge(k), pt);
					  }
					
					list_point.push_back(pt);
					list_weight.push_back(coef*(t1-t0)*valE(0));
				      }
				  }
			      }
		    }
		  else if (this->type_space_charge_limited_emission == this->FULL_INTERVAL)
		    {
		      // computations of En on regular subdivisions
		      VectReal_wp val_reg_En(order+1); Real_wp mean_value(0);
		      for (int i = 0; i <= order; i++)
			{
			  val_reg_En(i) = 0;
			  Real_wp x = Real_wp(i)/order;
			  for (int j = 0; j < nb_dof_edge; j++)
			    val_reg_En(i) += lob1D.EvaluatePhi(j, x)*value_En(j);
			  
			  mean_value += val_reg_En(i);
			}
		      
		      // DISP(value_En); DISP(val_reg_En);
		      // mean value
		      mean_value /= (order+1); val_reg_En.Fill(mean_value);
		      // loop on each subdivision
		      if (mean_value*En_breakdown > 0)
			if (abs(mean_value) > abs(En_breakdown))
			  {
			    
			    // number of particles to be created
			    Real_wp real_part = this->breakdown_currents(num).nb_part_per_unit
                            *length_edge;
			    
                            // true number
			    int nb_part = toInteger(round(real_part));
			    // coefficient for the weight
			    coef = 1.0/(nb_part*q0)*PhysicalConstant::epsilon0_permittivity;
			    
			    if (nb_part > 0)
			      {
				VectReal_wp gamma(nb_part);
				this->var_random.GenerateRandomNumbers(nb_part, gamma);
				if (this->breakdown_currents(num).uniform_spatial)
				  for (int n = 0; n < nb_part; n++)
				    gamma(n) = (2.0*Real_wp(n)+1)/(2.0*nb_part);
				
				for (int n = 0; n < nb_part; n++)
				  {
				    // position of the new particle
				    Real_wp lambda = gamma(n);
								    
				    // weight
				    valE(0) = 0; pt.Zero();
				    for (int k = 0; k < nb_dof_edge; k++)
				      {
					val = lob1D.EvaluatePhi(k, lambda);
					valE(0) += val*value_DsEn(k);
					Add(val, pointEdge(k), pt);
				      }
				    
				    list_point.push_back(pt);
				    list_weight.push_back(coef*valE(0));
				  }
			      }
			  }
		    }
		  
		  Copy(list_point, NewCoor);
		  Copy(list_weight, NewWeight);
		  
		  // creations of particles with randomly small height and velocity
		  if (NewCoor.GetM() > 0)
		    {
		      int nb_part = NewCoor.GetM();
		      
		      Real_wp hmin = this->breakdown_currents(num).height_min;
		      Real_wp hmax = this->breakdown_currents(num).height_max;
		      Real_wp vmin = this->breakdown_currents(num).gamv_normal_min;
		      Real_wp vmax = this->breakdown_currents(num).gamv_normal_max;
		      Real_wp vt_min = this->breakdown_currents(num).gamv_tangential_min;
		      Real_wp vt_max = this->breakdown_currents(num).gamv_tangential_max;
		      VectReal_wp height_particle(nb_part); height_particle.Fill(hmax);
		      VectR2 velocity_particle(nb_part);
		      if (hmin < hmax)
			{
			  this->var_random.GenerateRandomNumbers(nb_part, height_particle);
			  for (int n = 0; n < nb_part; n++)
			    {
			      height_particle(n) = hmin + (hmax-hmin)*height_particle(n); 
			      Add(-height_particle(n), normale(0), NewCoor(n));
			    }
			}
		      
		      if (vmin < vmax)
			{
			  VectReal_wp gamma(nb_part), gammab(nb_part);
			  this->var_random.GenerateRandomNumbers(nb_part, gamma);
			  this->var_random.GenerateRandomNumbers(nb_part, gammab);
			  for (int n = 0; n < nb_part; n++)
			    {
			      Real_wp gamv = (vmin+gamma(n)*(vmax-vmin));
			      Real_wp vt_gamv = (vt_min+gammab(n)*(vt_max-vt_min));
			      velocity_particle(n)(0)
                              = -normale(0)(0)*gamv + tangente(0)(0)*vt_gamv;
			      velocity_particle(n)(1)
                              = -normale(0)(1)*gamv + tangente(0)(1)*vt_gamv;
			    }
			}
		      
		      this->AddParticles(nb_part, q0, coef_qm, velocity_particle,
					 NewCoor, NewWeight, offset, Y);
		      
		      for (int i = 0; i < NewCoor.GetM(); i++)
			{
			  list_all_point.push_back(NewCoor(i));
			  list_all_weight.push_back(NewWeight(i)*q0);
			}
		    }
		}
              */
	    }
    abort();
    //Copy(list_all_point, CoorParticle);
    //Copy(list_all_weight, WeightParticle);
  }
  
  
  //! Treatment of a new iterate given by a time scheme
  /*!
    \param[in] nb_iter number of the iteration
    \param[in] tn time
    \param[in,out] Y the new iterate
  */
  template<class TypeElement,class TypeEquation>
  void VlasovMaxwell_2D<TypeElement,TypeEquation>
  ::GiveIterate(int nb_iter, const Real_wp& tn, VectReal_wp& Y)
  {
    HyperbolicProblem<TypeElement, TypeEquation>& var_leaf = this->GetLeafClass();
    int offset = this->GetOffset_PosParticle(); int nb = 0;
    
    // DISP(type_correction_divergence); DISP(CORRECTION_BORIS); DISP(Norm2(Y));
    
    Mesh<Dimension2>& mesh = this->var_harmonic.mesh;
    MeshNumbering<Dimension2>& mesh_num = this->var_harmonic.mesh_num;
    GhostIf<TypeEquationStationary::FormulationDG> dg_form;
    int size_E = this->GetNbScalarUnknowns(dg_form); 
    int size_H = this->GetNbVectorialUnknowns(dg_form);
    
    VectReal_wp En(size_E), Hn(size_H); R2 pk;
    // Real_wp* x_ptr = Y.GetData();
    // En.SetData(size_E, x_ptr); x_ptr += size_E;
    // Hn.SetData(size_H, x_ptr);
    for (int i = 0; i < size_E; i++)
      En(i) = Y(i);
    
    for (int i = 0; i < size_H; i++)
      Hn(i) = Y(i+size_E);
    
    CurrentSource_Vlasov<TypeElement,TypeEquation> Jsrc(var_leaf);
    
    // for inhomogeneous Dirichlet condition
    this->SetDirichletCondition(tn, 0, Y);
    
    // computing densities of destructed/created particles
    int Npart = this->GetNumberParticles();
    
    // we get the coordinates of all the particles to be localized
    GridInterpolation<Dimension2> & grid = this->loc_particle;
    grid.GlobalCoord.Reallocate(Npart);
    Npart = 0; nb = offset;
    for (int i = 0; i < this->nb_particles.GetM(); i++)
      {
	for (int j = 0; j < this->nb_particles(i); j++)
	  {
	    grid.GlobalCoord(Npart)(0) = Y(nb+4*j);
	    grid.GlobalCoord(Npart)(1) = Y(nb+4*j+1);
	    Npart++;
	  }
	nb += 4*this->nb_max_particles(i);
      }
    // DISP(grid.GlobalCoord);
    
    // particles which have crossed periodic boundaries are moved inside the domain
    this->MovePoints_Periodicity(grid.GlobalCoord);
    Npart = 0; nb = offset; R2 vk;
    for (int i = 0; i < this->nb_particles.GetM(); i++)
      {
	for (int j = 0; j < this->nb_particles(i); j++)
	  {
	    Y(nb+4*j) = grid.GlobalCoord(Npart)(0);
	    Y(nb+4*j+1) = grid.GlobalCoord(Npart)(1);
	    if (this->mirror_right_side)
	      if ((Y(nb+4*j) > this->var_harmonic.GetXmax())&&(Y(nb+4*j+2) > 0))
		{
		  // Y(nb+4*j) = -Y(nb+4*j) + 2.0*this->var_harmonic.xmax;
		  Y(nb+4*j+2) = -Y(nb+4*j+2);
		  grid.GlobalCoord(Npart)(0) = Y(nb+4*j);
		}
		
	    Npart++;
	  }
	nb += 4*this->nb_max_particles(i);
      }
    
    // localization of particles
    // grid.point_treated.Fill(false);
    // grid.LocalizePoints(this->var_harmonic.mesh, Fb, this->var_harmonic.nb_max_vertices);
    IVect NumBoxGrid_Point; Vector<IVect> ListPoints_Grid; IVect NbPoints_Grid;
    grid.PreLocalizePoints(NumBoxGrid_Point, ListPoints_Grid, NbPoints_Grid, false);
    
    for (int n = 0; n < grid.GlobalCoord.GetM(); n++)
      {
	// int num_elem = grid.ElementInterp(n);
	// if (num_elem < 0)
	int num_elem = NumBoxGrid_Point(n);
	// DISP(n); DISP(grid.GlobalCoord(n)); DISP(num_elem);
	// DISP(this->CloudElem_QuadraturePoint(num_elem));
	if ((num_elem < 0) || (this->CloudElem_QuadraturePoint(num_elem).GetM() == 0))
	  {
	    // DISP(grid.ElementInterp(n));
	    // DISP(grid.CoorInterp(n)); DISP(grid.GlobalCoord(n));
	  }
      }
    
    // for each particle outside the domain, we do specific treatment
    // destruction of the particle is the only treatment until now
    nb = offset; Npart = 0;
    for (int i = 0; i < this->nb_particles.GetM(); i++)
      {
	int new_number = 0;
	for (int j = 0; j < this->nb_particles(i); j++)
	  {
	    // int num_elem = grid.ElementInterp(Npart);
	    int num_elem = NumBoxGrid_Point(Npart);
	    if ((num_elem >= 0)&&(this->CloudElem_QuadraturePoint(num_elem).GetM() > 0))
	      {
		this->weight_particle(i)(new_number) = this->weight_particle(i)(j);
		for (int k = 0; k < 4; k++)
		  Y(nb+4*new_number+k) = Y(nb+4*j+k);
		
		new_number++;
		
	      }
	    else
	      {
		// particle will be destructed, updating rho
		// Jsrc.UpdateCurrent(-weight_particle(i)(j), vk,);
	      }
	    Npart++;
	  }
	
	if (new_number < this->nb_particles(i))
	  {
	    // some particles have been removed, setting to 0 removed particles
	    for (int j = new_number; j < this->nb_particles(i); j++)
	      {
		this->weight_particle(i)(j) = 0;
		for (int k = 0; k < 4; k++)
		  Y(nb+4*j+k) = 0;
	      }
	    
	    this->nb_particles(i) = new_number;
	  }
	nb += 4*this->nb_max_particles(i);
      }
    
    // creation of particles if the electric field is larger than breakdown field
    grid.GlobalCoord.Clear();
    VectReal_wp WeighNewParticles;
    this->CreateParticlesBreakdownField(En, Hn, nb_iter, Y, offset,
                                        WeighNewParticles, grid.GlobalCoord);
    
    // creation of particles for current beams
    // DISP(this->Liste_beam_current.GetM());
    // DISP(this->Liste_beam_current(0).delta_iter);
    for (int nbeam = 0; nbeam < this->Liste_beam_current.GetM(); nbeam++)
      if (nb_iter%this->Liste_beam_current(nbeam).delta_iter == 0)
	{
	  // area of the beam ?
	  Real_wp dA = this->Liste_beam_current(nbeam).ptA.
            Distance(this->Liste_beam_current(nbeam).ptB);
	  
	  // intensity
	  Real_wp Jb = this->Liste_beam_current(nbeam).I0;
	  
	  // charge of particles
	  Real_wp q0 = this->Liste_beam_current(nbeam).charge;
	  Real_wp coef_qm = this->Liste_beam_current(nbeam).ratio_q0_over_m0;
	  
	  // velocity of particles multiplied by relativistic factor
	  R2 v0 = this->Liste_beam_current(nbeam).v0; 
	  // now computing quantity of movement
	  R2 p0 = v0; // DISP(v0); DISP(p0);
	  
	  // number of particles to be created
	  Real_wp real_part = this->Liste_beam_current(nbeam).nb_part_per_unit*dA;
	  // true number
	  int nb_part = toInteger(ceil(real_part));
	  // DISP(nb_part);
	  Real_wp poids = dA*Jb*this->deltat*this->Liste_beam_current(nbeam).delta_iter
            /(nb_part*q0);
	  
	  if ((nb_part > 0)&&(poids > 0))
	    {
	      // particles uniformly distributed
	      R2 ptA = this->Liste_beam_current(nbeam).ptA;
	      R2 ptB = this->Liste_beam_current(nbeam).ptB;
	      VectR2 NewCoor(nb_part);  VectReal_wp NewWeight(nb_part); NewWeight.Fill(poids);
	      VectR2 velocity_particle(nb_part); velocity_particle.Fill(p0);
	      for (int n = 0; n < nb_part; n++)
		{
		  Real_wp lambda = (2.0*Real_wp(n)+1)/(2.0*nb_part);
		  Add(1.0-lambda, ptA, NewCoor(n)); Add(lambda, ptB, NewCoor(n));
		  Real_wp eval = this->Liste_beam_current(nbeam).GetFunction(lambda);
		  NewWeight(n) *= eval;
		  // DISP(n); DISP(eval*poids);
		}
	      
	      // if (nb_iter == 0)
	      // {
	      // this->var_deriv->ComputeWeightsBeam(this->Liste_beam_current(nbeam),
              // NewCoor, poids, NewWeight, Jsrc);
	      // this->Liste_beam_current(nbeam).WeightBeam = NewWeight;
	      // DISP(poids); DISP(NewWeight); exit(0);
	      // }
	      // else
	      // NewWeight = this->Liste_beam_current(nbeam).WeightBeam;
	      
	      // DISP(nb_part); DISP(poids); DISP(NewCoor); DISP(NewWeight); DISP(q0);
	      this->AddParticles(nb_part, q0, coef_qm, velocity_particle,
				 NewCoor, NewWeight, offset, Y);
	    }
	}
    
    // if new particles, correction on the electric field
    if (grid.GlobalCoord.GetM() > 0)
      if ((this->type_correction_divergence != this->CORRECTION_BORIS)&&(nb_iter>0))
	{
	  // localization of new particles
	  // grid.LocalizePoints(this->var_laplace.mesh, this->var_laplace.FaceBasis,
          // this->var_laplace.nb_max_vertices);
	  IVect NumBoxGrid_Point; Vector<IVect> ListPoints_Grid; IVect NPoints_Grid;
	  grid.PreLocalizePoints(NumBoxGrid_Point, ListPoints_Grid, NbPoints_Grid, false);
	  
	  // computation of rho
	  VectReal_wp b_rhs(this->var_laplace.GetNbDof()), xsol(this->var_laplace.GetNbDof());
          R2 pt_loc;
	  Jsrc.evalRho.Fill(0); b_rhs.Fill(0);
	  for (int i = 0; i < grid.GlobalCoord.GetM(); i++)
	    {
	      // pt_loc = grid.CoorInterp(i);
	      // num_elem = grid.ElementInterp(i);
	      int num_elem = NumBoxGrid_Point(i);
	      if (num_elem >= 0)
		{
		  // particle inside the domain
		  // updating value of rho
		  Real_wp poids = WeighNewParticles(i);
		  Jsrc.UpdateCurrent(poids, vk, grid.GlobalCoord(i),
				     num_elem, pt_loc);
		}
	    }
	  
	  if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_NO)
	    Mlt(this->var_harmonic.invEpsilon0(0,0), Jsrc.evalRho);
	  
	  // computation of the source (using evalRho and evalEn
	  this->EvaluateSourceBoris(b_rhs, Jsrc.evalRho);
	  
	  // computation of the solution
	  GlobalMatrixH1<Real_wp> nat_mat;
	  this->laplace_solver.ComputeSolution(b_rhs, xsol, nat_mat);
	  
	  // now computing the gradient of phi -> evalEn
	  // updating values of electric field
	  this->AddVolumetricProjection_PhiNodal(-1.0, xsol, Y, dg_form); // DISP(Norm2(Y));
	}
    
    VectReal_wp evalHnBoundary; VectR2 evalEnBoundary;
    this->EvaluateE_H_QuadratureBoundary(En, Hn, evalEnBoundary, evalHnBoundary, dg_form);
    // now correcting the electric field to satisfy div E = rho
    if (this->type_correction_divergence >= 0)
      {
	// we want to correct divergence error
        // ( so that for the new vector, the divergence will be null)
	// if ((type_correction_divergence == CORRECTION_BORIS)&&(nb_iter==1))
	if ((this->type_correction_divergence == this->CORRECTION_BORIS)||(nb_iter==0))
	  {
	    // DISP(this->var_laplace.var_grid);
	    // DISP(this->var_laplace.output_grid_param);
	    // this->var_laplace.InitVarGrid();
	    // this->var_laplace.ComputeVarGrid();
	    
	    // resolution of a Laplacian
	    // computation of the source
	    VectReal_wp b_rhs(this->var_laplace.GetNbDof()), xsol(this->var_laplace.GetNbDof());
	    // evaluation of rho
	    this->EvaluateRho(Y, b_rhs, Jsrc);
	    b_rhs.Fill(0);
	    
	    VectR2 evalEn; VectReal_wp evalHn;
	    this->EvaluateE_H_Quadrature(En, Hn, evalEn, evalHn, false, dg_form);
	    // computation of the source (using evalRho and evalEn
	    this->EvaluateSourceBoris(b_rhs, Jsrc.evalRho, evalEn, evalEnBoundary);
	    // DISP(Norm2(b_rhs));
	    // this->var_laplace.WriteOutputFile(b_rhs, "source_boris");
	    
	    // computation of the solution
	    GlobalMatrixH1<Real_wp> nat_mat;
	    this->laplace_solver.ComputeSolution(b_rhs, xsol, nat_mat);
	    // DISP(Norm2(xsol));
	    // this->var_laplace.WriteOutputFile(xsol, "solution_boris");
	    if (PhysicalConstant::adimensionalization == PhysicalConstant::ADIM_NO)
	      Mlt(this->var_harmonic.invEpsilon0(0,0), xsol);
	    
	    // updating values of electric field
	    VectReal_wp Prod_En(size_E+size_H); Prod_En.Fill(0);
	    this->AddVolumetricProjection_PhiNodal(-1.0, xsol, Prod_En, dg_form);
	    // this->var_deriv->AddVolumetricIntegral_PhiNodal(-1.0, xsol, Prod_En);
	    // this->var_deriv->AddSurfacicIntegral_PhiNodal(1.0, xsol, Prod_En);
	    // this->var_deriv->Solve_MassScalar(Prod_En);
	    // this->var_harmonic.WriteOutputFile(Prod_En, "gradient_boris");
	    for (int i = 0; i < Prod_En.GetM(); i++)
	      Y(i) += Prod_En(i);
	    
	  }
      }
	
    if (this->uniform_magnetic_field)
      this->SetUniform_MagneticField(this->value_uniform_magnetic, Y, dg_form);
     
    if (this->output_rho_param.SnapshotToStore(tn))
      {
	// writing rho
	// evaluation of rho
	VectReal_wp U0_rho; this->EvaluateRho(Y, U0_rho, Jsrc);
	this->output_rho_param.WriteNewSnapshot(nb_iter, tn, U0_rho, this->var_laplace,
                                                this->var_laplace.var_grid);
	
	// evaluation of E
	VectR2 evalEn; VectReal_wp evalHn;
	this->EvaluateE_H_Quadrature(En, Hn, evalEn, evalHn, false, dg_form);
	// computing rho - div (\eps0 E)		
	VectReal_wp b_rhs(this->var_laplace.GetNbDof()); b_rhs.Fill(0);
	// this->var_deriv->EvaluateSourceBoris(b_rhs, Jsrc.evalRho);
	this->EvaluateSourceBoris(b_rhs, Jsrc.evalRho, evalEn, evalEnBoundary);
	
	// inverting by the mass matrix
	for (int i = 0; i < this->var_laplace.GetNbDof(); i++)
	  b_rhs(i) *= this->inv_mass_matrix(i);
	
	this->output_rho_param.WriteSnapshot(nb_iter, tn, b_rhs, this->var_laplace,
                                             this->var_laplace.var_grid, 0);
      }
    
    if (this->computation_energy)
      {		
	VectReal_wp evalHn; VectR2 evalEn;
	this->EvaluateE_H_Quadrature(En, Hn, evalEn, evalHn, true, dg_form);
	// computation of the electromagnetic energy
	Real_wp energy_elec(0), energy_Ex(0), energy_Ey(0), energy_H(0);
	for (int i = 0; i < mesh.GetNbElt(); i++)
	  {
	    int ref = this->var_harmonic.mesh.Element(i).GetReference();
	    Real_wp eps = this->var_harmonic.ref_epsilon(ref)(0,0);
	    Real_wp mu = this->var_harmonic.ref_mu(ref);
	    int N = this->var_harmonic.GetNbPointsQuadratureInside(i);
	    for (int j = 0; j < N; j++)
	      {
		energy_Ex += eps*this->var_harmonic.Glob_jacobian(i)(j)*square(evalEn(i*N+j)(0));
		energy_Ey += eps*this->var_harmonic.Glob_jacobian(i)(j)*square(evalEn(i*N+j)(1));
		energy_H += mu*this->var_harmonic.Glob_jacobian(i)(j)*evalHn(i*N+j)*evalHn(i*N+j);
	      }
	  }
	
	energy_elec = energy_Ex+energy_Ey+energy_H;
	
	// computation of the cinetic energy
	Real_wp energy_cinetic(0); nb = offset;
	Real_wp c02, invC02, gamma;
	c02 = PhysicalConstant::speed_light; c02 *= c02; invC02 = 1.0/c02;
	for (int i = 0; i < this->nb_particles.GetM(); i++)
	  {
	    Real_wp mass = this->charge_species(i)/this->ratio_qm_species(i);
	    for (int j = 0; j < this->nb_particles(i); j++)
	      {
		// we get the velocity of the particle
		// quantity of movment
		pk(0) = Y(nb+4*j+2); pk(1) = Y(nb+4*j+3);
		
		// incrementing cinetic energy
		if (this->relativist_particle)
		  {
		    gamma = sqrt(1.0 + DotProd(pk,pk)*invC02);
		    energy_cinetic += mass*c02*(gamma-1.0)*this->weight_particle(i)(j);
		  }
		else
		  energy_cinetic += 0.5*mass*this->weight_particle(i)(j)*DotProd(pk,pk);
		
	      }
	    nb += 4*this->nb_max_particles(i);
	  }
	
	ofstream file_out(this->file_output_energy.data(), ios::app);
	file_out<<tn<<" "<<energy_Ex<<" "<<energy_Ey<<" "<<energy_H<<" "<<energy_elec<<" "
                <<energy_cinetic<<" "<<energy_cinetic+energy_elec<<endl;
	file_out.close();
	cout<<"Energy equal to "<<energy_cinetic+energy_elec<<endl;
	
      }
    
    if (this->output_particle_param.SnapshotToStore(tn))
      {
	nb = offset; Npart = 0;
	if (this->nb_particles.GetM() > this->output_particle_param.GetNature())
	  {
	    int nb_part_species = this->nb_particles(this->output_particle_param.GetNature());
	    int max_part = min(this->output_particle_param.GetComponent(), nb_part_species);
	    Matrix<float> positions(max_part,4); positions.Fill(0);
	    for (int i = 0; i < this->output_particle_param.GetNature(); i++)
	      nb += 4*this->nb_max_particles(i);
	    
	    for (int j = 0; j < nb_part_species; j++)
	      if (Npart < max_part)
		{
		  positions(Npart,0) = Y(nb+4*j);
		  positions(Npart,1) = Y(nb+4*j+1);
		  positions(Npart,2) = Y(nb+4*j+2);
		  positions(Npart,3) = Y(nb+4*j+3);
		  Npart++;
		}
	    
	    string numero = to_str(this->output_particle_param.IncrementSnapshot());
	    string file_name = GetBaseString(this->output_particle_param.GetFileName(1))
              + "_T" + numero + ".dat"; 
	    positions.Write(file_name);
	  }
      }
       
    Vector<bool> compute_snap(this->var_harmonic.output_grid_param.GetM());
    compute_snap.Fill(false);
    for (int i = 0; i < this->var_harmonic.output_grid_param.GetM(); i++)
      compute_snap(i) = this->var_harmonic.output_grid_param(i).SnapshotToStore(tn);
	
    this->WriteSnapshot(nb_iter, tn, Y);
   
    for (int i = 0; i < this->var_harmonic.output_grid_param.GetM(); i++)
      if (compute_snap(i))
	{
	  VectReal_wp b_rhs(this->var_laplace.GetNbDof()); b_rhs.Fill(0);
	  int order = mesh_num.GetOrder();
	  int Nh = (order+1)*(order+1);
	  for (int n = 0; n < mesh.GetNbElt(); n++)
	    for (int j = 0; j < Nh; j++)
	      b_rhs(this->var_laplace.mesh_num.Element(n).GetNumberDof(j)) = Hn(n*Nh + j);
	  
	  this->var_harmonic.output_grid_param(i).
	    WriteSnapshot(nb_iter, tn, b_rhs, this->var_laplace, this->var_laplace.var_grid, 2);
	}
    // DISP(this->var_harmonic.compteur_temps.GetSecondsChrono(this->var_harmonic.chrono_rigid));
    // DISP(this->var_harmonic.compteur_temps.GetSecondsChrono(this->var_harmonic.chrono_flux));
    // DISP(this->var_harmonic.compteur_temps
    // .GetSecondsChrono(this->var_harmonic.chrono_prod_mat_vect));
    // DISP(this->var_harmonic.compteur_temps.GetSecondsChrono(this->var_harmonic.chrono_mass));
    // DISP(nb_particles); 
  }
  
  
  //! adding particles in the simulation
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_2D<TypeElement, TypeEquation>
  ::AddParticles(int nb_part, Real_wp& q0, Real_wp& ratio_q0_m0,
                 VectR2& p0, VectR2& NewCoor, VectReal_wp& poids, int offset, VectReal_wp& Y)
  {
    int nb_p = 0, nb; R2 point;
    // DISP(nb_part); DISP(poids); DISP(v0);
    // DISP(q0); DISP(m0); DISP(Jb); DISP(dA); DISP(nb_particles);
    while (nb_part > 0)
      {
	// loop on all the species, to find largest available space
	int num_species = -1, size_available = 0, nb_slot = 0; nb = offset;
	for (int i = 0; i < this->nb_particles.GetM(); i++)
	  {
	    if (abs(this->charge_species(i)-q0) <= epsilon_machine)
	      if (abs(this->ratio_qm_species(i)-ratio_q0_m0) <= epsilon_machine)
		{
		  int size_slot = this->nb_max_particles(i) - this->nb_particles(i);
		  if (size_slot > size_available)
		    {
		      num_species = i;
		      size_available = size_slot;
		      nb_slot = nb;
		    }
		}
	    nb += 4*this->nb_max_particles(i);
	  }
	
	// if no species found or no available space, we create a new species
	if (num_species == -1)
	  {
	    num_species = this->nb_particles.GetM();
	    this->CreateNewSpecies(q0, ratio_q0_m0, Y);
	    size_available = this->nb_max_particles(num_species);
	    nb_slot = nb;
	  }
	
	// initializing the particles
	for (int j = 0; j < min(size_available, nb_part); j++)
	  {
	    // number of the particle
	    int num = this->nb_particles(num_species) + j;
	    
	    // position
	    point = NewCoor(nb_p);
	    
	    Y(nb_slot + 4*num) = point(0);
	    Y(nb_slot + 4*num+1) = point(1);
		
	    // same velocity for all particles of the beam
	    Y(nb_slot + 4*num + 2) = p0(nb_p)(0);
	    Y(nb_slot + 4*num + 3) = p0(nb_p)(1);
	    
	    // weight
	    this->weight_particle(num_species)(num) = poids(nb_p);
	    
	    nb_p++;
	  }
	
	// less particles to create now ...
	this->nb_particles(num_species) += min(size_available,nb_part);
	nb_part -= min(size_available,nb_part);
      }
  }
  
  
  //! computation of weights to fit the beam shape
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_2D<TypeElement, TypeEquation>
  ::ComputeWeightsBeam(Vlasov_BeamCurrent<Dimension2>& current, VectR2& NewCoor,
                       Real_wp& poids, VectReal_wp& NewWeight,
                       CurrentSource_Vlasov<TypeElement,TypeEquation>& Jsrc)
  {
    // "Galerkin" solution to compute weights of the beam
    int N = NewCoor.GetM();
    
    R2 pt; VectReal_wp b_rhs(N); b_rhs.Fill(0);
    // All_MatrixLU mat_lu; Matrix<Real_wp, General, ArrayRowSparse> A(N,N);
    Real_wp dist, length;
    length = current.ptA.Distance(current.ptB);
    int nb_interv = max(1, toInteger(length/Jsrc.cut_off_Sshape));
    Real_wp L0 = -Jsrc.cut_off_Sshape/length, L1 = 1.0 - L0, lambda;
    Real_wp x0, x1; VectReal_wp evalP(N); VectBool interac_null(N);
    Globatto<Real_wp>& lob = this->var_harmonic.FaceBasis->lob1D;
    // loop over quadrature point
    for (int i = 0; i < (nb_interv+2); i++)
      {
	if (i == 0)
	  {
	    x0 = L0;
	    x1 = 0.0;
	  }
	else if (i == nb_interv+1)
	  {
	    x0 = 1.0;
	    x1 = L1;
	  }
	else
	  {
	    x0 = Real_wp(i-1)/nb_interv;
	    x1 = Real_wp(i)/nb_interv;
	  }
	
	for (int k = 0; k <= lob.order; k++)
	  {
	    // position of the quadrature point
	    lambda = (1.0-lob.Points(k))*x0 + lob.Points(k)*x1;
	    pt.Zero(); Add(lambda, current.ptB, pt); Add(1.0-lambda, current.ptA, pt);
	    
	    // weight
	    Real_wp weight = (x1-x0)*length*lob.Weights(k);
	    
	    // evaluation of the beam distribution
	    Real_wp eval = current.GetFunction(lambda);
	    
	    // evalution of basis functions on this quadrature point
	    // and modification of right handside
	    evalP.Fill(0); interac_null.Fill(true);
	    for (int p = 0; p < N; p++)
	      {
		dist = NewCoor(p).Distance(pt);
		if (dist < Jsrc.cut_off_Sshape)
		  {
		    interac_null(p) = false;
		    evalP(p) = Jsrc.GetShapeFunction(dist);
		    b_rhs(p) += weight*evalP(p)*eval;
		  }
	      }
	    
	    // double loop over "basis" functions
	    //for (int p = 0; p < N; p++)
	    // if (!interac_null(p))
	    // for (int q = 0; q < N; q++)
	    // if (!interac_null(q))
	    // A.AddInteraction(p, q, weight*evalP(p)*evalP(q));
	    
	  }
      }
    // WriteMatrix_Matlab("Ah.dat", A); b_rhs.WriteText("rhs.dat");
    // DISP(b_rhs); DISP(A); DISP(NewCoor);
    
    // resolution
    // GetLU(A, mat_lu, Direct_Mumps);
    // SolveLU(mat_lu, b_rhs);
    // DISP(b_rhs); b_rhs.WriteText("sol.dat");
    
    // now deducing the weights
    Real_wp sum(0);
    for (int i = 0; i < N; i++)
      sum += b_rhs(i);
    
    sum /= N;
    
    for (int i = 0; i < N; i++)
      NewWeight(i) = poids/sum*b_rhs(i);
    
    // DISP(poids); DISP(NewWeight);
  }
  

  //! shifting points so that they belong to computational domain
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_2D<TypeElement,TypeEquation>::MovePoints_Periodicity(VectR2& Points)
  {
    if (this->periodic_vector.GetM() > 0)
      {
	Real_wp xmin = this->var_harmonic.GetXmin();
	Real_wp xmax = this->var_harmonic.GetXmax();
	Real_wp ymin = this->var_harmonic.GetYmin();
	Real_wp ymax = this->var_harmonic.GetYmax(); R2 point;
	for (int i = 0; i < Points.GetM(); i++)
	  if ((Points(i)(0) < xmin)||(Points(i)(0) > xmax)||(Points(i)(1) < ymin)
              ||(Points(i)(1) > ymax))
	    {
	      // DISP(i); DISP(Points(i));
	      for (int m = 1; m < this->periodic_vector.GetM(); m++)
		{
		  Add(Points(i), this->periodic_vector(m), point); // DISP(m); DISP(point);
		  if ((point(0) >= xmin)&&(point(0) <= xmax)&&(point(1) >= ymin)
                      &&(point(1) <= ymax))
		    {
		      Points(i) = point; // cout<<" point moved "<<point<<endl;
		      break;
		    }
		}
	    }
      }
    
  }
  
  // DIVERGENCE CORRECTION //
  ///////////////////////////
  
  ///////////////////////////
  // CREATION OF PARTICLES //
  
  
  //! initial fields and position particles (preload of particles)
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_2D<TypeElement,TypeEquation>
  ::SetInitialVector(const Real_wp& t_begin, VectReal_wp& Y)
  {
    int offset = this->GetOffset_PosParticle();
    // cout<<" on passe la "<<endl;
    // initialization for the electric field
    VarInstationary<TypeElement,TypeEquation>::SetInitialVector(t_begin, Y);
    
    // setting electric field
    R2 x0, v0, p0; Real_wp coef_qm, charge;
    VectReal_wp Rx, Ry, Rnu, Rteta; int Npart = 0;
    for (int num_p = 0; num_p < this->initial_particles_param.GetM(); num_p++)
      {
	VectString& parameters = this->initial_particles_param(num_p);
	int nb = 1; // DISP(parameters(0));
	if (!parameters(0).compare("ONE"))
	  {
	    // only one particle, we read the mass, charge, initial position and speed
	    this->nb_particles.PushBack(1);
	    this->nb_max_particles.PushBack(this->increment_particles);
	    VectReal_wp weight_zero(this->increment_particles); weight_zero.Fill(0);
	    Vlasov_BeamCurrent<Dimension2>::SetInputSpecies(charge, coef_qm, parameters, nb); 
	    this->charge_species.PushBack(charge); this->ratio_qm_species.PushBack(coef_qm);
	    x0(0) = to_num<Real_wp>(parameters(nb++)); x0(1) = to_num<Real_wp>(parameters(nb++));
	    // v0 is the velocity multiplied by the relativistic factor
	    v0(0) = to_num<Real_wp>(parameters(nb++)); v0(1) = to_num<Real_wp>(parameters(nb++));
	    p0 = v0;
	    // weight
	    weight_zero(0) = to_num<Real_wp>(parameters(nb++));
	    this->weight_particle.PushBack(weight_zero);
	    
	    Y(offset) = x0(0); Y(offset+1) = x0(1);
	    Y(offset+2) = p0(0); Y(offset+3) = p0(1);
	    
	    offset += 4*this->increment_particles;
	  }
	else if (!parameters(0).compare("CIRCLE"))
	  {
	    int n = to_num<int>(parameters(nb++));
	    if (n > this->increment_particles)
	      {
		cout<<"Too many particles to add "<<endl;
		abort();
	      }
	    R2 center;
	    center(0) = to_num<Real_wp>(parameters(nb++));
	    center(1) = to_num<Real_wp>(parameters(nb++));
	    Real_wp radius = to_num<Real_wp>(parameters(nb++));
	    this->nb_particles.PushBack(n);
	    this->nb_max_particles.PushBack(this->increment_particles);
	    VectReal_wp weight_zero(this->increment_particles); weight_zero.Fill(0);
	    for (int i = 0; i < n; i++)
	      weight_zero(i) = 1.0;
	    
	    this->weight_particle.PushBack(weight_zero);
	    Vlasov_BeamCurrent<Dimension2>::SetInputSpecies(charge, coef_qm, parameters, nb);
	    // speed is the velocity multiplied by the relativistic factor
	    Real_wp speed = to_num<Real_wp>(parameters(nb++));
	    this->charge_species.PushBack(charge); this->ratio_qm_species.PushBack(coef_qm);
	    
	    for (int i = 0; i < n; i++)
	      {
		Real_wp teta = Real_wp(i)/n*2.0*pi_wp;
		Y(offset+4*i) = center(0)+radius*cos(teta);
		Y(offset+4*i+1) = center(1)+radius*sin(teta);
		v0(0) = speed*cos(teta);
		v0(1) = speed*sin(teta);
		Y(offset+4*i+2) = v0(0);
		Y(offset+4*i+3) = v0(1);
	      }
	    
	    offset += 4*this->increment_particles;
	  }
	else if (!parameters(0).compare("PLASMA"))
	  {
	    this->var_random.SetInputRandomGenerator(parameters(nb)); nb++;
            int type_position = this->var_random.type_generator;
	    this->var_random.SetInputRandomGenerator(parameters(nb)); nb++;
            int type_velocity = this->var_random.type_generator;
	    
	    int Ntot = to_num<int>(parameters(nb++));
	    // vth is the velocity multiplied by the relativistic factor
	    Real_wp vth_x = to_num<Real_wp>(parameters(nb++));
	    Real_wp vth_y = to_num<Real_wp>(parameters(nb++));
	    Real_wp alpha = to_num<Real_wp>(parameters(nb++));
	    Real_wp kx = to_num<Real_wp>(parameters(nb++));
	    Vlasov_BeamCurrent<Dimension2>::SetInputSpecies(charge, coef_qm, parameters, nb);
	    Real_wp rho0 = to_num<Real_wp>(parameters(nb++));
	    Real_wp xmin_ = to_num<Real_wp>(parameters(nb++));
	    Real_wp xmax_ = to_num<Real_wp>(parameters(nb++));
	    Real_wp ymin_ = to_num<Real_wp>(parameters(nb++));
	    Real_wp ymax_ = to_num<Real_wp>(parameters(nb++));
	    
	    int nx(0), ny(0);
	    
	    Real_wp alphab = (ymax_-ymin_)/(xmax_-xmin_);
	    VectReal_wp Rx_glob, Ry_glob, Rnu_glob, Rteta_glob;
	    if (type_position == VarRandomGenerator::UNIFORM)
	      {
		this->var_random.RoundToSquare(Ntot, alphab, nx, ny);
		Rx_glob.Reallocate(Ntot); Ry_glob.Reallocate(Ntot);
		
		Real_wp dx = Real_wp(1)/nx, x0 = 0.5*dx;
		Real_wp dy = Real_wp(1)/ny, y0 = 0.5*dy;
		for (int i = 0; i < nx; i++)
		  for (int j = 0; j < ny; j++)
		    {
		      Rx_glob(i*ny+j) = x0 + i*dx;
		      Ry_glob(i*ny+j) = y0 + j*dy;
		    }
		// DISP(dx); DISP(dy); DISP(x0); DISP(y0); DISP(Rx_glob); DISP(Ry_glob);
	      }
	    else if (type_position == VarRandomGenerator::BIT_REVERSED)
	      {
		this->var_random.GenerateBitReversedNumbers(Ntot, Rx_glob);
		this->var_random.GenerateBitReversedNumbers(Ntot, Ry_glob);
	      }
		    
	    if (type_velocity == VarRandomGenerator::UNIFORM)
	      {
		Rnu_glob.Reallocate(Ntot); Rteta_glob.Reallocate(Ntot);
		
		Real_wp dx = Real_wp(1)/nx, x0 = 0.5*dx;
		Real_wp dy = Real_wp(1)/ny, y0 = 0.5*dy;
		for (int i = 0; i < nx; i++)
		  for (int j = 0; j < ny; j++)
		    {
		      Rnu_glob(i*ny+j) = x0 + i*dx;
		      Rteta_glob(i*ny+j) = y0 + j*dy;
		    }
		
		this->var_random.ApplyRandomPermutation(Rnu_glob);
		this->var_random.ApplyRandomPermutation(Rteta_glob);
	      }
	    else if (type_velocity == VarRandomGenerator::BIT_REVERSED)
	      {
		this->var_random.GenerateBitReversedNumbers(Ntot, Rnu_glob);
		this->var_random.GenerateBitReversedNumbers(Ntot, Rteta_glob);
	      }
	    
	    // multipication of rho0 by the area of the domain
	    rho0 *= (xmax_-xmin_)*(ymax_-ymin_);
	    
	    int nb_div = (Ntot-1)/this->increment_particles + 1, n = 0;
	    Npart = 0; int Npart_speed = 0;
	    for (int k = 0; k < nb_div; k++)
	      {
		if (k < (nb_div-1))
		  n = this->increment_particles;
		else
		  n = Ntot-(nb_div-1)*this->increment_particles;
		
		if (type_position == VarRandomGenerator::RANDOM)
		  {
		    this->var_random.GenerateRandomNumbers(n, Rx);
		    this->var_random.GenerateRandomNumbers(n, Ry);
		  }
		else
		  {
		    Rx.Reallocate(n); Ry.Reallocate(n);
		    for (int i = 0; i < n; i++)
		      {
			Rx(i) = Rx_glob(Npart);
			Ry(i) = Ry_glob(Npart);
			Npart++;
		      }
		  }
		
		if (type_velocity == VarRandomGenerator::RANDOM)
		  {
		    this->var_random.GenerateRandomNumbers(n, Rnu);
		    this->var_random.GenerateRandomNumbers(n, Rteta);
		  }
		else
		  {
		    Rnu.Reallocate(n); Rteta.Reallocate(n);
		    for (int i = 0; i < n; i++)
		      {
			Rnu(i) = Rnu_glob(Npart_speed);
			Rteta(i) = Rteta_glob(Npart_speed);
			Npart_speed++;
		      }
		  }
		
		this->nb_particles.PushBack(n);
		this->nb_max_particles.PushBack(this->increment_particles);
		VectReal_wp weight_zero(this->increment_particles); weight_zero.Fill(0);
		for (int i = 0; i < n; i++)
		  weight_zero(i) = rho0*(xmax_-xmin_)*(ymax_-ymin_)/Ntot;
		
		this->weight_particle.PushBack(weight_zero);
		
		this->charge_species.PushBack(charge); this->ratio_qm_species.PushBack(coef_qm);
		
		Real_wp coef = alpha/kx;
		for (int i = 0; i < n; i++)
		  {
		    Real_wp nu = sqrt(-2.0*log(1-Rnu(i)));
		    Real_wp teta = 2.0*pi_wp*Rteta(i);
		    // inverting distribution function with a Newton method
		    if (abs(alpha) > 1e-10)
		      {
			Real_wp val = (xmax_-xmin_)*Rx(i), test_prec = 2e20, test = 1e20;
			x0(0) = val; Real_wp  feval(0), fderiv(0);
			while (test < test_prec)
			  {
			    test_prec = test;
			    feval = x0(0) + coef*sin(kx*x0(0))-val;
			    test = abs(feval);
			    fderiv = 1.0 + alpha*cos(kx*x0(0));
			    x0(0) -= feval/fderiv;
			  }
		      }
		    else
		      x0(0) = (xmax_-xmin_)*Rx(i);
		    
		    Y(offset+4*i) = xmin_+x0(0);
		    Y(offset+4*i+1) = ymin_+(ymax_-ymin_)*Ry(i);
		    v0(0) = nu*vth_x*cos(teta);
		    v0(1) = nu*vth_y*sin(teta);
		    // DISP(v0); DISP(xmin_); DISP(ymin_);
                    // DISP(Y(offset+4*i)); DISP(Y(offset+4*i+1));
		    Y(offset+4*i+2) = v0(0);
		    Y(offset+4*i+3) = v0(1);
		  }
		
		offset += 4*this->increment_particles;
	      }
	  }
      }
  }
  
  // CREATION OF PARTICLES //
  ///////////////////////////
  
}

#define MONTJOIE_FILE_VLASOV_MAXWELL2D_CXX
#endif
