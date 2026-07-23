#ifndef MONTJOIE_FILE_VLASOV_MAXWELL3D_CXX

#include "VlasovMaxwell3D.hxx"

namespace Montjoie
{
  ////////////////////
  // INITIALIZATION //
  
  //! constructor for 3-D Vlasov-Maxwell problem
  template<class TypeElement, class TypeEquation>
  VlasovMaxwell_3D<TypeElement,TypeEquation>::VlasovMaxwell_3D()
    : TimeMaxwell_3D<TypeElement,TypeEquation>(),
      VlasovMaxwell_Base<TypeElement, TypeEquation>(this->GetLeafClass())
  {
    uniform_magnetic_field = false;
    value_uniform_magnetic.Fill(0); 
  }
  
  
  //! Reading of an input file line
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_3D<TypeElement,TypeEquation>::
  SetInputData(const string& description_field, const VectString& parameters)
  {
    TimeMaxwell_3D<TypeElement,TypeEquation>::SetInputData(description_field, parameters);
    VlasovMaxwell_Base<TypeElement,TypeEquation>::SetInputData(description_field, parameters);
    
    if (!description_field.compare("UniformMagneticField"))
      {
	for (int k = 0; k < 3; k++)
	  this->value_uniform_magnetic(k) = to_num<Real_wp>(parameters(k));
	
	this->uniform_magnetic_field = true;
      }
  }
  
  
  template<class TypeElement,class TypeEquation>
  void VlasovMaxwell_3D<TypeElement,TypeEquation>::InitTimeIterations()
  {
    TimeMaxwell_3D<TypeElement, TypeEquation>::InitTimeIterations();
    VlasovMaxwell_Base<TypeElement, TypeEquation>::InitTimeIterations();
  }
  
  
  template<class TypeElement,class TypeEquation>
  void VlasovMaxwell_3D<TypeElement,TypeEquation>::
  ComputeProjectors(VectR3& PtsInterpol3D, VectR3& PtsQuad3D, VectReal_wp& Weights3D)
  {
    //const Mesh<Dimension3>& mesh = this->var_harmonic.mesh;
    const MeshNumbering<Dimension3>& mesh_num = this->var_harmonic.mesh_num;
    int r = mesh_num.GetOrder();
    //const TypeElement& Fb = this->var_harmonic.GetFaceBasis(r);
    
    // constructing quadrature points
    Globatto<Real_wp> gauss;
    // DISP(this->order_integration_current);
    gauss.ConstructQuadrature(this->order_integration_current);
    VectReal_wp PtsQuad1D = gauss.Points();
    HexahedronQuadrature::
      ConstructQuadrature(2*this->order_integration_current, PtsQuad3D, Weights3D);
    
    // we assume that dofs are defined on Gauss-Lobatto points
    // (edge finite element or discontinuous Galerkin)
    // const Array3D<int>& NumNodes3D = Fb.GetNumNodes3D();

    abort();
    //this->proj_nodal_to_quadrature.Init(phi_dof, NumNodes3D, PtsQuad1D, PtsQuad3D);
    
    VectR2 PtsQuad2D;
    QuadrangleQuadrature::ConstructQuadrature(2*r, PtsQuad2D, WeightsBoundary);    
    
    if (this->order_interpolation_current > 0)
      {
	// interpolation points are defined on regular points (or lobatto)
	// with a given number of subdivisions
	SubdivGlobatto phi_interp;
	phi_interp.Init(this->regular_interpolation_current,
                        this->nb_subdiv_interpolation, this->order_interpolation_current);
	int ri = phi_interp.GetOrder();
	PtsInterpol3D.Reallocate((ri+1)*(ri+1)*(ri+1));
	Array3D<int> NumNodes3D_int(ri+1, ri+1, ri+1);
	for (int i = 0; i <= ri; i++)
	  for (int j = 0; j <= ri; j++)
	    for (int k = 0; k <= ri; k++)
	      {
		NumNodes3D_int(i, j, k) = k + j*(ri+1) + i*(ri+1)*(ri+1);
		PtsInterpol3D(NumNodes3D_int(i, j, k))
		  .Init(phi_interp.Points(i), phi_interp.Points(j), phi_interp.Points(k));
	      }
	
	this->proj_interpolate_to_quadrature.Init(phi_interp, NumNodes3D_int,
                                                  PtsQuad1D, PtsQuad3D);
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
    \param[out] Cloud_pt list of interpolation points (local number in the element)
    to scan for each subdivision of the grid
    \param[out] grid interpolation grid
   */
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_3D<TypeElement,TypeEquation>::
  PreComputeCloudsGrid(VectR3& pts, const Real_wp& radius,
                       Vector<IVect>& Cloud_elem, GridInterpolation<Dimension3>& grid)
  {
    Mesh<Dimension3>& mesh = this->var_harmonic.mesh;
    MeshNumbering<Dimension3>& mesh_num = this->var_harmonic.mesh_num;
    GeneratePeriodicDirections();

    // is there a periodicity condition for x or y axis ?
    this->GeneratePeriodicDirections();

    // is there a periodicity condition for x or y axis ?
    bool periodic_x = false, periodic_y = false, periodic_z = false;
    for (int num = 0; num < mesh.GetNbPeriodicReferences(); num++)
      {
	if ((abs(mesh_num.GetTranslationPeriodicBoundary(num)(1)) < 1e-5)
	    &&(abs(mesh_num.GetTranslationPeriodicBoundary(num)(2)) < 1e-5))
	  periodic_x = true;
	
	if ((abs(mesh_num.GetTranslationPeriodicBoundary(num)(0)) < 1e-5)
	    &&(abs(mesh_num.GetTranslationPeriodicBoundary(num)(2)) < 1e-5))
	  periodic_y = true;
	
	if ((abs(mesh_num.GetTranslationPeriodicBoundary(num)(0)) < 1e-5)
	    &&(abs(mesh_num.GetTranslationPeriodicBoundary(num)(1)) < 1e-5))
	  periodic_z = true;
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
    
    if (!periodic_z)
      {
	grid.SetZmin(mesh.GetZmin() - 1.5*radius);
	grid.SetZmax(mesh.GetZmax() + 1.5*radius);
      }
    else
      {
	grid.SetZmin(mesh.GetZmin());
	grid.SetZmax(mesh.GetZmax());
      }
    
    int nbx = max(toInteger(ceil( (grid.GetXmax()-grid.GetXmin())/delta )), 1);
    int nby = max(toInteger(ceil( (grid.GetYmax()-grid.GetYmin())/delta )), 1);
    int nbz = max(toInteger(ceil( (grid.GetZmax()-grid.GetZmin())/delta )), 1);
    Real_wp step_x(0), step_y(0), step_z(0);
    step_x = (grid.GetXmax()-grid.GetXmin())/nbx;
    step_y = (grid.GetYmax()-grid.GetYmin())/nby;
    step_z = (grid.GetZmax()-grid.GetZmin())/nbz;
    grid.SetNbSubdivisions(nbx, nby, nbz);
    grid.SetSubdivisionStep(step_x, step_y, step_z);
    
    // second step : all the quadrature points are localized on this grid
    grid.GlobalCoord = pts;
    
    // NumBoxGrid_Point(i) : number of the small box of the grid where the point i is
    // ListPoints_Grid(i) : list of points contained by the small box i
    // NbPoints_Grid(i) : amount of points contained by the small box i
    IVect NumBoxGrid_Point; Vector<IVect> ListPoints_Grid; IVect NbPoints_Grid;
    grid.PreLocalizePoints(NumBoxGrid_Point, ListPoints_Grid, NbPoints_Grid, true);
    
    // third step  for each small box of the regular grid,
    // we look for neighbouring quadrature points
    VectR3 pts_box(8); Real_wp dist; R3 point, pt_glob;
    Cloud_elem.Reallocate(nbx*nby*nbz);
    int imin(0), imax(0), jmin(0), jmax(0), kmin(0), kmax(0);
    // DISP(step_x); DISP(step_y); DISP(nbx); DISP(nby);
    // DISP(NbPoints_Grid);
    for (int ix = 0; ix < nbx; ix++)
      for (int iy = 0; iy < nby; iy++)
	for (int iz = 0; iz < nbz; iz++)
	  {
	    int num_box = iz*nbx*nby + iy*nbx + ix;
	    int size_res = 1000;
	    Cloud_elem(num_box).Reallocate(size_res);
	  
	    // eight vertices of the box
	    pts_box(0)(0) = grid.GetXmin() + ix*step_x;
            pts_box(0)(1) = grid.GetYmin() + iy*step_y;
            pts_box(0)(2) = grid.GetZmin() + iz*step_z;
	    
            pts_box(1)(0) = pts_box(0)(0) + step_x;
            pts_box(1)(1) = pts_box(0)(1);
            pts_box(1)(2) = pts_box(0)(2); 
	    
            pts_box(2)(0) = pts_box(0)(0) + step_x;
            pts_box(2)(1) = pts_box(0)(1)+step_y; pts_box(2)(2) = pts_box(0)(2); 
	    pts_box(3)(0) = pts_box(0)(0);
            pts_box(3)(1) = pts_box(0)(1)+step_y; pts_box(3)(2) = pts_box(0)(2);
	    pts_box(4)(0) = pts_box(0)(0);
            pts_box(4)(1) = pts_box(0)(1); pts_box(4)(2) = pts_box(0)(2) + step_z; 
	    pts_box(5)(0) = pts_box(4)(0) + step_x;
            pts_box(5)(1) = pts_box(4)(1); pts_box(5)(2) = pts_box(4)(2); 
	    pts_box(6)(0) = pts_box(4)(0) + step_x;
            pts_box(6)(1) = pts_box(4)(1)+step_y; pts_box(6)(2) = pts_box(4)(2); 
	    pts_box(7)(0) = pts_box(4)(0);          
            pts_box(7)(1) = pts_box(4)(1)+step_y; pts_box(7)(2) = pts_box(4)(2);
	    
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
	  
	    if (periodic_z)
	      {
		kmin = iz-3;
		kmax = iz+3;
	      }
	    else
	      {
		kmin = max(0,iz-3);
		kmax = min(iz+3,nbz-1);
	      }
	    
	    // DISP(imin); DISP(imax); DISP(jmin); DISP(jmax);
	    for (int i1 = imin; i1 <= imax; i1++)
	      for (int i2 = jmin; i2 <= jmax; i2++)
		for (int i3 = kmin; i3 <= kmax; i3++)
		  {
		    int i1b = i1%nbx;
		    int i2b = i2%nby;
		    int i3b = i3%nbz;
		    if (i1b < 0)
		      i1b += nbx;
		
		    if (i2b < 0)
		      i2b += nby;
		    
		    if (i3b < 0)
		      i3b += nbz;
		    
		    int nbox = i3b*nbx*nby + i2b*nbx + i1b;
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
  void VlasovMaxwell_3D<TypeElement,TypeEquation>::GeneratePeriodicDirections()
  {
    Mesh<Dimension3>& mesh = this->var_harmonic.mesh;
    MeshNumbering<Dimension3>& mesh_num = this->var_harmonic.mesh_num;
    
    // periodicity vectors
    R3 vec, point, vec2, vec3;
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
    else if (mesh.GetNbPeriodicReferences() == 3)
      {
	this->periodic_vector.Reallocate(27);
	this->periodic_vector(0).Zero();
	vec = mesh_num.GetTranslationPeriodicBoundary(0);
	vec2 = mesh_num.GetTranslationPeriodicBoundary(1);
	vec3 = mesh_num.GetTranslationPeriodicBoundary(2);
	
	this->periodic_vector(1) = vec;
	Mlt(-1.0, vec); this->periodic_vector(2) = vec;
	this->periodic_vector(3) = vec2;
	Mlt(-1.0, vec2); this->periodic_vector(4) = vec2;
	this->periodic_vector(5) = vec3;
	Mlt(-1.0, vec3); this->periodic_vector(6) = vec2;
	
	Add(vec, vec2, this->periodic_vector(7));
	Subtract(vec, vec2, this->periodic_vector(8));
	Mlt(-1.0, vec);
	Add(vec, vec2, this->periodic_vector(9));
	Subtract(vec, vec2, this->periodic_vector(10));
	
	Add(vec, vec3, this->periodic_vector(11));
	Subtract(vec, vec3, this->periodic_vector(12));
	Mlt(-1.0, vec);
	Add(vec, vec3, this->periodic_vector(13));
	Subtract(vec, vec3, this->periodic_vector(14));
	
	Add(vec3, vec2, this->periodic_vector(15));
	Subtract(vec3, vec2, this->periodic_vector(16));
	Mlt(-1.0, vec3);
	Add(vec3, vec2, this->periodic_vector(17));
	Subtract(vec3, vec2, this->periodic_vector(18));
	
	Add(vec, vec2, point); Add(point, vec3, this->periodic_vector(19));
	Mlt(-1.0, vec3); Add(point, vec3, this->periodic_vector(20));
	Mlt(-1.0, vec2);
	Add(vec, vec2, point); Add(point, vec3, this->periodic_vector(21));
	Mlt(-1.0, vec3); Add(point, vec3, this->periodic_vector(22));
	Mlt(-1.0, vec);
	Add(vec, vec2, point); Add(point, vec3, this->periodic_vector(23));
	Mlt(-1.0, vec3); Add(point, vec3, this->periodic_vector(24));
	Mlt(-1.0, vec2);
	Add(vec, vec2, point); Add(point, vec3, this->periodic_vector(25));
	Mlt(-1.0, vec3); Add(point, vec3, this->periodic_vector(26));

      }
  }
  
  
  //! returns offset in the vector of unknowns to access to the parameters of the first particle
  template<class TypeElement, class TypeEquation>
  int VlasovMaxwell_3D<TypeElement,TypeEquation>::GetOffset_PosParticle() const
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
  void VlasovMaxwell_3D<TypeElement,TypeEquation>::
  EvaluateE_H_Quadrature(const VectReal_wp& En, const VectReal_wp& Hn,
			 VectR3& En_quad, VectR3& Hn_quad,
			 bool normal_quadrature, GhostIf<false>& dg_form)
  {
    const Mesh<Dimension3>& mesh = this->var_harmonic.mesh;
    const MeshNumbering<Dimension3>& mesh_num = this->var_harmonic.mesh_num;
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
    VectReal_wp Eloc(N);
    VectR3 Enode(Nh), Equad(Nquad), Hnode(Nh);
    R3 vec_u, vec_v; Matrix3_3 mat_dfj; VectR3 Hquad(Nquad); 
    
    // loop over elements
    for (int i = 0; i < nb_elt; i++)
      {
	// Eloc -> values of E on local dofs
	abort();
        /* for (int j = 0; j < N; j++)
	  {
	    int num_dof = mesh_num.Element(i).GetNumberDof(j);
	    if (var.IsNegativeSign(i, j))
	      Eloc(j) = -En(num_dof);
	    else
	      Eloc(j) = En(num_dof);
              }*/
	
	// we get values of H on the element i
	for (int j = 0; j < Nh; j++)
	  {
	    Hnode(j)(0) = Hn(3*(i*Nh + j));
	    Hnode(j)(1) = Hn(3*(i*Nh + j)+1);
	    Hnode(j)(2) = Hn(3*(i*Nh + j)+2);
	  }
	
	// computation of E on nodal points
	Fb.ComputeNodalValuesRef(Eloc, Enode);
	
	// projection on quadrature points
	if (!normal_quadrature)
	  {
	    this->proj_nodal_to_quadrature.Project(Enode, Equad);
	    this->proj_nodal_to_quadrature.Project(Hnode, Hquad);
	  }
	else
	  {
	    Equad = Enode;
	    Hquad = Hnode;
	  }
	
	// multiplying by DF_i^{*-1}
	if (normal_quadrature)
	  {
	    if (this->var_harmonic.Glob_DFjm1(i).GetM() == 1)
	      {
		mat_dfj = this->var_harmonic.Glob_DFjm1(i)(0);
		jacob = sqrt(Det(mat_dfj)); jacob = 1.0/jacob;
		for (int j = 0; j < Nquad; j++)
		  {
		    int num = i*Nquad+j;
		    MltTrans(mat_dfj, Equad(j), En_quad(num));
		    Mlt(jacob, En_quad(num));
		    Hn_quad(num) = Hquad(j);
		  }
	      }
	    else
	      for (int j = 0; j < Nquad; j++)
		{
		  int num = i*Nquad+j;
		  mat_dfj = this->var_harmonic.Glob_DFjm1(i)(j);
		  jacob = sqrt(Det(mat_dfj)); jacob = 1.0/jacob;
		  MltTrans(mat_dfj, Equad(j), En_quad(num));
		  Mlt(jacob, En_quad(num));
		  Hn_quad(num) = Hquad(j);
		}
	  }
	else
	  for (int j = 0; j < Nquad; j++)
	    {
	      int num = i*Nquad+j;
	      mat_dfj = this->Glob_DFjm1(i)(j);
	      jacob = sqrt(Det(mat_dfj)); jacob = 1.0/jacob;
	      MltTrans(mat_dfj, Equad(j), En_quad(num));
	      Mlt(jacob, En_quad(num));
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
  void VlasovMaxwell_3D<TypeElement,TypeEquation>::
  EvaluateE_H_Quadrature(const VectReal_wp& En, const VectReal_wp& Hn,
			 VectR3& En_quad, VectR3& Hn_quad,
			 bool normal_quadrature, GhostIf<true>& dg_form)
  {
    const Mesh<Dimension3>& mesh = this->var_harmonic.mesh;
    const MeshNumbering<Dimension3>& mesh_num = this->var_harmonic.mesh_num;
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
    
    EllipticProblem<TypeElement,TypeEquationStationary>& var = this->var_harmonic;
    Real_wp jacob; int Nall = this->PointsQuadrature.GetM();
    En_quad.Reallocate(Nall); Hn_quad.Reallocate(Nall);
    VectR3 Hloc(Nh);
    VectR3 Enode(Nh), Equad(Nquad), Hquad(Nquad); 
    
    // loop over elements
    for (int i = 0; i < nb_elt; i++)
      {
	// Eloc -> values of E on local dofs
	for (int j = 0; j < Nh; j++)
	  {
	    int num = i*Nh + j;
	    Enode(j)(0) = En(num);
	    Enode(j)(1) = En(num + Nvol);
	    Enode(j)(2) = En(num + 2*Nvol);
	    Hloc(j)(0) = Hn(num);
	    Hloc(j)(1) = Hn(num + Nvol);
	    Hloc(j)(2) = Hn(num + 2*Nvol);
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
  void VlasovMaxwell_3D<TypeElement,TypeEquation>::
  EvaluateE_H_QuadratureBoundary(const VectReal_wp& En, const VectReal_wp& Hn,
				 VectR3& Equad, VectR3& Hquad, GhostIf<true>& dg_form)
  {
    const Mesh<Dimension3>& mesh = this->var_harmonic.mesh;
    const MeshNumbering<Dimension3>& mesh_num = this->var_harmonic.mesh_num;
    int num_edge, num_elem, nb_pts_quad_edge, nb_dof_edge, num_loc;
    int r = mesh_num.GetOrder();
    const TypeElement& Fb = this->var_harmonic.GetFaceBasis(r);
    Equad.Reallocate((r+1)*(r+1)*mesh.GetNbBoundaryRef());
    Hquad.Reallocate((r+1)*(r+1)*mesh.GetNbBoundaryRef());
    int Nh = Fb.GetNbPointsQuadratureInside();
    int Nvol = mesh_num.GetNbDof();
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      {
	this->var_harmonic.GetBoundaryValues(mesh, i, num_edge, num_elem,
					     nb_pts_quad_edge, nb_dof_edge, num_loc);
	
	for (int k = 0; k < nb_pts_quad_edge; k++)
	  {
	    int num = Fb.GetQuadNumber(num_loc, k) + num_elem*Nh;
	    int num_quad = i*nb_dof_edge + k;
	    Equad(num_quad)(0) = En(num);
	    Equad(num_quad)(1) = En(num+Nvol);
	    Equad(num_quad)(2) = En(num+2*Nvol);
	    Hquad(num_quad)(0) = Hn(num);
	    Hquad(num_quad)(1) = Hn(num+Nvol);
	    Hquad(num_quad)(2) = Hn(num+2*Nvol);
	  }
      }
  }
  
  
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_3D<TypeElement,TypeEquation>::
  EvaluateE_H_QuadratureBoundary(const VectReal_wp& En, const VectReal_wp& Hn,
				 VectR3& Equad, VectR3& Hquad, GhostIf<false>& dg_form)
  {
    const Mesh<Dimension3>& mesh = this->var_harmonic.mesh;
    /*int num_edge, num_elem, nb_pts_quad_edge, nb_dof_edge, num_loc;
    int r = mesh_num.GetOrder();
    const TypeElement& Fb = this->var_harmonic.GetFaceBasis(r);
    VectReal_wp Uloc(Fb.GetNbDof());
    VectR3 Unode(Fb.GetNbPointsNodalElt());
    Equad.Reallocate((r+1)*(r+1)*mesh.GetNbBoundaryRef());
    Hquad.Reallocate((r+1)*(r+1)*mesh.GetNbBoundaryRef());
    int Nh = Fb.GetNbPointsQuadratureInside();
    Real_wp jacob; R3 vec_u; */
    for (int i = 0; i < mesh.GetNbBoundaryRef(); i++)
      {
        abort();
        /*
	this->var_harmonic.GetBoundaryValues(mesh, i, num_edge, num_elem,
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
	
	for (int k = 0; k < nb_pts_quad_edge; k++)
	  {
	    int num = Fb.GetQuadNumber(num_loc, k);
	    int num_quad = i*nb_pts_quad_edge + k;
	    int num2 = 0;
	    if (this->var_harmonic.Glob_DFjm1(num_elem).GetM() > 1)
	      num2 = num;
	    
	    jacob = sqrt(Det(this->var_harmonic.Glob_DFjm1(num_elem)(num2)));
	    MltTrans(this->var_harmonic.Glob_DFjm1(num_elem)(num2), Unode(num), vec_u);
	    Mlt(1.0/jacob, vec_u);
	    Equad(num_quad) = vec_u;
	    Hquad(num_quad)(0) = Hn(3*(num + num_elem*Nh));
	    Hquad(num_quad)(1) = Hn(3*(num + num_elem*Nh)+1);
	    Hquad(num_quad)(2) = Hn(3*(num + num_elem*Nh)+2);
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
  void VlasovMaxwell_3D<TypeElement,TypeEquation>::
  EvaluateE_H_Interpolate(const VectR3& En_quad, const VectR3& Hn_quad,
			  VectR3& En_interpol, VectR3& Hn_interpol)
  {
    Vector<VectReal_wp>& Jacobian = this->Glob_jacobian;
    int nb_pts_interpol = this->nb_points_interpolation_current;
    int Nall = this->PointsInterpolation.GetM();
    En_interpol.Reallocate(Nall); Hn_interpol.Reallocate(Nall);
    int N = this->nb_points_quadrature_current;
    VectR3 Einterpol(nb_pts_interpol), Eloc(N);
    VectR3 Hinterpol(nb_pts_interpol), Hloc(N);
    
    for (int i = 0; i < this->var_harmonic.mesh.GetNbElt(); i++)
      {
	for (int j = 0; j < N; j++)
	  {
	    int nquad = i*N + j;
	    Eloc(j) = En_quad(nquad);
	    Mlt(Jacobian(i)(j), Eloc(j));
	    Hloc(j) = Hn_quad(nquad);
	    Mlt(Jacobian(i)(j), Hloc(j));
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
  void VlasovMaxwell_3D<TypeElement,TypeEquation>::
  AddIntegralCurrentJ(const Real_wp& alpha, const VectR3& Jn_quad,
                      VectReal_wp& Prod_En, GhostIf<false>& dg_form)
  {
    const Mesh<Dimension3>& mesh = this->var_harmonic.mesh;
    const MeshNumbering<Dimension3>& mesh_num = this->var_harmonic.mesh_num;
    int nb_elt = mesh.GetNbElt();
    const TypeElement& Fb = this->var_harmonic.GetFaceBasis(mesh_num.GetOrder());
    int N = Fb.GetNbDof();
    int Nquad = this->nb_points_quadrature_current;
    int Nh = Fb.GetNbPointsQuadratureInside();
    
    VectR3 Jquad(Nquad), Jloc(Nh); VectReal_wp Edof(N);
    R3 vec_u, vec_v;
    for (int i = 0; i < nb_elt; i++)
      {
	  for (int j = 0; j < Nquad; j++)
	    {
	      int num = i*Nquad + j;
	      vec_u = Jn_quad(num);
	      Mlt(this->Glob_DFjm1(i)(j), vec_u, vec_v);
	      Jquad(j) = vec_v;
	      Mlt(this->WeightsQuadrature(j), Jquad(j));
	    }
	  
	  this->proj_nodal_to_quadrature.TransposeProject(Jquad, Jloc);
	
	  // integration 
	  //Fb.ComputeIntegralRef(Jloc, Edof);
	
	  // summation
          abort();
	  /* for (int j = 0; j < N; j++)
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
  void VlasovMaxwell_3D<TypeElement,TypeEquation>::
  AddIntegralCurrentJ(const Real_wp& alpha, const VectR3& Jn_quad,
                      VectReal_wp& Prod_En, GhostIf<true>& dg_form)
  {
    const Mesh<Dimension3>& mesh = this->var_harmonic.mesh;
    const MeshNumbering<Dimension3>& mesh_num = this->var_harmonic.mesh_num;
    int nb_elt = mesh.GetNbElt();
    const TypeElement& Fb = this->var_harmonic.GetFaceBasis(mesh_num.GetOrder());
    int Nquad = this->nb_points_quadrature_current;
    int Nh = Fb.GetNbPointsQuadratureInside();
    int Nvol = mesh_num.GetNbDof();
    
    VectR3 Jquad(Nquad), Jloc(Nh);
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
	      int num_dof = this->var_harmonic.mesh_num.GetNumberDof(i,j);
	      Prod_En(num_dof) += Jloc(j)(0);
	      Prod_En(num_dof+Nvol) += Jloc(j)(1);
	      Prod_En(num_dof+2*Nvol) += Jloc(j)(2);
	    }
      }
  }
  
  
  //! Interpolation of J on quadrature points
  /*!
    \param[in] Jsrc evaluation of J on interpolation points
    \param[in,out] Jquad evaluation of J on quadrature points
   */
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_3D<TypeElement,TypeEquation>::
  EvaluateCurrentQuadrature(const VectR3& Jn_interpol, VectR3& Jn_quad)
  {
    const Mesh<Dimension3>& mesh = this->var_harmonic.mesh;
    int nb_elt = mesh.GetNbElt();
    int Nquad = this->nb_points_quadrature_current;
    int nb_points_interpol = this->nb_points_interpolation_current;
    
    VectR3 Jquad(Nquad), Jloc(nb_points_interpol);
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
  
  
  //! Computation of \f$ \int rho \varphi \f$ (for Boris correction)
  /*!
    \param[in,out] b_rhs vector to which \f$ \int rho \varphi \f$ is added
    \param[in] evalRho evaluation of rho on quadrature/interpolation points
   */
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_3D<TypeElement,TypeEquation>
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
  void VlasovMaxwell_3D<TypeElement,TypeEquation>::
  EvaluateSourceBoris(VectReal_wp& b_rhs, VectReal_wp& evalRho, VectR3& evalEn,
		      VectR3& evalEnBoundary)
  {
    const Mesh<Dimension3>& mesh = this->var_harmonic.mesh;
    
    // part \int_\Omega \rho \phi
    EvaluateSourceBoris(b_rhs, evalRho);
    
    // part \int_\Omega E \nabla \phi
    int nb_pts_quad = this->nb_points_quadrature_current;
    int N = this->var_laplace.mesh_num.GetNbLocalDof(0);
    Real_wp coef; VectR3 feval(nb_pts_quad), Edof(N); R3 vec_u, vec_v;
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
	//this->var_laplace.GetFaceBasis(r).ComputeIntegralGradientRef(Edof, contrib);
	abort();
        
	for (int j = 0; j < N; j++)
	  b_rhs(this->var_laplace.mesh_num.Element(i).GetNumberDof(j)) += contrib(j);
      }
    
    // part -\int_\Gamma E.n \varphi
    //int num_edge, num_elem, nb_pts_quad_edge, nb_dof_edge, num_loc;
    for (int num_edge_ref = 0; num_edge_ref < mesh.GetNbBoundaryRef(); num_edge_ref++)
      {
        abort();
        /*
	this->var_laplace.GetBoundaryValues(mesh, num_edge_ref, num_edge,
					    num_elem, nb_pts_quad_edge, nb_dof_edge, num_loc);
	
	int ref = mesh.Element(num_elem).GetReference();
	Real_wp epsilon = this->var_harmonic.ref_epsilon(ref)(0,0);
	R3 normale;
	
	if (mesh.Boundary(num_edge).GetNbElements() < 2)
	  {
	    for (int k = 0; k < nb_pts_quad_edge; k++)
	      {
		int num = num_edge_ref*nb_pts_quad_edge + k;
		vec_u = evalEnBoundary(num);
		normale = this->Glob_normale(num_edge_ref)(k);
		
		coef = -DotProd(normale, vec_u);
		coef *= epsilon*this->WeightsBoundary(k);
		
		int num_dof = this->var_laplace.mesh_num.Element(num_elem).
                GetNumberDof(this->var_laplace.GetLocalNumber(num_elem, num_loc, k));
		
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
  void VlasovMaxwell_3D<TypeElement,TypeEquation>::
  AddSurfacicIntegral_PhiNodal(const Real_wp& alpha, const VectReal_wp& phi,
                               VectReal_wp& Prod_En, GhostIf<false>& dg_form)
  {
    //int nb_pts_quad_edge, nb_dof_edge, num_loc, num_edge, num_elem, num_dof;
    const Mesh<Dimension3>& mesh = this->var_laplace.mesh;
    //nb_dof_edge = (mesh.GetOrder()+1)*(mesh.GetOrder()+1);
    //int nb_dof_loc = this->var_laplace.mesh.GetNbLocalDof(0);
    //VectR3 feval(nb_dof_edge); VectReal_wp res(nb_dof_loc);
    //const TypeElement& Fb = this->var_harmonic.GetFaceBasis(mesh.GetOrder());
    //R3 vec_u, vec_v; Real_wp val;
    for (int num_edge_ref = 0; num_edge_ref < mesh.GetNbBoundaryRef(); num_edge_ref++)
      {
        abort();
        /*
	this->var_laplace.GetBoundaryValues(mesh, num_edge_ref, num_edge,
					    num_elem, nb_pts_quad_edge, nb_dof_edge, num_loc);
	
	if (mesh.Boundary(num_edge).GetNbElements() < 2)
	  {
	    feval.Fill(R3(0,0));
	    
	    for (int k = 0; k < nb_pts_quad_edge; k++)
	      {
		int num_dof_loc = this->var_laplace.GetLocalNumber(num_elem, num_loc, k);
		num_dof = mesh_num.Element(num_elem).GetNumberDof(num_dof_loc);
		val = phi(num_dof);
		
		int num = num_edge*nb_pts_quad_edge + k;
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
  void VlasovMaxwell_3D<TypeElement,TypeEquation>::
  AddSurfacicIntegral_PhiNodal(const Real_wp& alpha, const VectReal_wp& phi,
                               VectReal_wp& Prod_En, GhostIf<true>& dg_form)
  {
    int nb_pts_quad_edge, nb_dof_edge, num_loc, num_edge, num_elem, num_dof;
    const Mesh<Dimension3>& mesh = this->var_laplace.mesh;
    const MeshNumbering<Dimension3>& mesh_num = this->var_laplace.mesh_num;
    nb_dof_edge = (mesh_num.GetOrder()+1)*(mesh_num.GetOrder()+1);
    int nb_dof_loc = this->var_laplace.GetNbLocalDof(0);
    VectReal_wp fx(nb_dof_edge), fy(nb_dof_edge), fz(nb_dof_edge);
    VectReal_wp res_x(nb_dof_loc), res_y(nb_dof_loc), res_z(nb_dof_edge);
    const TypeElement& Fb = this->var_harmonic.GetFaceBasis(mesh_num.GetOrder());
    R3 vec_u;
    int Nvol = mesh_num.GetNbDof();
    for (int num_edge_ref = 0; num_edge_ref < mesh.GetNbBoundaryRef(); num_edge_ref++)
      {
	this->var_laplace.GetBoundaryValues(mesh, num_edge_ref, num_edge,
					    num_elem, nb_pts_quad_edge, nb_dof_edge, num_loc);
	
	if (mesh.Boundary(num_edge).GetNbElements() < 2)
	  {
	    fx.Fill(0); fy.Fill(0);
	    
	    for (int k = 0; k < nb_pts_quad_edge; k++)
	      {
		int num_dof_loc = this->var_laplace.GetLocalNumber(num_elem, num_loc, k);
		num_dof = mesh_num.Element(num_elem).GetNumberDof(num_dof_loc);
		Real_wp val = phi(num_dof);
		
		int num = num_edge*nb_pts_quad_edge + k;
		vec_u = this->Glob_normale(num_edge)(k);
		Mlt(alpha*WeightsBoundary(k)*val, vec_u);
		
		fx(k) = vec_u(0);
		fy(k) = vec_u(1);
		fz(k) = vec_u(2);
	      }
	    
	    Fb.ComputeIntegralSurface_Ref(fx, res_x, num_loc);
	    Fb.ComputeIntegralSurface_Ref(fy, res_y, num_loc);
	    Fb.ComputeIntegralSurface_Ref(fz, res_z, num_loc);
	    
	    for (int j = 0; j < nb_dof_loc; j++)
	      {
		num_dof = num_elem*nb_dof_loc + j;
		Prod_En(num_dof) += res_x(j); 
		Prod_En(Nvol+num_dof) += res_y(j); 
		Prod_En(2*Nvol+num_dof) += res_z(j); 
	      }
	  }
      }
  }
  
  
  //! adding \f$ \int \nabla \phi \cdot \varphi \f$
  /*!
    \param[in] alpha coefficient
    \param[in] phi vector phi which lies on the nodal space H^1
    \param[in,out] Prod_En vector to which integral is added
   */
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_3D<TypeElement,TypeEquation>::
  AddVolumetricIntegral_PhiNodal(const Real_wp& alpha, const VectReal_wp& phi,
                                 VectReal_wp& Prod_En, GhostIf<false>& dg_form)
  {
    // number of elements and quadrature points on each element
    const Mesh<Dimension3>& mesh = this->var_harmonic.mesh;
    const MeshNumbering<Dimension3>& mesh_num = this->var_harmonic.mesh_num;
    const TypeElement& Fb = this->var_harmonic.GetFaceBasis(mesh_num.GetOrder());
    const HexahedronLobatto& Fb_laplace = this->var_laplace.GetFaceBasis(mesh_num.GetOrder());
    int nb_elt = mesh.GetNbElt();
    
    // variables for interpolation
    int r = mesh_num.GetOrder(), N = r+1;
    int nb_dof = N*N*N;
    VectReal_wp Uloc(nb_dof);
    Vector<R3> Enode(nb_dof), grad_phi(nb_dof);
    VectReal_wp contrib(Fb.GetNbDof());
    R3 vec_u, vec_v, vec_w; Matrix3_3 dfjm1;
    
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
	    Real_wp jacob = sqrt(Det(dfjm1));
	    vec_u = grad_phi(j);
	    MltTrans(dfjm1, vec_u, vec_w);
	    Mlt(dfjm1, vec_w, vec_v);
	    jacob = alpha*this->var_harmonic.Glob_jacobian(i)(j)/(jacob*jacob);
	    Enode(j) = vec_v; Mlt(jacob, Enode(j));
	  }
	
	//Fb.ComputeIntegralRef(Enode, contrib);
        
        abort();
        /*
	for (int j = 0; j < Fb.GetNbDof(); j++)
	  {
	    int num_dof = this->var_harmonic.mesh_num.Element(i).GetNumberDof(j);
	    if (this->var_harmonic.IsNegativeSign(i,j))
	      Prod_En(num_dof) -= contrib(j);
	    else
	      Prod_En(num_dof) += contrib(j);
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
  void VlasovMaxwell_3D<TypeElement,TypeEquation>::
  AddVolumetricIntegral_PhiNodal(const Real_wp& alpha, const VectReal_wp& phi,
                                 VectReal_wp& Prod_En, GhostIf<true>& dg_form)
  {
    // number of elements and quadrature points on each element
    const Mesh<Dimension3>& mesh = this->var_harmonic.mesh;
    const MeshNumbering<Dimension3>& mesh_num = this->var_harmonic.mesh_num;
    const TypeElement& Fb = this->var_harmonic.GetFaceBasis(mesh_num.GetOrder());
    const HexahedronLobatto& Fb_laplace = this->var_laplace.GetFaceBasis(mesh_num.GetOrder());
    int nb_elt = mesh.GetNbElt();
    
    // variables for interpolation
    int r = mesh_num.GetOrder(), N = r+1;
    int nb_dof = N*N*N;
    VectReal_wp Uloc(nb_dof);
    Vector<R3> Enode(nb_dof), grad_phi(nb_dof);
    R3 vec_u, vec_v; Matrix3_3 dfjm1;
    int Nvol = mesh_num.GetNbDof();
    
    for (int i = 0; i < nb_elt; i++)
      {
	// we get values of phi on the element i
	for (int j = 0; j < nb_dof; j++)
	  Uloc(j) = phi(this->var_laplace.mesh_num.Element(i).GetNumberDof(j));

	// derivation	
	Fb_laplace.ComputeNodalGradientRef(Uloc, grad_phi);
	
	// multiplying by DF_i^{*-1}
	for (int j = 0; j < nb_dof; j++)
	  {
	    dfjm1 = this->var_harmonic.Glob_DFjm1(i)(j);
	    Real_wp jacob = sqrt(Det(dfjm1));
	    vec_u = grad_phi(j);
	    MltTrans(dfjm1, vec_u, vec_v);
	    jacob = alpha*Fb.Weights_ND(j);
	    Mlt(jacob, vec_v); Enode(j) = vec_v;
	  }
	
	for (int j = 0; j < Fb.GetNbDof(); j++)
	  {
	    int num_dof = this->var_harmonic.mesh_num.Element(i).GetNumberDof(j);
	    Prod_En(num_dof) += Enode(j)(0);
	    Prod_En(Nvol+num_dof) += Enode(j)(1);
	    Prod_En(2*Nvol+num_dof) += Enode(j)(2);
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
  void VlasovMaxwell_3D<TypeElement,TypeEquation>::
  AddVolumetricProjection_PhiNodal(const Real_wp& alpha, const VectReal_wp& phi,
				VectReal_wp& Prod_En, GhostIf<false>& dg_form)
  {
    // number of elements and quadrature points on each element
    Mesh<Dimension3>& mesh = this->var_harmonic.mesh;
    MeshNumbering<Dimension3>& mesh_num = this->var_harmonic.mesh_num;
    const TypeElement& Fb = this->var_harmonic.GetFaceBasis(mesh_num.GetOrder());
    const HexahedronLobatto& Fb_laplace = this->var_laplace.GetFaceBasis(mesh_num.GetOrder());
    int nb_elt = mesh.GetNbElt();
    
    VectBool DofUsed(this->var_harmonic.mesh_num.GetNbDof()); DofUsed.Fill(false);
    VectReal_wp contrib(Fb.GetNbDof()); contrib.Fill(0);
    VectR3 grad_phi(Fb_laplace.GetNbDof()), Edof(Fb_laplace.GetNbDof());
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
  void VlasovMaxwell_3D<TypeElement,TypeEquation>::
  AddVolumetricProjection_PhiNodal(const Real_wp& alpha, const VectReal_wp& phi,
				VectReal_wp& Prod_En, GhostIf<true>& dg_form)
  {
    // number of elements and quadrature points on each element
    Mesh<Dimension3>& mesh = this->var_harmonic.mesh;
    MeshNumbering<Dimension3>& mesh_num = this->var_harmonic.mesh_num;
    const TypeElement& Fb = this->var_harmonic.GetFaceBasis(mesh_num.GetOrder());
    const HexahedronLobatto& Fb_laplace = this->var_laplace.GetFaceBasis(mesh_num.GetOrder());
    int nb_elt = mesh.GetNbElt();
    
    VectReal_wp Uloc(Fb_laplace.GetNbDof());
    VectR3 grad_phi(Fb_laplace.GetNbDof());
    R3 vec_u; int Nvol = mesh_num.GetNbDof();
    for (int i = 0; i < nb_elt; i++)
      {
	// we get values of phi on the element i
	for (int j = 0; j < Fb_laplace.GetNbDof(); j++)
	  Uloc(j) = phi(this->var_laplace.mesh_num.Element(i).GetNumberDof(j));
	
	// computation of gradient
	Fb_laplace.ComputeNodalGradientRef(Uloc, grad_phi);
		
	for (int j = 0; j < Fb.GetNbDof(); j++)
	  {
	    Real_wp jacob = sqrt(Det(this->var_harmonic.Glob_DFjm1(i)(j)));
	    MltTrans(this->var_harmonic.Glob_DFjm1(i)(j), grad_phi(j), vec_u);
	    int num_dof = mesh_num.Element(i).GetNumberDof(j);
	    Prod_En(num_dof) += alpha*vec_u(0);
	    Prod_En(num_dof+Nvol) += alpha*vec_u(1);
	    Prod_En(num_dof+2*Nvol) += alpha*vec_u(2);
	  }
      }
  }
  
  
  //! projection of rho on nodal dofs
  /*!
    \param[out] phi result of projection
    \param[in] evalRho evaluation of rho on quadrature/interpolate points
  */
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_3D<TypeElement,TypeEquation>
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
  void VlasovMaxwell_3D<TypeElement,TypeEquation>::
  SetUniform_MagneticField(const R3& valH, VectReal_wp& Y, GhostIf<false>& dg_form)
  {
    int N = this->var_harmonic.GetNbPointsQuadratureInside(0);
    // uniform magnetic field
    for (int i = 0; i < this->var_harmonic.mesh.GetNbElt(); i++)
      for (int j = 0; j < N; j++)
	{
	  int offset = this->var_harmonic.mesh_num.GetNbDof() + 3*(i*N + j);
	  Y(offset) = valH(0);
	  Y(offset+1) = valH(1);
	  Y(offset+2) = valH(2);
	}
  }
  
  
  //! setting the magnetic field to an uniform value
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_3D<TypeElement,TypeEquation>::
  SetUniform_MagneticField(const R3& valH, VectReal_wp& Y, GhostIf<true>& dg_form)
  {
    int N = this->var_harmonic.GetNbPointsQuadratureInside(0);
    int Nvol = this->var_harmonic.mesh_num.GetNbDof();
    // uniform magnetic field
    for (int i = 0; i < this->var_harmonic.mesh.GetNbElt(); i++)
      for (int j = 0; j < N; j++)
	{
	  Y(3*Nvol + i*N + j) = valH(0);
	  Y(4*Nvol + i*N + j) = valH(1);
	  Y(5*Nvol + i*N + j) = valH(2);
	}
  }
  
  
  //! not used
  template<class TypeElement,class TypeEquation>
  Real_wp VlasovMaxwell_3D<TypeElement,TypeEquation>
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
  void VlasovMaxwell_3D<TypeElement,TypeEquation>::RunTimeIterations()
  {
    HyperbolicProblem<TypeElement, TypeEquation>& var_leaf = this->GetLeafClass();
    CurrentSource_Vlasov<TypeElement,TypeEquation> Jsrc(var_leaf);
    GhostIf<TypeEquationStationary::FormulationDG> dg_form;
    IVect NumBoxGrid_Point;
    int size_E = this->GetNbScalarUnknowns(dg_form);
    int size_H = this->GetNbVectorialUnknowns(dg_form);
    int offset = this->GetOffset_PosParticle();
    
    int Npart, num_elem, nb_iter = 0;
    R3 pt_loc, pk, vk, pa, va, vb, valE, valH, valHs;
    Real_wp poids_acc, poids_J, t = this->initial_time;
    VectR3 evalEn, evalEnBoundary; VectR3 evalHn, evalHnBoundary;
    Real_wp dt = this->deltat, delta;
    int nb; Real_wp q2m_dt, gamma, gamma2;
    Real_wp one(1), zero(0);
    
    GridInterpolation<Dimension3>& grid = this->loc_particle;
    VectR3 Jn_quad, En_quad, Hn_quad;
    if (this->var_time_scheme.GetTimeSchemeType() == TimeSchemeEnum::LEAP_FROG_SYSTEM)
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
	    
	    Npart = 0; Jsrc.evalJ.Fill(R3()); Jsrc.evalRho.Fill(0);
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
			pk(0) = X(nb+6*j+3); pk(1) = X(nb+6*j+4); pk(2) = X(nb+6*j+5);
			// deducing the velocity
			this->GetVelocity_FromMomentum(pk, vk); 
			
			// updating value of J
			poids_J = this->weight_particle(i)(j)*this->charge_species(i);
			Jsrc.UpdateCurrent(poids_J, vk, grid.GlobalCoord(Npart),
					   num_elem, pt_loc);
		      }
		    Npart++;
		  }
		nb += 6*this->nb_max_particles(i);
	      }
	    
	    // term of hyperbolic correction S_h \phi^{n+1/2}
	    if (this->type_correction_divergence == this->CORRECTION_HYPERBOLIC)
	      {
		Phin.SetData(this->var_laplace.mesh_num.GetNbDof(), &X(size_E+size_H));
		
		//var_leaf.AddVolumetricIntegral_PhiNodal(-one, Phin, Prod_En, dg_form);
		//var_leaf.AddSurfacicIntegral_PhiNodal(one, Phin, Prod_En, dg_form);
                abort();
                
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
			pk(0) = X(nb+6*j+3); pk(1) = X(nb+6*j+4); pk(2) = X(nb+6*j+5);
			// deducing the velocity
			this->GetVelocity_FromMomentum(pk, vk); 
			
			// updating value of J
			poids_J = this->weight_particle(i)(j)*this->charge_species(i);
			Jsrc.UpdateCurrent(poids_J, vk, grid.GlobalCoord(Npart),
					   num_elem, pt_loc);
			
			valE.Zero(); valH.Zero();
			for (int num_point = 0; num_point <
                               this->CloudElem_QuadraturePoint(num_elem).GetM(); num_point++)
			  {
			    Real_wp coefS = Jsrc.CoefficientShape(num_point);
			    if (coefS != Real_wp(0))
			      {
				int num = this->CloudElem_QuadraturePoint(num_elem)(num_point);
				Add(coefS, evalEn(num), valE);
				Add(coefS, evalHn(num), valH);
			      }
			    
			  }
			
			valH *= this->var_harmonic.mu0(0,0);
			valHs(0) = valH(0)*valH(0); valHs(1) = valH(1)*valH(1);
                        valHs(2) = valH(2)*valH(2);
			
			// advancing momentum p
			// (p^{n+3/2} - p^{n+1/2}) = q/m (E^{n+1} + v \times mu_0 H^{n+1})
			
			// p_a = p^{n+1/2} + q/2m dt E^{n+1} (advancing of E of dt/2)
			pa = pk; Add(q2m_dt, valE, pa);
			
			// now solving dp/dt = q/m v \times B (rotation)
			// (p_b - p_a)/dt = q/m (v_a + v_b)/2 \times B
			this->GetVelocity_FromMomentum(pa, va, gamma, gamma2);
			delta = 1.0/(gamma2+valHs(0)+valHs(1)+valHs(2));
			vb(0) = (gamma2 + valHs(0) - valHs(1) - valHs(2))*va(0)
			  + 2.0*(valH(0)*valH(1)+gamma*valH(2))*va(1)
			  + 2.0*(valH(0)*valH(2)-gamma*valH(1))*va(2);
			vb(1) = 2.0*(valH(0)*valH(1)-gamma*valH(2))*va(0)
			  + (gamma2 + valHs(1) - valHs(0) - valHs(2))*va(1)
			  + 2.0*(valH(1)*valH(2)+gamma*valH(0))*va(2);
			vb(2) = 2.0*(valH(0)*valH(2)+gamma*valH(1))*va(0)
			  + 2.0*(valH(1)*valH(2)-gamma*valH(0))*va(1)
			  +(gamma2+valHs(2) - valHs(0) - valHs(1))*va(2);
			
			// p^{n+3/2} = p_b + q/2m dt E^{n+1} (advancing of E of dt/2)
			pk = vb; Mlt(delta*gamma, pk); Add(q2m_dt, valE, pk);
			X(nb+6*j+3) = pk(0); X(nb+6*j+4) = pk(1); X(nb+6*j+5) = pk(2);
			
			// (x^{n+2}-x^{n+1})/dt = v^{n+3/2}
			this->GetVelocity_FromMomentum(pk, vk);
			X(nb+6*j) += this->deltat*vk(0); 
			X(nb+6*j+1) += this->deltat*vk(1);
			X(nb+6*j+2) += this->deltat*vk(2);			
		      }
		    Npart++;
		  }
		nb += 6*this->nb_max_particles(i);
	      }
	    
	    // term of hyperbolic correction (\phi^{n+3/2} - \phi^{n+1/2})/dt
            // = S_h^* E^{n+1} + \rho^{n+1}
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
      TimeMaxwell_3D<TypeElement,TypeEquation>::RunTimeIterations();
    
  }
  
  //! evaluation of operator G, so that time-evolution system reas as dY/dt = G(t,Y)
  /*!
    \param[in] tn time t
    \param[in] nb_deriv n-th derivative of G with respect to t (nb_deriv = n)
    \param[in] X input vector
    \param[out] Y output vector containing the result of G(t,X)
  */
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_3D<TypeElement,TypeEquation>
  ::EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv, const VectReal_wp& X,
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
    GridInterpolation<Dimension3> & grid = this->loc_particle;
    this->LocalizeParticles(X, NumBoxGrid_Point);
        
    // evaluation of E and H on all the quadrature points
    VectR3 evalEn, evalHn;
    if (this->order_interpolation_current > 0)
      {
	VectR3 En_quad, Hn_quad;
	this->EvaluateE_H_Quadrature(En, Hn, En_quad, Hn_quad, false, dg_form);
	this->EvaluateE_H_Interpolate(En_quad, Hn_quad, evalEn, evalHn);
      }
    else
      this->EvaluateE_H_Quadrature(En, Hn, evalEn, evalHn, false, dg_form);

    // now we compute J at the quadrature points
    CurrentSource_Vlasov<TypeElement,TypeEquation> Jsrc(var_leaf);
    int Npart = 0; R3 vk, pk, pt_loc; int nb = offset;
    Real_wp poids_J, poids_acc;
    R3 valE, valH;    VectR3 Jn_quad;
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
		pk(0) = X(nb+6*j+3); pk(1) = X(nb+6*j+4); pk(2) = X(nb+6*j+5);
		// deducing the velocity
		this->GetVelocity_FromMomentum(pk, vk); 
		// DISP(pk); DISP(vk);
		
		// updating value of J
		poids_J = this->weight_particle(i)(j)*this->charge_species(i);
		Jsrc.UpdateCurrent(poids_J, vk, grid.GlobalCoord(Npart), num_elem, pt_loc);
		
		// computation of the mean value of E and H
		valE.Zero(); valH.Zero();
		for (int num_point = 0; num_point <
                       this->CloudElem_QuadraturePoint(num_elem).GetM(); num_point++)
		  {
		    Real_wp coefS = Jsrc.CoefficientShape(num_point);
		    if (coefS != Real_wp(0))
		      {
			int num = this->CloudElem_QuadraturePoint(num_elem)(num_point);
			Add(coefS, evalEn(num), valE);
			Add(coefS, evalHn(num), valH);
		      }
		    // DISP(evalHn(i_elt,jloc));
		  }
		
		// multypling by mu0 to obtain B and not H
		valH *= this->var_harmonic.mu0(0,0);
		// DISP(Npart); DISP(vk); DISP(X(nb+4*j)); DISP(val_H); DISP(val_Ex); DISP(val_Ey);
		
		// val_H = 10;
		// DISP(val_H); DISP(val_Ex); DISP(val_Ey); DISP(vk); DISP(pk);
		// moving the particle
		// speed part dx_k/dt = v_k
		Y(nb+6*j) = vk(0); Y(nb+6*j+1) = vk(1); Y(nb+6*j+2) = vk(2);
		
		// acceleration part dp_k/dt = q ( E + v \times B)
		Y(nb+6*j+3) = poids_acc*(valE(0) + vk(1)*valH(2) - vk(2)*valH(1));
		Y(nb+6*j+4) = poids_acc*(valE(1) - vk(0)*valH(2) + vk(2)*valH(0));
		Y(nb+6*j+5) = poids_acc*(valE(2) + vk(0)*valH(1) - vk(1)*valH(0));
		// DISP( poids_acc*(val_Ex + vk(1)*val_H));
	      }
	    
	    Npart++;
	  }
	nb += 6*this->nb_max_particles(i);
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
	VectR3 evalHnBoundary, evalEnBoundary;
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
  void VlasovMaxwell_3D<TypeElement,TypeEquation>::
  EvaluateDerivativeFunction(const Real_wp& tn, int nb_deriv, const Real_wp& alpha, int level,
                             const VectReal_wp& X, VectReal_wp& Y, bool invert_mass, bool source)
  {
    abort();
  }
  

  //! localization of particles on the grid
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_3D<TypeElement,TypeEquation>
  ::LocalizeParticles(const VectReal_wp& X, IVect& NumBoxGrid_Point)
  {
    int Npart = this->GetNumberParticles();
    
    // we get the coordinates of all the particles to be localized
    GridInterpolation<Dimension3> & grid = this->loc_particle;
    grid.GlobalCoord.Reallocate(Npart);
    int offset = this->GetOffset_PosParticle(); Npart = 0; int nb = offset;
    for (int i = 0; i < this->nb_particles.GetM(); i++)
      {
	for (int j = 0; j < this->nb_particles(i); j++)
	  {
	    grid.GlobalCoord(Npart)(0) = X(nb+6*j);
	    grid.GlobalCoord(Npart)(1) = X(nb+6*j+1);
	    grid.GlobalCoord(Npart)(2) = X(nb+6*j+2);
	    Npart++;
	  }
	nb += 6*this->nb_max_particles(i);
      }
    
    // if periodicity, points can be shifted
    this->MovePoints_Periodicity(grid.GlobalCoord); // DISP(grid.GlobalCoord);
    Vector<IVect> ListPoints_Grid; IVect NbPoints_Grid;
    grid.PreLocalizePoints(NumBoxGrid_Point, ListPoints_Grid, NbPoints_Grid, false);
    
  }

  
  //! Evaluation of rho
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_3D<TypeElement,TypeEquation>
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
    GridInterpolation<Dimension3> & grid = this->loc_particle;
    
    //Real_wp charge;
    int Npart = 0; int nb = 0; R3 vk, pt_loc;
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
	
	nb += 6*this->nb_max_particles(i);
      }
    
    this->ComputeValue_PhiNodal(bsrc, Jsrc.evalRho);
  }
  
  
  // RUNNING TIME STEPS //
  ////////////////////////
  
  ///////////////////////////
  // DIVERGENCE CORRECTION //
  
  //! writing snapshots if necessary, and position/velocity of the first particle
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_3D<TypeElement,TypeEquation>
  ::WriteSnapshot(int nb_iter, const Real_wp& t, const VectReal_wp& Uh)
  {
    abort();
    //TimeMaxwell_3D<TypeElement,TypeEquation>::WriteSnapshot(nb_iter, t, Uh);
    int offset = this->GetOffset_PosParticle();
    if (nb_iter%10 == 0)
      cout<<"Number of particles "<<this->GetNumberParticles()<<endl;
    
    // if there is at least one particle, we write position and velocity of this particle
    if (this->GetNumberParticles() > 0)
      {
	// DISP(offset); DISP(Uh(offset));
	ofstream file_out(this->name_file_position.data(),ios::app);
	file_out<<Uh(offset)<<" "<<Uh(offset+1)<<"  "<<Uh(offset+2)<<"  "
                <<Uh(offset+3)<<"  "<<Uh(offset+4)<<"  "<<Uh(offset+5)<<endl;
	file_out.close();
	
      }
    
    // int test_input; cout<<"we wait"<<endl; cin>>test_input;
  }
  

  template<class TypeElement,class TypeEquation>
  void VlasovMaxwell_3D<TypeElement,TypeEquation>::
  CreateParticlesBreakdownField(const VectReal_wp& En, const VectReal_wp& Hn, int nb_iter,
				VectReal_wp& Y, int offset, VectReal_wp& WeightParticle,
                                VectR3& CoorParticle)
  {
    // to do
  }
  
  
  //! Treatment of a new iterate given by a time scheme
  /*!
    \param[in] nb_iter number of the iteration
    \param[in] tn time
    \param[in,out] Y the new iterate
  */
  template<class TypeElement,class TypeEquation>
  void VlasovMaxwell_3D<TypeElement,TypeEquation>
  ::GiveIterate(int nb_iter, const Real_wp& tn, VectReal_wp& Y)
  {
    HyperbolicProblem<TypeElement, TypeEquation>& var_leaf = this->GetLeafClass();
    int offset = this->GetOffset_PosParticle(); int nb = 0;
    
    // DISP(type_correction_divergence); DISP(CORRECTION_BORIS); DISP(Norm2(Y));
    
    Mesh<Dimension3>& mesh = this->var_harmonic.mesh;
    GhostIf<TypeEquationStationary::FormulationDG> dg_form;
    int size_E = this->GetNbScalarUnknowns(dg_form); 
    int size_H = this->GetNbVectorialUnknowns(dg_form);
    
    VectReal_wp En(size_E), Hn(size_H); R3 pk;
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
    GridInterpolation<Dimension3> & grid = this->loc_particle;
    grid.GlobalCoord.Reallocate(Npart);
    Npart = 0; nb = offset;
    for (int i = 0; i < this->nb_particles.GetM(); i++)
      {
	for (int j = 0; j < this->nb_particles(i); j++)
	  {
	    grid.GlobalCoord(Npart)(0) = Y(nb+6*j);
	    grid.GlobalCoord(Npart)(1) = Y(nb+6*j+1);
	    grid.GlobalCoord(Npart)(2) = Y(nb+6*j+2);
	    Npart++;
	  }
	nb += 6*this->nb_max_particles(i);
      }
    // DISP(grid.GlobalCoord);
    
    // particles which have crossed periodic boundaries are moved inside the domain
    this->MovePoints_Periodicity(grid.GlobalCoord);
    Npart = 0; nb = offset; R3 vk;
    for (int i = 0; i < this->nb_particles.GetM(); i++)
      {
	for (int j = 0; j < this->nb_particles(i); j++)
	  {
	    Y(nb+6*j) = grid.GlobalCoord(Npart)(0);
	    Y(nb+6*j+1) = grid.GlobalCoord(Npart)(1);
	    Y(nb+6*j+2) = grid.GlobalCoord(Npart)(2);
	    Npart++;
	  }
	nb += 6*this->nb_max_particles(i);
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
		for (int k = 0; k < 6; k++)
		  Y(nb+6*new_number+k) = Y(nb+6*j+k);
		
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
		for (int k = 0; k < 6; k++)
		  Y(nb+6*j+k) = 0;
	      }
	    
	    this->nb_particles(i) = new_number;
	  }
	nb += 6*this->nb_max_particles(i);
      }
    
    // creation of particles if the electric field is larger than breakdown field
    grid.GlobalCoord.Clear();
    VectReal_wp WeighNewParticles;
    this->CreateParticlesBreakdownField(En, Hn, nb_iter, Y, offset,
                                        WeighNewParticles, grid.GlobalCoord);
    
    // creation of particles for current beams
    // DISP(this->Liste_beam_current.GetM());
    // DISP(this->Liste_beam_current(0).delta_iter);
        // creation of particles for current beams
    for (int nbeam = 0; nbeam < this->Liste_beam_current.GetM(); nbeam++)
      if (nb_iter%this->Liste_beam_current(nbeam).delta_iter == 0)
	{
	  // area of the beam ?
	  Real_wp dA = this->Liste_beam_current(nbeam).ptA.
            Distance(this->Liste_beam_current(nbeam).ptB);
	  
          dA *= this->Liste_beam_current(nbeam).ptC.Distance(this->Liste_beam_current(nbeam).ptD);
	  
	  // intensity
	  Real_wp Jb = this->Liste_beam_current(nbeam).I0;
	  
	  // charge of particles
	  Real_wp q0 = this->Liste_beam_current(nbeam).charge;
	  Real_wp coef_qm = this->Liste_beam_current(nbeam).ratio_q0_over_m0;
	  
	  // velocity of particles multiplied by relativistic factor
	  R3 v0 = this->Liste_beam_current(nbeam).v0; 
	  // now computing quantity of movement
	  R3 p0 = v0; // DISP(v0); DISP(p0);
	  
	  // number of particles to be created
	  Real_wp real_part = square(this->Liste_beam_current(nbeam).nb_part_per_unit)*dA;
	  // true number
	  int nb_part = toInteger(ceil(real_part)), nb_px, nb_py;
	  this->var_random.RoundToSquare(nb_part, 1.0, nb_px, nb_py);
	  
	  Real_wp poids = dA*Jb*this->deltat*this->Liste_beam_current(nbeam).delta_iter
            /(nb_part*q0);
	  
	  if ((nb_part > 0)&&(poids > 0))
	    {
	      // particles uniformly distributed
	      R3 ptA = this->Liste_beam_current(nbeam).ptA;
	      R3 ptB = this->Liste_beam_current(nbeam).ptB;
	      R3 ptC = this->Liste_beam_current(nbeam).ptC;
	      R3 ptD = this->Liste_beam_current(nbeam).ptD, ptM;
	      VectR3 NewCoor(nb_part);  VectReal_wp NewWeight(nb_part); NewWeight.Fill(poids);
	      VectR3 velocity_particle(nb_part); velocity_particle.Fill(p0);
	      for (int nx = 0; nx < nb_px; nx++)
		for (int ny = 0; ny < nb_py; ny++)
		  {
		    Real_wp xv = (2.0*Real_wp(nx)+1)/(2.0*nb_px);
		    Real_wp yv = (2.0*Real_wp(ny)+1)/(2.0*nb_py);
		    ptM.Zero(); Add(1.0-xv, ptA, ptM); Add(xv, ptB, ptM);
		    Add(1.0-yv, ptC, ptM); Add(yv, ptD, ptM);
                    Add(-1.0, this->Liste_beam_current(nbeam).center, ptM);
		    NewCoor(nb_py*nx+ny) = ptM;
		    Real_wp lambda = ptM.Distance(this->Liste_beam_current(nbeam).center);
		    lambda = 0.5 + 0.5*lambda;
		    Real_wp eval = this->Liste_beam_current(nbeam).GetFunction(lambda);
		    NewWeight(nb_py*nx+ny) *= eval;
		    // DISP(n); DISP(eval*poids);
		  }
	      // DISP(ptA); DISP(ptB); DISP(ptC); DISP(ptD);
	      // DISP(NewCoor); DISP(NewWeight); DISP(this->Liste_beam_current(nbeam).center);
	      
	      // if (nb_iter == 0)
	      // {
	      // this->var_deriv->ComputeWeightsBeam(this->Liste_beam_current(nbeam),
              // NewCoor, poids, NewWeight, Jsrc);
	      // this->Liste_beam_current(nbeam).WeightBeam = NewWeight;
	      // DISP(poids); DISP(NewWeight); exit(0);
	      // }
	      // else
	      // NewWeight = this->Liste_beam_current(nbeam).WeightBeam;
	      
	      // DISP(nb_part); DISP(poids); DISP(NewCoor); DISP(NewWeight); DISP(q0); DISP(m0);
	      this->AddParticles(nb_part, q0, coef_qm, velocity_particle,
				 NewCoor, NewWeight, offset, Y);
	      
	      Mlt(q0, NewWeight);
	      // Append(grid.GlobalCoord, NewCoor);
	      // Append(WeighNewParticles, NewWeight);
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
          R3 pt_loc;
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
    
    VectR3 evalHnBoundary; VectR3 evalEnBoundary;
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
	    
	    VectR3 evalEn, evalHn;
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
	VectR3 evalEn, evalHn;
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
	VectR3 evalHn, evalEn;
	this->EvaluateE_H_Quadrature(En, Hn, evalEn, evalHn, true, dg_form);
	// computation of the electromagnetic energy
	Real_wp energy_elec(0), energy_Ex(0), energy_Ey(0), energy_Ez(0),
	  energy_Hx(0), energy_Hy(0), energy_Hz(0);
	for (int i = 0; i < mesh.GetNbElt(); i++)
	  {
	    int ref = this->var_harmonic.mesh.Element(i).GetReference();
	    Real_wp eps = this->var_harmonic.ref_epsilon(ref)(0,0);
	    Real_wp mu = this->var_harmonic.ref_mu(ref)(0,0);
	    int N = this->var_harmonic.GetNbPointsQuadratureInside(i);
	    for (int j = 0; j < N; j++)
	      {
		energy_Ex += eps*this->var_harmonic.Glob_jacobian(i)(j)*square(evalEn(i*N+j)(0));
		energy_Ey += eps*this->var_harmonic.Glob_jacobian(i)(j)*square(evalEn(i*N+j)(1));
		energy_Ez += eps*this->var_harmonic.Glob_jacobian(i)(j)*square(evalEn(i*N+j)(2));
		energy_Hx += mu*this->var_harmonic.Glob_jacobian(i)(j)*square(evalHn(i*N+j)(0));
		energy_Hy += mu*this->var_harmonic.Glob_jacobian(i)(j)*square(evalHn(i*N+j)(1));
		energy_Hz += mu*this->var_harmonic.Glob_jacobian(i)(j)*square(evalHn(i*N+j)(2));
		
	      }
	  }
	
	energy_elec = energy_Ex + energy_Ey + energy_Ez + energy_Hx + energy_Hy + energy_Hz;
	
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
		pk(0) = Y(nb+6*j+3); pk(1) = Y(nb+6*j+4); pk(2) = Y(nb+6*j+5);
		
		// incrementing cinetic energy
		if (this->relativist_particle)
		  {
		    gamma = sqrt(1.0 + DotProd(pk,pk)*invC02);
		    energy_cinetic += mass*c02*(gamma-1.0)*this->weight_particle(i)(j);
		  }
		else
		  energy_cinetic += 0.5*mass*this->weight_particle(i)(j)*DotProd(pk,pk);
		
	      }
	    nb += 6*this->nb_max_particles(i);
	  }
	
	ofstream file_out(this->file_output_energy.data(), ios::app);
	file_out<<tn<<" "<<energy_Ex<<" "<<energy_Ey<<" "<<energy_Ez<<" ";
	file_out<<energy_Hx<<" "<<energy_Hy<<" "<<energy_Hz<<" "<<energy_elec<<" "
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
	    Matrix<float> positions(max_part, 6); positions.Fill(0);
	    for (int i = 0; i < this->output_particle_param.GetNature(); i++)
	      nb += 6*this->nb_max_particles(i);
	    
	    for (int j = 0; j < nb_part_species; j++)
	      if (Npart < max_part)
		{
		  positions(Npart,0) = Y(nb+6*j);
		  positions(Npart,1) = Y(nb+6*j+1);
		  positions(Npart,2) = Y(nb+6*j+2);
		  positions(Npart,3) = Y(nb+6*j+3);
		  positions(Npart,4) = Y(nb+6*j+4);
		  positions(Npart,5) = Y(nb+6*j+5);
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
	  /* VectR3 Hdof(this->var_laplace.GetNbDof());
  	  int order = mesh_num.GetOrder();
	  int Nh = (order+1)*(order+1)*(order+1);
	  for (int n = 0; n < mesh.GetNbElt(); n++)
	    for (int j = 0; j < Nh; j++)
	      {
		Hdof(this->var_laplace.mesh.GetNumberDof(n, j))(0) = Hn(3*(n*Nh + j));
		Hdof(this->var_laplace.mesh.GetNumberDof(n, j))(1) = Hn(3*(n*Nh + j)+1);
		Hdof(this->var_laplace.mesh.GetNumberDof(n, j))(2) = Hn(3*(n*Nh + j)+2);
	      }
	  
	  this->var_harmonic.output_grid_param(i).
	    WriteSnapshot(nb_iter, tn, Hdof, this->var_laplace, this->var_laplace.var_grid, 2);
	  */
	}
    // DISP(this->var_harmonic.compteur_temps.GetSecondsChrono(this->var_harmonic.chrono_rigid));
    // DISP(this->var_harmonic.compteur_temps.GetSecondsChrono(this->var_harmonic.chrono_flux));
    // DISP(this->var_harmonic.compteur_temps.
    // GetSecondsChrono(this->var_harmonic.chrono_prod_mat_vect));
    // DISP(this->var_harmonic.compteur_temps.GetSecondsChrono(this->var_harmonic.chrono_mass));
    // DISP(nb_particles); 
  }
  
  
  //! adding particles in the simulation
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_3D<TypeElement, TypeEquation>
  ::AddParticles(int nb_part, Real_wp& q0, Real_wp& ratio_q0_m0,
                 VectR3& p0, VectR3& NewCoor, VectReal_wp& poids, int offset, VectReal_wp& Y)
  {
    int nb_p = 0, nb; R3 point;
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
	    nb += 6*this->nb_max_particles(i);
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
	    
	    Y(nb_slot + 6*num) = point(0);
	    Y(nb_slot + 6*num+1) = point(1);
	    Y(nb_slot + 6*num+2) = point(2);
		
	    // same velocity for all particles of the beam
	    Y(nb_slot + 6*num + 3) = p0(nb_p)(0);
	    Y(nb_slot + 6*num + 4) = p0(nb_p)(1);
	    Y(nb_slot + 6*num + 5) = p0(nb_p)(2);
	    
	    // weight
	    this->weight_particle(num_species)(num) = poids(nb_p);
	    
	    nb_p++;
	  }
	
	// less particles to create now ...
	this->nb_particles(num_species) += min(size_available,nb_part);
	nb_part -= min(size_available,nb_part);
      }
  }
  
  
  //! shifting points so that they belong to computational domain
  template<class TypeElement, class TypeEquation>
  void VlasovMaxwell_3D<TypeElement,TypeEquation>::MovePoints_Periodicity(VectR3& Points)
  {
    if (this->periodic_vector.GetM() > 0)
      {
	Real_wp xmin = this->var_harmonic.GetXmin();
	Real_wp xmax = this->var_harmonic.GetXmax();
	Real_wp ymin = this->var_harmonic.GetYmin();
	Real_wp ymax = this->var_harmonic.GetYmax();
	Real_wp zmin = this->var_harmonic.GetZmin();
	Real_wp zmax = this->var_harmonic.GetZmax();
	R3 point;
	for (int i = 0; i < Points.GetM(); i++)
	  if ((Points(i)(0) < xmin)||(Points(i)(0) > xmax)
              ||(Points(i)(1) < ymin)||(Points(i)(1) > ymax)
	      ||(Points(i)(2) < zmin)||(Points(i)(2) > zmax))
	    {
	      // DISP(i); DISP(Points(i));
	      for (int m = 1; m < this->periodic_vector.GetM(); m++)
		{
		  Add(Points(i), this->periodic_vector(m), point); // DISP(m); DISP(point);
		  if ((point(0) >= xmin)&&(point(0) <= xmax)
                      &&(point(1) >= ymin)&&(point(1) <= ymax)
		      &&(point(2) >= zmin)&&(point(2) <= zmax))
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
  void VlasovMaxwell_3D<TypeElement,TypeEquation>
  ::SetInitialVector(const Real_wp& t_begin, VectReal_wp& Y)
  {
    int offset = this->GetOffset_PosParticle();
    // cout<<" on passe la "<<endl;
    // initialization for the electric field
    VarInstationary<TypeElement,TypeEquation>::SetInitialVector(t_begin, Y);
    
    R3 x0, v0, p0; Real_wp coef_qm, charge;
    VectReal_wp Rx, Ry, Rz, Rnu, Rteta, Rphi; int Npart = 0;
    Real_wp xmin = this->var_harmonic.GetXmin(); Real_wp xmax = this->var_harmonic.GetXmax();
    Real_wp ymin = this->var_harmonic.GetYmin(); Real_wp ymax = this->var_harmonic.GetYmax();
    Real_wp zmin = this->var_harmonic.GetZmin(); Real_wp zmax = this->var_harmonic.GetZmax();
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
	    Vlasov_BeamCurrent<Dimension3>::SetInputSpecies(charge, coef_qm, parameters, nb); 
	    this->charge_species.PushBack(charge); this->ratio_qm_species.PushBack(coef_qm);
	    x0(0) = to_num<Real_wp>(parameters(nb++));
            x0(1) = to_num<Real_wp>(parameters(nb++)); x0(2) = to_num<Real_wp>(parameters(nb++));
	    // v0 is the velocity multiplied by the relativistic factor
	    v0(0) = to_num<Real_wp>(parameters(nb++));
            v0(1) = to_num<Real_wp>(parameters(nb++)); v0(2) = to_num<Real_wp>(parameters(nb++));
	    p0 = v0;
	    // weight
	    weight_zero(0) = to_num<Real_wp>(parameters(nb++));
	    this->weight_particle.PushBack(weight_zero);
	    
	    Y(offset) = x0(0); Y(offset+1) = x0(1); Y(offset+2) = x0(2);
	    Y(offset+3) = p0(0); Y(offset+4) = p0(1); Y(offset+5) = p0(2);
	    
	    offset += 6*this->increment_particles;
	  }
	else if (!parameters(0).compare("CIRCLE"))
	  {
	    int n = to_num<int>(parameters(nb++));
	    if (n > this->increment_particles)
	      {
		cout<<"Too many particles to add "<<endl;
		abort();
	      }
	    R3 center;
	    center(0) = to_num<Real_wp>(parameters(nb++));
	    center(1) = to_num<Real_wp>(parameters(nb++));
            center(2) = to_num<Real_wp>(parameters(nb++));
	    Real_wp radius = to_num<Real_wp>(parameters(nb++));
	    this->nb_particles.PushBack(n);
	    this->nb_max_particles.PushBack(this->increment_particles);
	    VectReal_wp weight_zero(this->increment_particles); weight_zero.Fill(0);
	    for (int i = 0; i < n; i++)
	      weight_zero(i) = 1.0;
	    
	    this->weight_particle.PushBack(weight_zero);
	    Vlasov_BeamCurrent<Dimension3>::SetInputSpecies(charge, coef_qm, parameters, nb);
	    // speed is the velocity multiplied by the relativistic factor
	    Real_wp speed = to_num<Real_wp>(parameters(nb++));
	    this->charge_species.PushBack(charge); this->ratio_qm_species.PushBack(coef_qm);
	    
	    for (int i = 0; i < n; i++)
	      {
		Real_wp teta = Real_wp(i)/n*2.0*pi_wp;
		Y(offset+6*i) = center(0)+radius*cos(teta);
		Y(offset+6*i+1) = center(1)+radius*sin(teta);
		Y(offset+6*i+2) = center(2);
		v0(0) = speed*cos(teta);
		v0(1) = speed*sin(teta); v0(2) = 0;
		Y(offset+6*i+3) = v0(0);
		Y(offset+6*i+4) = v0(1); Y(offset+6*i+5) = v0(2);
	      }
	    
	    offset += 6*this->increment_particles;
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
	    Real_wp vth_z = to_num<Real_wp>(parameters(nb++));
	    Real_wp alpha = to_num<Real_wp>(parameters(nb++));
	    Real_wp kx = to_num<Real_wp>(parameters(nb++));
	    Vlasov_BeamCurrent<Dimension3>::SetInputSpecies(charge, coef_qm, parameters, nb);
	    Real_wp rho0 = to_num<Real_wp>(parameters(nb++));
	    xmin = to_num<Real_wp>(parameters(nb++));
	    xmax = to_num<Real_wp>(parameters(nb++));
	    ymin = to_num<Real_wp>(parameters(nb++));
	    ymax = to_num<Real_wp>(parameters(nb++));
	    zmin = to_num<Real_wp>(parameters(nb++));
	    zmax = to_num<Real_wp>(parameters(nb++));
	    
	    int nx(0), ny(0), nz(0);
	    Real_wp alpha_ = (ymax-ymin)/(xmax-xmin);
	    Real_wp alphaz = (zmax-zmin)/(xmax-xmin);
	    VectReal_wp Rx_glob, Ry_glob, Rz_glob, Rnu_glob, Rteta_glob, Rphi_glob;
	    if (type_position == VarRandomGenerator::UNIFORM)
	      {
		this->var_random.RoundToSquare(Ntot, alpha_, alphaz, nx, ny, nz);
		Rx_glob.Reallocate(Ntot); Ry_glob.Reallocate(Ntot); Rz_glob.Reallocate(Ntot);
		
		Real_wp dx = Real_wp(1)/nx, x0 = 0.5*dx;
		Real_wp dy = Real_wp(1)/ny, y0 = 0.5*dy;
		Real_wp dz = Real_wp(1)/nz, z0 = 0.5*dz;
		int ind = 0;
		for (int i = 0; i < nx; i++)
		  for (int j = 0; j < ny; j++)
		    for (int k = 0; k < nz; k++)
		      {
			Rx_glob(ind) = x0 + i*dx;
			Ry_glob(ind) = y0 + j*dy;
			Rz_glob(ind) = z0 + k*dz; ind++;
		      }
	      }
	    else if (type_position == VarRandomGenerator::BIT_REVERSED)
	      {
		this->var_random.GenerateBitReversedNumbers(Ntot, Rx_glob);
		this->var_random.GenerateBitReversedNumbers(Ntot, Ry_glob);
		this->var_random.GenerateBitReversedNumbers(Ntot, Rz_glob);
	      }
		    
	    if (type_velocity == VarRandomGenerator::UNIFORM)
	      {
		Rnu_glob.Reallocate(Ntot); Rteta_glob.Reallocate(Ntot); Rphi_glob.Reallocate(Ntot);
		
		Real_wp dx = Real_wp(1)/nx, x0 = 0.5*dx;
		Real_wp dy = Real_wp(1)/ny, y0 = 0.5*dy;
		Real_wp dz = Real_wp(1)/ny, z0 = 0.5*dz;
		int ind = 0;
		for (int i = 0; i < nx; i++)
		  for (int j = 0; j < ny; j++)
		    for (int k = 0; k < nz; k++)
		      {
			Rnu_glob(ind) = x0 + i*dx;
			Rteta_glob(ind) = y0 + j*dy;
			Rphi_glob(ind) = z0 + k*dz; ind++;
		      }
		
		this->var_random.ApplyRandomPermutation(Rnu_glob);
		this->var_random.ApplyRandomPermutation(Rteta_glob);
		this->var_random.ApplyRandomPermutation(Rphi_glob);
	      }
	    else if (type_velocity == VarRandomGenerator::BIT_REVERSED)
	      {
		this->var_random.GenerateBitReversedNumbers(Ntot, Rnu_glob);
		this->var_random.GenerateBitReversedNumbers(Ntot, Rteta_glob);
		this->var_random.GenerateBitReversedNumbers(Ntot, Rphi_glob);
	      }
	    
	    // multipication of rho0 by the area of the domain
	    rho0 *= (xmax-xmin)*(ymax-ymin)*(zmax-zmin);
	    
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
		    this->var_random.GenerateRandomNumbers(n, Rz);
		  }
		else
		  {
		    Rx.Reallocate(n); Ry.Reallocate(n); Rz.Reallocate(n);
		    for (int i = 0; i < n; i++)
		      {
			Rx(i) = Rx_glob(Npart);
			Ry(i) = Ry_glob(Npart);
			Rz(i) = Rz_glob(Npart);
			Npart++;
		      }
		  }
		
		if (type_velocity == VarRandomGenerator::RANDOM)
		  {
		    this->var_random.GenerateRandomNumbers(n, Rnu);
		    this->var_random.GenerateRandomNumbers(n, Rteta);
		    this->var_random.GenerateRandomNumbers(n, Rphi);
		  }
		else
		  {
		    Rnu.Reallocate(n); Rteta.Reallocate(n); Rphi.Reallocate(n);
		    for (int i = 0; i < n; i++)
		      {
			Rnu(i) = Rnu_glob(Npart_speed);
			Rteta(i) = Rteta_glob(Npart_speed);
			Rphi(i) = Rphi_glob(Npart_speed);
			Npart_speed++;
		      }
		  }
		
		this->nb_particles.PushBack(n);
		this->nb_max_particles.PushBack(this->increment_particles);
		VectReal_wp weight_zero(this->increment_particles); weight_zero.Fill(0);
		for (int i = 0; i < n; i++)
		  weight_zero(i) = rho0*(xmax-xmin)*(ymax-ymin)*(zmax-zmin)/Ntot;
		
		this->weight_particle.PushBack(weight_zero);
		
		this->charge_species.PushBack(charge); this->ratio_qm_species.PushBack(coef_qm);
		
		Real_wp coef = alpha/kx;
		for (int i = 0; i < n; i++)
		  {
		    Real_wp nu = sqrt(-2.0*log(1-Rnu(i)));
		    Real_wp teta = 2.0*pi_wp*Rteta(i);
		    Real_wp phi = pi_wp*Rphi(i);
		    // inverting distribution function with a Newton method
		    if (abs(alpha) > 1e-10)
		      {
			Real_wp val = (xmax-xmin)*Rx(i), test_prec = 2e20, test = 1e20;
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
		      x0(0) = (xmax-xmin)*Rx(i);
		    
		    Y(offset+6*i) = xmin+x0(0);
		    Y(offset+6*i+1) = ymin+(ymax-ymin)*Ry(i);
		    Y(offset+6*i+2) = zmin+(zmax-zmin)*Rz(i);
		    v0(0) = nu*vth_x*cos(teta)*cos(phi);
		    v0(1) = nu*vth_y*sin(teta)*cos(phi);
		    v0(2) = nu*vth_z*sin(phi);
		    Y(offset+6*i+3) = v0(0);
		    Y(offset+6*i+4) = v0(1);
		    Y(offset+6*i+5) = v0(2);
		  }
		
		offset += 6*this->increment_particles;
	      }
	  }
      }
  }
  
  // CREATION OF PARTICLES //
  ///////////////////////////
  
}

#define MONTJOIE_FILE_VLASOV_MAXWELL3D_CXX
#endif
