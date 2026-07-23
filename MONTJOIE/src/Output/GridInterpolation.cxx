#ifndef MONTJOIE_FILE_GRID_INTERPOLATION_CXX

namespace Montjoie
{
  /**************************
   * GridInterpolation_Base *
   **************************/
  
  //! default constructor
  template<class Dimension>
  GridInterpolation_Base<Dimension>::GridInterpolation_Base() 
  {
    dfjm1_store = true;
    nb_sections_cyclic = 0;
    nb_sections_x = 1; nb_sections_y = 1; nb_sections_z = 1;
    nbPoints_grid = 0; 
    xmin = 1e300; ymin = 1e300; zmin = 1e300;
    xmax = -1e300; ymax = -1e300; zmax = -1e300;
    nb_subdiv_grid_x = 0; nb_subdiv_grid_y = 0; nb_subdiv_grid_z = 0;
    step_subdiv_grid_x = 1.0; step_subdiv_grid_y = 1.0; step_subdiv_grid_z = 1.0;
#ifdef SELDON_WITH_MPI
    comm_ = MPI_COMM_SELF;
#endif
  }
  
  
  //! returns memory used by the object in bytes
  template<class Dimension>
  size_t GridInterpolation_Base<Dimension>::GetMemorySize() const
  {
    size_t taille = sizeof(*this) + GlobalCoord.GetMemorySize()
      + CoorInterp.GetMemorySize() + TetaInterp.GetMemorySize() + dfjm1.GetMemorySize()
      + ElementInterp.GetMemorySize()
      + Epsilon_InsideElement.GetMemorySize() + CyclicSection.GetMemorySize()
      + point_treated.GetMemorySize() + Seldon::GetMemorySize(ListeBoxGrid_Element);
    
    return taille; 
  }
  

  template<class Dimension>
  void GridInterpolation_Base<Dimension>::ResizeElements(int Nb)
  {
    ElementInterp.Resize(Nb);
    CoorInterp.Resize(Nb);
    if (TetaInterp.GetM() > 0)
      TetaInterp.Resize(Nb);
    if (CyclicSection.GetM() > 0)
      CyclicSection.Resize(Nb);
    
    GlobalCoord.Resize(Nb);
    if (dfjm1_store)
      {
	Epsilon_InsideElement.Resize(Nb);
	dfjm1.Resize(Nb);
      }
    
    point_treated.Resize(Nb);
  }
  
  
  //! returns the number of sections (when computations are made on a cyclic domain)
  template<class Dimension>
  int GridInterpolation_Base<Dimension>::GetNbCyclicSections() const
  {
    return nb_sections_cyclic;
  }

  
  //! returns the number of cells along x, y, z (for computation on periodized pattern)
  template<class Dimension>
  void GridInterpolation_Base<Dimension>::GetNbCartesianSections(int& nx, int& ny, int& nz) const
  {
    nx = nb_sections_x;
    ny = nb_sections_y;
    nz = nb_sections_z;
  }

  
  //! sets the number of sections (when computations are made on a cyclic domain)
  template<class Dimension>
  void GridInterpolation_Base<Dimension>::SetNbCyclicSections(int n)
  {
    nb_sections_cyclic = n;
  }


  //! sets the number of cells along x, y, z (for computation on periodized pattern)
  template<class Dimension>
  void GridInterpolation_Base<Dimension>::SetNbCartesianSections(int nx, int ny, int nz)
  {
    nb_sections_x = nx;
    nb_sections_y = ny;
    nb_sections_z = nz;
  }

  
  //! adding points to the list of points that need to be localized
  /*!
    \param[in] Points_ND list of points to be added
    \param[out] number identification numbers of these points in the object
    this functions appends list of points to the previous list of points stored in the object
   */
  template<class Dimension>
  void GridInterpolation_Base<Dimension>::
  Append(const VectR_N& Points_ND, const VectReal_wp& teta, const IVect& SectorNumber, IVect& number)
  {
    int N = Points_ND.GetM();
    number.Reallocate(N);
    int Nold = GlobalCoord.GetM();
    number.Fill();
    for (int i = 0; i < number.GetM(); i++)
      number(i) += Nold;

    GlobalCoord.PushBack(Points_ND);
    TetaInterp.PushBack(teta);
    CyclicSection.PushBack(SectorNumber); 
  }


  //! Clears array used to prelocalize points
  template<class Dimension>
  void GridInterpolation_Base<Dimension>::ClearPrelocalizationArrays()
  {
    ListeBoxGrid_Element.Clear();
  }
  

  //! compress the interpolation grid
  template<class Dimension>
  void GridInterpolation_Base<Dimension>
  ::CompressGrid(Vector<GridInterpolationFull<Dimension> >& grid)
  {
    // objects in grid are compressed first
    int N, nb_pts;
    for (int n = 0; n < grid.GetM(); n++)
      {
	// in grid(n), points which are not localized are removed
	N = grid(n).GetNbPoints();
	nb_pts = 0;
	for (int i = 0; i < N; i++)
	  {
	    int p = grid(n).GetPointNumber(i);
	    if (this->GetElementNumber(p) >= 0)
	      nb_pts++;
	  }
	
	if (nb_pts != N)
	  {
	    IVect old_list(grid(n).list_points);
	    IVect OldIndex(grid(n).IndexPoints);
	    grid(n).list_points.Reallocate(nb_pts);
	    grid(n).IndexPoints.Reallocate(nb_pts);
	    nb_pts = 0;
	    for (int i = 0; i < N; i++)
	      {
		int p = old_list(i);
		if (this->GetElementNumber(p) >= 0)
		  {
		    grid(n).list_points(nb_pts) = p;
		    grid(n).IndexPoints(nb_pts) = OldIndex(i);
		    nb_pts++;
		  }
	      }
	    
	    N = nb_pts;
	  }
	
      }

    // then in the current object, points that are not localized are removed
    N = ElementInterp.GetM();
    nb_pts = 0;
    for (int i = 0; i < N; i++)
      if (ElementInterp(i) >= 0)
	nb_pts++;

    VectR_N OldCoordGlob(GlobalCoord);
    IVect OldElement(ElementInterp);
    VectR_N OldCoor(CoorInterp);
    VectReal_wp OldEpsilon(Epsilon_InsideElement);    
    Vector<bool> old_point(point_treated);
    VectReal_wp OldTeta(TetaInterp);
    IVect OldSection(CyclicSection);
    GlobalCoord.Reallocate(nb_pts);
    TetaInterp.Reallocate(nb_pts);
    CoorInterp.Reallocate(nb_pts);
    ElementInterp.Reallocate(nb_pts);
    Epsilon_InsideElement.Reallocate(nb_pts);
    CyclicSection.Reallocate(nb_pts);
    point_treated.Reallocate(nb_pts);
    nb_pts = 0;
    for (int i = 0; i < N; i++)
      if (OldElement(i) >= 0)
	{
	  GlobalCoord(nb_pts) = OldCoordGlob(i);
	  TetaInterp(nb_pts) = OldTeta(i);
	  CoorInterp(nb_pts) = OldCoor(i);
	  ElementInterp(nb_pts) = OldElement(i);
	  Epsilon_InsideElement(nb_pts) = OldEpsilon(i);
	  CyclicSection(nb_pts) = OldSection(i);
	  point_treated(nb_pts) = old_point(i);
	  nb_pts++;
	}
    
    OldCoordGlob.Clear(); OldCoor.Clear(); OldEpsilon.Clear();
    old_point.Clear(); OldTeta.Clear(); OldSection.Clear();
    
    if (this->dfjm1_store)
      {
	Vector<MatrixN_N> old_dfjm1(dfjm1);
	dfjm1.Reallocate(nb_pts);
	nb_pts = 0;
	for (int i = 0; i < N; i++)
	  if (OldElement(i) >= 0)
	    {
	      dfjm1(nb_pts) = old_dfjm1(i);
	      nb_pts++;
	    }
      }
    
    // numbers in grid are updated
    Vector<int> permut(N);
    permut.Fill(-1);
    nb_pts = 0;
    for (int i = 0; i < N; i++)
      if (OldElement(i) >= 0)
	permut(i) = nb_pts++;
    
    OldElement.Clear();
    for (int n = 0; n < grid.GetM(); n++)
      {
	N = grid(n).GetNbPoints();
	for (int i = 0; i < N; i++)
	  grid(n).list_points(i) = permut(grid(n).list_points(i));
	
	grid(n).list_points_proc.Reallocate(1);
	grid(n).list_points_proc(0) = grid(n).IndexPoints;
      }
	  
#ifdef SELDON_WITH_MPI
    const MPI_Comm& comm = comm_;
    int nb_proc; MPI_Comm_size(comm, &nb_proc);
    int rank_proc; MPI_Comm_rank(comm, &rank_proc);
    if (nb_proc == 1)
      return;
    
    // now each point is conserved in one processor and only one
    for (int n = 0; n < grid.GetM(); n++)
      {
	N = grid(n).GetNbPoints();
	
	// distances to the boundary of elements are sent to the 
	// master processor
	VectReal_wp eps(N);
	for (int i = 0; i < N; i++)
	  {
	    int p = grid(n).GetPointNumber(i);
	    eps(i) = this->GetEpsilonInsideElement(p);
	  }
	
	Vector<VectReal_wp> all_eps;
	Vector<int64_t> xtmp; 	    MPI_Status status;
	if (rank_proc != 0)
	  {
	    MpiSsend(comm, &N, xtmp, 1, 0, 92);
	    MpiSsend(comm, eps, xtmp, N, 0, 93);
	    MpiSsend(comm, grid(n).IndexPoints, xtmp, N, 0, 94);
	  }
	else
	  {
	    all_eps.Reallocate(nb_proc);
	    all_eps(0) = eps;
	    grid(n).list_points_proc.Reallocate(nb_proc);
	    grid(n).list_points_proc(0) = grid(n).IndexPoints;
	    for (int i = 1; i < nb_proc; i++)
	      {
		MpiRecv(comm, &nb_pts, xtmp, 1, i, 92, status);
		all_eps(i).Reallocate(nb_pts);
		grid(n).list_points_proc(i).Reallocate(nb_pts);
		MpiRecv(comm, all_eps(i), xtmp, nb_pts, i, 93, status);
		MpiRecv(comm, grid(n).list_points_proc(i), xtmp, nb_pts, i, 94, status);
	      }
	  }
	
	IVect PointsToRemove;
	if (rank_proc == 0)
	  {
	    // the best processor is selected (positive distance)
	    Vector<int> best_proc(grid(n).nb_global_grid_points);
	    Vector<Real_wp> dist_min(grid(n).nb_global_grid_points);
	    best_proc.Fill(-1); dist_min.Fill(-Real_wp(100));
	    for (int i = 0; i < nb_proc; i++)
	      {
		for (int j = 0; j < all_eps(i).GetM(); j++)
		  {
		    int p = grid(n).list_points_proc(i)(j);
		    if (all_eps(i)(j) > dist_min(p))
		      {
			best_proc(p) = i;
			dist_min(p) = all_eps(i)(j);
		      }
		  }
	      }
	    
	    // sending back the list of points to be removed
	    for (int i = nb_proc-1; i >= 0; i--)
	      {
		int nb_pts_to_remove = 0;
		for (int j = 0; j < all_eps(i).GetM(); j++)
		  {
		    int p = grid(n).list_points_proc(i)(j);
		    if (best_proc(p) != i)
		      nb_pts_to_remove++;
		  }
		
		if (i > 0)
		  MpiSsend(comm, &nb_pts_to_remove, xtmp, 1, i, 96);
		
		if (nb_pts_to_remove > 0)
		  {
		    IVect OldList(grid(n).list_points_proc(i));
		    nb_pts = all_eps(i).GetM() - nb_pts_to_remove;
		    grid(n).list_points_proc(i).Reallocate(nb_pts);
		    PointsToRemove.Reallocate(nb_pts_to_remove);
		    nb_pts_to_remove = 0; nb_pts = 0;
		    for (int j = 0; j < all_eps(i).GetM(); j++)
		      {
			int p = OldList(j);
			if (best_proc(p) != i)
			  {
			    PointsToRemove(nb_pts_to_remove) = j;
			    nb_pts_to_remove++;
			  }
			else
			  {
			    grid(n).list_points_proc(i)(nb_pts) = p;
			    nb_pts++;
			  }
		      }
		    
		    if (i > 0)
		      MpiSsend(comm, PointsToRemove, xtmp, nb_pts_to_remove, i, 97);
		  }
		else
		  PointsToRemove.Clear();
	      }
	  }
	else
	  {
	    int nb_pts_to_remove = 0;
	    MpiRecv(comm, &nb_pts_to_remove, xtmp, 1, 0, 96, status);
	    if (nb_pts_to_remove > 0)
	      {
		PointsToRemove.Reallocate(nb_pts_to_remove);
		MpiRecv(comm, PointsToRemove, xtmp, nb_pts_to_remove, 0, 97, status);
	      }
	  }
	
	if (PointsToRemove.GetM() > 0)
	  {
	    IVect old_list(grid(n).list_points);
	    IVect old_index(grid(n).IndexPoints);
	    nb_pts = old_list.GetM() - PointsToRemove.GetM();
	    grid(n).list_points.Reallocate(nb_pts);
	    grid(n).IndexPoints.Reallocate(nb_pts);
	    nb_pts = 0; int k = 0;
	    for (int i = 0; i < old_list.GetM(); i++)
	      {
		if ((k < PointsToRemove.GetM()) && (i == PointsToRemove(k)))
		  k++;
		else
		  {
		    grid(n).list_points(nb_pts) = old_list(i);
		    grid(n).IndexPoints(nb_pts) = old_index(i);
		    nb_pts++;
		  }
	      }
	  }
      }
#endif

  }
  
  
  //! displaying details of class GridInterpolation_Base
  template<class Dimension>
  ostream& operator<<(ostream& out, const GridInterpolation_Base<Dimension>& g)
  {
    out<<"Grid with "<<g.nbPoints_grid<<endl;
    out<<"Grid included in the box ["<<g.xmin<<","<<g.xmax<<"]x["<<
      g.ymin<<","<<g.ymax<<"]x["<<g.zmin<<","<<g.zmax<<"]"<<endl;
    
    out<<"Center "<<g.center<<endl;
    return out;
  }
  
  
  /******************************
   * GridInterpolationFull_Base *
   ******************************/
  
  
  //! default constructor
  template<class Dimension>
  GridInterpolationFull_Base<Dimension>::GridInterpolationFull_Base()
  {
    var_interp = NULL;
    nb_global_grid_points = 0;
    nbPoints_x = 1; nbPoints_y = 1; nbPoints_z = 1;
    grid_interval_to_be_computed = false;
    xmin = 0.0; xmax = 0.0; ymin = 0.0; ymax = 0.0; zmin = 0.0; zmax = 0.0;
    xmin0 = 0.0; xmax0 = 0.0; ymin0 = 0.0; ymax0 = 0.0; zmin0 = 0.0; zmax0 = 0.0;
    z0_adim = 1.0;
    type_output_file = 0;    
  }
  

  //! returns size of memory used by this object in bytes
  template<class Dimension>
  size_t GridInterpolationFull_Base<Dimension>::GetMemorySize() const
  {
    size_t taille = sizeof(*this) + list_points.GetMemorySize();
    taille += IndexPoints.GetMemorySize() + Seldon::GetMemorySize(list_points_proc);
    taille += Seldon::GetMemorySize(PointsInputFile);
    taille += GlobalPoints2D.GetMemorySize() + GlobalPoints3D.GetMemorySize();
    return taille;
  }
  
  
  //! sets an adimensionalization parameter
  template<class Dimension>
  void GridInterpolationFull_Base<Dimension>::AdimSpace(const Real_wp& z)
  {
    z0_adim = z;
    xmin /= z0_adim;
    xmax /= z0_adim;
    ymin /= z0_adim;
    ymax /= z0_adim;
    zmin /= z0_adim;
    zmax /= z0_adim;
    xmin0 /= z0_adim;
    xmax0 /= z0_adim;
    ymin0 /= z0_adim;
    ymax0 /= z0_adim;
    zmin0 /= z0_adim;
    zmax0 /= z0_adim;
    Mlt(1.0/z0_adim, center);
  }
  
  
  //! displays details of class GridInterpolationFull_Base
  template<class Dimension>
  ostream& operator<<(ostream& out, const GridInterpolationFull_Base<Dimension>& g)
  {
    out<<"Grid with "<<g.nbPoints_x<<" x "<<g.nbPoints_y<<" x "
       <<g.nbPoints_z<<" points"<<endl;
    
    out<<"Grid included in the box ["<<g.xmin<<","<<g.xmax<<"]x["
       <<g.ymin<<","<<g.ymax<<"]x["<<g.zmin<<","<<g.zmax<<"]"<<endl;
    
    out<<"Grid of type "<<g.type_output_file<<" boundaries automatically computed "
       <<g.grid_interval_to_be_computed<<endl;
    
    return out;
  }


  /*********************************
   * GridInterpolation<Dimension1> *
   *********************************/
  
  
#ifdef MONTJOIE_WITH_ONE_DIM
  //! default constructor
  GridInterpolation<Dimension1>::GridInterpolation()
  {
    nb_global_points = 0;
  }
  
  
  //! returns size of memory used by this object in bytes
  size_t GridInterpolation<Dimension1>::GetMemorySize() const
  {
    size_t taille = sizeof(Real_wp)*(Points.GetM()+CoorInterp.GetM());
    taille += sizeof(int)*ElemInterp.GetM();
    taille += Seldon::GetMemorySize(ProcNumberPoints);
    taille += TetaInterp.GetMemorySize() + PhiInterp.GetMemorySize() + sizeof(*this);
    
    return taille;
  }


  //! 1-D output on lines only
  int GridInterpolationFull<Dimension1>::GetType(const string& description)
  {
    if (!description.compare("FileOutputPlaneAxi"))
      return PLANE_AXI;
    else if (!description.compare("FileOutputLineAxi"))
      return LINE_AXI;
    else if (!description.compare("FileOutputPoint"))
      return POINT;
    else if (!description.compare("FileOutputPointAxi"))
      return POINT_AXI;
    else if (!description.compare("FileOutputGrille"))
      return THREE_PLANES_AXI;
    else if (!description.compare("FileOutputGrille3D"))
      return VOLUME_AXI;
    else if (!description.compare("FileOutputCircleAxi"))
      return CIRCLE_AXI;
    else if (!description.compare("FileOutputPointsFileAxi"))
      return POINTS_FILE_AXI;

    return LINE;
  }
    
  
  //! initialization of the interpolation grid with regular points
  void GridInterpolation<Dimension1>::Init(const Real_wp& xmin, const Real_wp& xmax, int N)
  {
    if (N <= 0)
      {
	Points.Clear();
	return;
      }
    
    Points.Reallocate(N);
    if (N == 1)
      {
        Points(0) = xmin;
      }
    else
      {
        Real_wp h = (xmax-xmin)/(N-1);
        for (int i = 0; i < N; i++)
          Points(i) = xmin + h*i;
      }
  }
  
  
  //! initialization of the interpolation grid with non-regular points
  void GridInterpolation<Dimension1>::Init(const VectReal_wp& pos)
  {
    Points = pos;
    // sorting points (if not already done)
    Sort(pos.GetM(), Points);
  }
  

  //! adding points to the list of points that need to be localized
  void GridInterpolation<Dimension1>::
  Append(const VectReal_wp& Points_ND, const VectReal_wp& teta, const VectReal_wp& phi, IVect& new_number)
  {
    int N = Points_ND.GetM();
    IVect number(N);
    int Nold = Points.GetM();
    number.Fill();

    VectReal_wp PointAdd(Points_ND);
    Sort(PointAdd, number);
    new_number.Reallocate(N);
    
    Points.Resize(Nold+N);
    TetaInterp.Resize(Nold + N);
    PhiInterp.Resize(Nold + N);
    for (int i = 0; i < number.GetM(); i++)
      {
        Points(i) = Points_ND(number(i));
        TetaInterp(i) = teta(number(i));
        PhiInterp(i) = phi(number(i));
        new_number(number(i)) = i;
      }
    
    for (int i = 0; i < number.GetM(); i++)
      new_number(i) += Nold;
  }

  
  //! an interpolation over a regular grid is computed
  /*!
    we need to know the element where is each interpolation point
    and its local coordinate on this element
  */
  void GridInterpolation<Dimension1>::LocalizePoints(Mesh<Dimension1>& mesh)
  {
    int N = Points.GetM();
    ElemInterp.Reallocate(N);
    CoorInterp.Reallocate(N);
    ElemInterp.Fill(-1);
    CoorInterp.Fill(0);
    
    if (mesh.GetNbVertices() <= 1)
      return;
   
    int ne = 0;
    Real_wp xe = mesh.Vertex(ne);
    // loop over points to find
    Real_wp Lmax = abs(mesh.Vertex(mesh.GetNbVertices()-1) - xe);
    for (int i = 0; i < N; i++)
      {	
	// we search the next point ne, so that Points(i) belongs to [vertex(ne-1), vertex(ne)]
	while ((ne < mesh.GetNbVertices())&&(Points(i) >= xe-100.0*epsilon_machine*Lmax))
	  {
	    ne++;
	    if (ne < mesh.GetNbVertices())
	      xe = mesh.Vertex(ne);
	  }
	// DISP(i); DISP(Points(i)); DISP(ne);
	
	if ((ne > 0)&&(ne < mesh.GetNbVertices()))
	  {
	    // the point is in the mesh
	    Real_wp lambda = (Points(i) - mesh.Vertex(ne-1))/(mesh.Vertex(ne) - mesh.Vertex(ne-1));
	    ElemInterp(i) = ne-1;
	    CoorInterp(i) = lambda;
	    // DISP(ElemInterp(i)); DISP(CoorInterp(i));
	  }
	
	if (ne == mesh.GetNbVertices())
	  {
	    if (abs(Points(i) - mesh.Vertex(ne-1)) <= 100.0*epsilon_machine*Lmax)
	      {
		ElemInterp(i) = ne-2;
		CoorInterp(i) = 1.0;
	      }
	    
	    ne++;
	  }
      }
  }
  
  
  /*************************************
   * GridInterpolationFull<Dimension1> *
   *************************************/
  

  //! not used
  void GridInterpolationFull<Dimension1>::
  SetInputData(const string& description_field,
	       const VectString& parameters)
  {
    int nb = 0;
    if (!description_field.compare("SismoPlaneAxi"))
      {
	if (parameters.GetM() <= 10)
	  {
	    cout << "In SetInputData of GridInterpolation" << endl;
	    cout << "SismoPlaneAxi needs more parameters, for instance :" << endl;
	    cout << "SismoPlaneAxi = xmin ymin zmin xmax ymax zmax x0 y0 z0 Nx Ny" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	type_output_file = PLANE_AXI;
	xmin = to_num<Real_wp>(parameters(nb++));
	ymin = to_num<Real_wp>(parameters(nb++));
	zmin = to_num<Real_wp>(parameters(nb++));
	xmax = to_num<Real_wp>(parameters(nb++));
	ymax = to_num<Real_wp>(parameters(nb++));
	zmax = to_num<Real_wp>(parameters(nb++));
	center(0) = to_num<Real_wp>(parameters(nb++));	
	center(1) = to_num<Real_wp>(parameters(nb++));
	center(2) = to_num<Real_wp>(parameters(nb++));
	nbPoints_x = to_num<int>(parameters(nb++));
	nbPoints_y = to_num<int>(parameters(nb++));
	nbPoints_z = 1;
      }
    else if (!description_field.compare("SismoPlane"))
      {
	if (parameters.GetM() <= 5)
	  {
	    cout << "In SetInputData of GridInterpolation" << endl;
	    cout << "SismoPlane needs more parameters, for instance :" << endl;
	    cout << "SismoPlane = xmin xmax ymin ymax Nx Ny" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	type_output_file = PLANE;
	xmin = to_num<Real_wp>(parameters(nb++));
	xmax = to_num<Real_wp>(parameters(nb++));
	ymin = to_num<Real_wp>(parameters(nb++));
	ymax = to_num<Real_wp>(parameters(nb++));
	nbPoints_x = to_num<int>(parameters(nb++));
	nbPoints_y = to_num<int>(parameters(nb++));
	nbPoints_z = 1;
      }
    else if (!description_field.compare("SismoGrille3D"))
      {
	if (parameters.GetM() <= 3)
	  {
	    cout << "In SetInputData of GridInterpolation" << endl;
	    cout << "SismoGrille3D needs more parameters, for instance :" << endl;
	    cout << "SismoGrille3D = AUTO Nx Ny Nz" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	type_output_file = VOLUME_AXI;
	if (!parameters(nb).compare("AUTO"))
	  {
	    grid_interval_to_be_computed = true;
	    nb++;
	  }
	else
	  {
	    if (parameters.GetM() <= 8)
	      {
		cout << "In SetInputData of GridInterpolation" << endl;
		cout << "SismoGrille3D needs more parameters, for instance :" << endl;
		cout << "SismoGrille3D = xmin xmax ymin ymax zmin zmax Nx Ny Nz" << endl;
		cout << "Current parameters are : " << endl << parameters << endl;
		abort();
	      }

	    xmin = to_num<Real_wp>(parameters(nb++));
	    xmax = to_num<Real_wp>(parameters(nb++));
	    ymin = to_num<Real_wp>(parameters(nb++));
	    ymax = to_num<Real_wp>(parameters(nb++));
	    zmin = to_num<Real_wp>(parameters(nb++));
	    zmax = to_num<Real_wp>(parameters(nb++));
	  }
        
	nbPoints_x = to_num<int>(parameters(nb++));
	nbPoints_y = to_num<int>(parameters(nb++));
	nbPoints_z = to_num<int>(parameters(nb++));
      }
    else if (!description_field.compare("SismoGrille"))
      {
	if (parameters.GetM() <= 6)
	  {
	    cout << "In SetInputData of GridInterpolation" << endl;
	    cout << "SismoGrille needs more parameters, for instance :" << endl;
	    cout << "SismoGrille = x0 y0 z0 AUTO Nx Ny Nz" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	// 3-D display for axisymmetric case
	type_output_file = THREE_PLANES_AXI;
	center(0) = to_num<Real_wp>(parameters(nb++));
	center(1) = to_num<Real_wp>(parameters(nb++));
	center(2) = to_num<Real_wp>(parameters(nb++));
	if (!parameters(nb).compare("AUTO"))
	  {
	    grid_interval_to_be_computed = true;
	    nb++;
	  }
	else
	  {
	    if (parameters.GetM() <= 11)
	      {
		cout << "In SetInputData of GridInterpolation" << endl;
		cout << "SismoGrille needs more parameters, for instance :" << endl;
		cout << "SismoGrille = x0 y0 z0 xmin xmax ymin ymax zmin zmax Nx Ny Nz" << endl;
		cout << "Current parameters are : " << endl << parameters << endl;
		abort();
	      }

	    xmin = to_num<Real_wp>(parameters(nb++));
	    xmax = to_num<Real_wp>(parameters(nb++));
	    ymin = to_num<Real_wp>(parameters(nb++));
	    ymax = to_num<Real_wp>(parameters(nb++));
	    zmin = to_num<Real_wp>(parameters(nb++));
	    zmax = to_num<Real_wp>(parameters(nb++));
	  }
	
	nbPoints_x = to_num<int>(parameters(nb++));
	nbPoints_y = to_num<int>(parameters(nb++));
	nbPoints_z = to_num<int>(parameters(nb++));
      }
    else if (!description_field.compare("SismoLine"))
      {
	if (parameters.GetM() <= 2)
	  {
	    cout << "In SetInputData of GridInterpolation" << endl;
	    cout << "SismoLine needs more parameters, for instance :" << endl;
	    cout << "SismoLine = xmin xmax Nx" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
	type_output_file = LINE;
	xmin = to_num<Real_wp>(parameters(nb++));
	xmax = to_num<Real_wp>(parameters(nb++));
	nbPoints_x = to_num<int>(parameters(nb++));
	nbPoints_y = 1;
      }
    else if (!description_field.compare("SismoLineAxi"))
      {
	if (parameters.GetM() <= 6)
	  {
	    cout << "In SetInputData of GridInterpolation" << endl;
	    cout << "SismoLineAxi needs more parameters, for instance :" << endl;
	    cout << "SismoLineAxi = xmin xmax ymin ymax zmin zmax N" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	type_output_file = LINE_AXI;
	xmin = to_num<Real_wp>(parameters(nb++));
	xmax = to_num<Real_wp>(parameters(nb++));
	ymin = to_num<Real_wp>(parameters(nb++));
	ymax = to_num<Real_wp>(parameters(nb++));
	zmin = to_num<Real_wp>(parameters(nb++));
	zmax = to_num<Real_wp>(parameters(nb++));
	nbPoints_x = to_num<int>(parameters(nb++));
	nbPoints_y = 1; nbPoints_z = 1;
      }
    else if (!description_field.compare("SismoCircleAxi"))
      {
	if (parameters.GetM() <= 8)
	  {
	    cout << "In SetInputData of GridInterpolation" << endl;
	    cout << "SismoCircleAxi needs more parameters, for instance :" << endl;
	    cout << "SismoCircleAxi = x0 y0 z0 nx ny nz rx ry N" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	type_output_file = CIRCLE_AXI;
	center(0) = to_num<Real_wp>(parameters(nb++));
	center(1) = to_num<Real_wp>(parameters(nb++));
	center(2) = to_num<Real_wp>(parameters(nb++));
        // normale to the plane where the circle lies
	xmin = to_num<Real_wp>(parameters(nb++));
	ymin = to_num<Real_wp>(parameters(nb++));
	zmin = to_num<Real_wp>(parameters(nb++));
	// radius of the circle
	xmax = to_num<Real_wp>(parameters(nb++));
	ymax = to_num<Real_wp>(parameters(nb++));	
	nbPoints_x = to_num<int>(parameters(nb++));
	nbPoints_y = 1;
      }
    else if (!description_field.compare("SismoPoint"))
      {
	if (parameters.GetM() <= 1)
	  {
	    cout << "In SetInputData of GridInterpolation" << endl;
	    cout << "SismoPoint needs more parameters, for instance :" << endl;
	    cout << "SismoPoint = x0 y0" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	type_output_file = POINT;
	xmin = to_num<Real_wp>(parameters(nb++));
	xmax = xmin;
      }
    else if (!description_field.compare("SismoPointAxi"))
      {
	if (parameters.GetM() <= 2)
	  {
	    cout << "In SetInputData of GridInterpolation" << endl;
	    cout << "SismoPointAxi needs more parameters, for instance :" << endl;
	    cout << "SismoPointAxi = x0 y0 z0" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
	type_output_file = POINT_AXI;
	xmin = to_num<Real_wp>(parameters(nb++));
	ymin = to_num<Real_wp>(parameters(nb++));
	zmin = to_num<Real_wp>(parameters(nb++));
	xmax = xmin; ymax = ymin; zmax = zmin;
	nbPoints_x = 1; nbPoints_y = 1; nbPoints_z = 1;
      }
    else if (!description_field.compare("SismoPointsFileAxi"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of GridInterpolation" << endl;
	    cout << "SismoPointsFileAxi needs more parameters, for instance :" << endl;
	    cout << "SismoPointsFileAxi = " << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	type_output_file = POINTS_FILE_AXI;
	PointsInputFile = parameters(0);
      }
  }
  
  
  //! initialization of 1-D output
  void GridInterpolationFull<Dimension1>
  ::Init(int type, const Real_wp& x0, const Real_wp& xN, int N)
  {
    this->type_output_file = type;
    xmin = x0;
    xmax = xN;
    nbPoints_x = N;
  }
  
  
  //! output on a 3-D plane
  void GridInterpolationFull<Dimension1>
  ::SetPlaneAxiOutput(const R3& ptO, const R3& ptA, const R3& ptB,
                      int nbx, int nby)
  {
    this->type_output_file = PLANE_AXI;
    this->xmin = ptO(0); this->ymin = ptO(1); this->zmin = ptO(2);
    this->xmax = ptA(0); this->ymax = ptA(2); this->zmax = ptA(2);
    this->center = ptB;
    this->nbPoints_x = nbx;
    this->nbPoints_y = nby;
    this->nbPoints_z = 1;
  }
  
  
  //! output on a 3-D line
  void GridInterpolationFull<Dimension1>
  ::SetLineAxiOutput(const R3& ptA, const R3& ptB, int nbx)
  {
    this->type_output_file = LINE_AXI;
    this->xmin = ptA(0); this->xmax = ptB(0);
    this->ymin = ptA(1); this->ymax = ptB(1);
    this->zmin = ptA(2); this->zmax = ptB(2);
    this->nbPoints_x = nbx;    
    this->nbPoints_y = 1;
    this->nbPoints_z = 1;
  }
  
  
  //! output on a 3-D point
  void GridInterpolationFull<Dimension1>::SetPointAxiOutput(const R3& ptA)
  {
    this->type_output_file = POINT_AXI;
    this->xmin = ptA(0); this->ymin = ptA(1); this->zmin = ptA(2);
    this->nbPoints_x = 1;
    this->nbPoints_y = 1;
    this->nbPoints_z = 1;
  }
    
  
  //! output on a 3-D circle
  void GridInterpolationFull<Dimension1>
  ::SetCircleAxiOutput(const R3& ptA, const R3& normale, 
                       const Real_wp& rx, const Real_wp& ry, int nbx)
  {
    this->type_output_file = CIRCLE_AXI;
    this->center = ptA;
    this->xmin = normale(0); this->ymin = normale(1); this->zmin = normale(2);
    this->xmax = rx; this->ymax = ry;
    this->nbPoints_x = nbx;
    this->nbPoints_y = 1;
    this->nbPoints_z = 1;
  }
  
  
  //! output on 3-D points
  void GridInterpolationFull<Dimension1>::SetPointsFileAxiOutput(const string& name)
  {
    this->type_output_file = POINTS_FILE_AXI;
    PointsInputFile = name;
  }
    
  
  //! output on three planes in 3-D
  void GridInterpolationFull<Dimension1>
  ::SetThreePlanesAxiOutput(const R3& ptO, const R3& ptMin, const R3& ptMax,
                            int nbx_, int nby_, int nbz_)
  {
    this->type_output_file = THREE_PLANES_AXI;
    this->center = ptO;
    this->xmin = ptMin(0); this->ymin = ptMin(1); this->zmin = ptMin(2);
    this->xmax = ptMax(0); this->ymax = ptMax(1); this->zmax = ptMax(2);
    nbPoints_x = nbx_; nbPoints_y = nby_; nbPoints_z = nbz_;
  }
  
  
  //! output on 3-D regular grid
  void GridInterpolationFull<Dimension1>
  ::SetVolumeAxiOutput(const R3& ptMin, const R3& ptMax,
                       int nbx_, int nby_, int nbz_)
  {
    this->type_output_file = VOLUME_AXI;
    this->xmin = ptMin(0); this->ymin = ptMin(1); this->zmin = ptMin(2);
    this->xmax = ptMax(0); this->ymax = ptMax(1); this->zmax = ptMax(2);
    nbPoints_x = nbx_; nbPoints_y = nby_; nbPoints_z = nbz_;
  }
  

  //! computes points of the grid (2-D outputs)
  void GridInterpolationFull<Dimension1>
  ::GenerateGridPoints(VectReal_wp& Points1D, VectReal_wp& TetaPoints)
  {
    VectR2 Points2D;
    if (type_output_file == PLANE)
      {
	// 2-D cartesian grid on a plane
	int nbPoints_grid = nbPoints_x*nbPoints_y;
	Real_wp step_x(0), step_y(0);
	if (nbPoints_x > 1)
	  step_x = (xmax-xmin)/(nbPoints_x-1);
	if (nbPoints_y > 1)
	  step_y = (ymax-ymin)/(nbPoints_y-1);

	Points2D.Reallocate(nbPoints_grid);
	R2 point;
	for (int i = 0; i < nbPoints_x; i++)
	  for (int j = 0; j < nbPoints_y; j++)
	    {
	      int nb = j*nbPoints_x + i;
	      point(0) = xmin + i*step_x;
	      point(1) = ymin + j*step_y;
	      Points2D(nb) = point;
	    }
      }
    else if (this->type_output_file == LINE)
      {
	// output on a 1-D line of extremities (xmin,ymin) and (xmax,ymax)
	Points1D.Reallocate(nbPoints_x);
	Real_wp step_x(0);
	if (nbPoints_x > 1)
          step_x = (xmax-xmin)/(nbPoints_x-1);
        
	for (int i = 0; i < nbPoints_x; i++)
          Points1D(i) = xmin + i*step_x;
	
        TetaPoints.Reallocate(nbPoints_x); TetaPoints.Zero();
      }

    if (Points2D.GetM() > 0)
      {
	TetaPoints.Reallocate(Points2D.GetM());
	Points1D.Reallocate(Points2D.GetM());
	for (int i = 0; i < Points2D.GetM(); i++)
	  {
	    Real_wp x = Points2D(i)(0), y = Points2D(i)(1);
	    Real_wp r, teta;
	    CartesianToPolar(x, y, r, teta);

	    Points1D(i) = r;
	    TetaPoints(i) = teta;
	  }
      }
  }
  

  //! computes points of the grid
  void GridInterpolationFull<Dimension1>
  ::GenerateGridPoints(VectReal_wp& Points1D, VectReal_wp& TetaPoints, VectReal_wp& PhiPoints)
  {
    VectR3 Points3D;
    if (type_output_file == PLANE_AXI)
      {
	// 2-D cartesian grid on a plane
	R3 vec_u, vec_v;
	vec_u(0) = xmax - xmin; vec_u(1) = ymax - ymin; vec_u(2) = zmax - zmin;
	vec_v(0) = center(0) - xmin; vec_v(1) = center(1) - ymin;
	vec_v(2) = center(2) - zmin;
	Real_wp step_x = 0.0, step_y = 0.0;
	if (nbPoints_x > 1)
	  step_x = 1.0/(nbPoints_x-1);
	if (nbPoints_y>1)
	  step_y = 1.0/(nbPoints_y-1);
	
	Points3D.Reallocate(nbPoints_x*nbPoints_y); R3 point;
	for (int i = 0; i < nbPoints_x; i++)
	  for (int j = 0; j < nbPoints_y; j++)
	    {
	      int nb = j*nbPoints_x + i;
	      point(0) = xmin + i*step_x*vec_u(0) + j*step_y*vec_v(0);
	      point(1) = ymin + i*step_x*vec_u(1) + j*step_y*vec_v(1);
	      point(2) = zmin + i*step_x*vec_u(2) + j*step_y*vec_v(2);
	      Points3D(nb) = point;
	    }
      }
    else if  (this->type_output_file == POINTS_FILE_AXI)
      {	
	xmin = 0; xmax = 1.0; nbPoints_x = 1;
	ymin = 0; ymax = 1.0; nbPoints_y = 1;
	zmin = 0; zmax = 1.0; nbPoints_z = 1;
	
	// we read 3-D points
	GlobalPoints3D.ReadText(PointsInputFile);
	Mlt(1.0/z0_adim, GlobalPoints3D);
	int nbPoints_grid = GlobalPoints3D.GetM();
	nbPoints_x = nbPoints_grid;	
        Points3D = GlobalPoints3D;
      }
    else if (type_output_file == THREE_PLANES_AXI)
      {
	// For an axisymmetric computation, output on 
	// three planes x = center(0), y = center(1) and z = center(2)
	// teta is computed for each point
	Real_wp step_x(0), step_y(0), step_z(0);
	if (nbPoints_x > 1)
	  step_x = (xmax-xmin)/(nbPoints_x-1);
	if (nbPoints_y > 1)
	  step_y = (ymax-ymin)/(nbPoints_y-1);
	if (nbPoints_z > 1)
	  step_z = (zmax-zmin)/(nbPoints_z-1);
    
	int nbPoints_grid = nbPoints_x*nbPoints_y + nbPoints_x*nbPoints_z + nbPoints_y*nbPoints_z;
	
	R3 point3D;
	Points3D.Reallocate(nbPoints_grid);
	for (int j = 0; j < nbPoints_y; j++)
	  for (int k = 0; k < nbPoints_z; k++)
	    {
	      int nb = j*nbPoints_z + k;
	      point3D(0) = center(0);
	      point3D(1) = ymin + j*step_y;
	      point3D(2) = zmin + k*step_z;
	      Points3D(nb) = point3D;
	    }
	
	for (int i = 0; i < nbPoints_x; i++)
	  for (int k = 0; k < nbPoints_z; k++)
	    {
	      int nb = nbPoints_z*nbPoints_y + i*nbPoints_z + k;
	      point3D(0) = xmin + i*step_x;
	      point3D(1) = center(1);
	      point3D(2) = zmin + k*step_z;
	      Points3D(nb) = point3D;
	    }
	
	for (int i = 0; i < nbPoints_x; i++)
	  for (int j = 0; j < nbPoints_y; j++)
	    {
	      int nb = nbPoints_z*(nbPoints_x+nbPoints_y) + i*nbPoints_y + j;
	      point3D(0) = xmin + i*step_x;
	      point3D(1) = ymin + j*step_y;
	      point3D(2) = center(2);
	      Points3D(nb) = point3D;
	    }
      }
    else if (type_output_file == VOLUME_AXI)
      {
        Real_wp step_x(0), step_y(0), step_z(0);
        if (nbPoints_x > 1)
          step_x = (xmax-xmin)/(nbPoints_x-1);
        if (nbPoints_y > 1)
          step_y = (ymax-ymin)/(nbPoints_y-1);
        if (nbPoints_z > 1)
          step_z = (zmax-zmin)/(nbPoints_z-1);

	// 3-D regular cartesian grid
	int nbPoints_grid = nbPoints_x*nbPoints_y*nbPoints_z;
	Points3D.Reallocate(nbPoints_grid); R3 point;
	for (int i = 0; i < nbPoints_x; i++)
	  for (int j = 0; j < nbPoints_y; j++)
	    for (int k = 0; k < nbPoints_z; k++)
	      {
		int nb = k*nbPoints_x*nbPoints_y + j*nbPoints_x + i;
		point(0) = xmin + i*step_x;
		point(1) = ymin + j*step_y;
		point(2) = zmin + k*step_z;
		Points3D(nb) = point;
	      }
      }
    else if (this->type_output_file == LINE)
      {
	// output on a 1-D line of extremities (xmin,ymin) and (xmax,ymax)
	Points1D.Reallocate(nbPoints_x);
	Real_wp step_x(0);
	if (nbPoints_x > 1)
          step_x = (xmax-xmin)/(nbPoints_x-1);
        
	for (int i = 0; i < nbPoints_x; i++)
          Points1D(i) = xmin + i*step_x;
	
        TetaPoints.Reallocate(nbPoints_x); TetaPoints.Fill(0);
        PhiPoints.Reallocate(nbPoints_x); PhiPoints.Fill(0);
      }
    else if (type_output_file == LINE_AXI)
      {
	// 1-D grid on a line
	Points3D.Reallocate(nbPoints_x);
	Real_wp step_x = 0.0, step_y = 0.0, step_z = 0.0; R3 point;
	if (nbPoints_x > 1)
	  {
	    step_x = (xmax-xmin)/(nbPoints_x-1);
	    step_y = (ymax-ymin)/(nbPoints_x-1);
	    step_z = (zmax-zmin)/(nbPoints_x-1);
	  }
	
	for (int i = 0; i < nbPoints_x; i++)
	  {
	    point(0) = xmin + i*step_x;
	    point(1) = ymin + i*step_y;
	    point(2) = zmin + i*step_z;
	    Points3D(i) = point;
	  }
      }
    else if (this->type_output_file == CIRCLE_AXI)
      {
        Points3D.Reallocate(nbPoints_x);
	Real_wp step_teta = 2.0*pi_wp/nbPoints_x;
	R3 normale, pt3D, vec_u, vec_v; Real_wp radius_x = xmax, radius_y = ymax;
	normale(0) = xmin; normale(1) = ymin; normale(2) = zmin;
	GetVectorPlane(normale, vec_u, vec_v);
	
	for (int i = 0; i < nbPoints_x; i++)
	  {
	    Real_wp teta = step_teta*i;
	    pt3D = center; Add(radius_x*cos(teta), vec_u, pt3D);
	    Add(radius_y*sin(teta), vec_v, pt3D);
	    Points3D(i) = pt3D;
	  }
      }
    else if (this->type_output_file == POINT)
      {
	// output on a single point
	Points1D.Reallocate(1);
	Points1D(0) = xmin;
	TetaPoints.Reallocate(1); TetaPoints(0) = 0;
	PhiPoints.Reallocate(1); PhiPoints(0) = 0;
       }
    else if (type_output_file == POINT_AXI)
      {
        Points3D.Reallocate(1);
        Points3D(0).Init(xmin, ymin, zmin);
      }        
    
    if (Points3D.GetM() > 0)
      {
        TetaPoints.Reallocate(Points3D.GetM());
        PhiPoints.Reallocate(Points3D.GetM());
	Points1D.Reallocate(Points3D.GetM());
	for (int i = 0; i < Points3D.GetM(); i++)
	  {
	    // r, teta, z from cartesian coordinates
	    Real_wp x = Points3D(i)(0), y = Points3D(i)(1), z = Points3D(i)(2), r, teta, phi;
            CartesianToSpherical(x, y, z, r, teta, phi);
            
	    Points1D(i) = r;
            PhiPoints(i) = phi;
	    TetaPoints(i) = teta;
	  }
      }

  }
  

  //! not used
  void GridInterpolationFull<Dimension1>::InitGrid(GridInterpolation<Dimension1>& grid,
						   int dim)
  {
    this->var_interp = &grid;
    if (this->grid_interval_to_be_computed)
      {
	xmin = xmin0; ymin = ymin0; zmin = zmin0;
	xmax = xmax0; ymax = ymax0; zmax = zmax0;
	
	if ( (this->type_output_file == THREE_PLANES_AXI)
             || (this->type_output_file == VOLUME_AXI) )
	  {
	    xmin = -xmax;
	    zmin = xmin; zmax = xmax;
	    ymin = xmin; ymax = xmax;
	  }
        
        xmin += 5.0*(xmax-xmin)*epsilon_machine;
        xmax -= 5.0*(xmax-xmin)*epsilon_machine;
        ymin += 5.0*(ymax-ymin)*epsilon_machine;
        ymax -= 5.0*(ymax-ymin)*epsilon_machine;
        zmin += 5.0*(zmax-zmin)*epsilon_machine;
        zmax -= 5.0*(zmax-zmin)*epsilon_machine;
      }
    
    if (dim == 3)
      {
	VectReal_wp Points1D, TetaPoints, PhiPoints;
	GenerateGridPoints(Points1D, TetaPoints, PhiPoints);
	
	grid.Append(Points1D, TetaPoints, PhiPoints, list_points);
      }
    else
      {
	VectReal_wp Points1D, TetaPoints;
	GenerateGridPoints(Points1D, TetaPoints);

	VectReal_wp PhiPoints(TetaPoints.GetM());
	PhiPoints.Zero();

	grid.Append(Points1D, TetaPoints, PhiPoints, list_points);
      }
  }
  

  //! writing header in binary
  template<class real>
  void GridInterpolationFull<Dimension1>::WriteBinary(ostream& FileStream, real& a) const
  {
    real tmp; int nb;
    switch (type_output_file)
      {
      case PLANE :
	// for 2-D regular grid, xmin, xmax, nbx, ymin, ymax, nby
	tmp = toDouble(xmin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(xmax*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&nbPoints_x)),sizeof(int));
	tmp = toDouble(ymin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp),sizeof(real));
	tmp = toDouble(ymax*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp),sizeof(real));
	FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&nbPoints_y)),sizeof(int));	
	break;
      case PLANE_AXI :
	// plane OAB with O = (xmin, ymin, zmin)
	// A = (xmax, ymax, zmax), B = center
	tmp = toDouble(xmin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(xmax*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(ymin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(ymax*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(zmin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(zmax*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));

	tmp = toDouble(center(0)*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(center(1)*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(center(2)*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	nb = nbPoints_x; FileStream.write(reinterpret_cast<char*>(&nb), sizeof(int));
	nb = nbPoints_y; FileStream.write(reinterpret_cast<char*>(&nb), sizeof(int));
	break;
      case VOLUME_AXI :
	tmp = toDouble(xmin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(xmax*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	nb = nbPoints_x; FileStream.write(reinterpret_cast<char*>(&nb), sizeof(int));
	tmp = toDouble(ymin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(ymax*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	nb = nbPoints_y; FileStream.write(reinterpret_cast<char*>(&nb), sizeof(int));
	tmp = toDouble(zmin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(zmax*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	nb = nbPoints_z; FileStream.write(reinterpret_cast<char*>(&nb), sizeof(int));
	break;
      case LINE :
        tmp = toDouble(xmin); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
        tmp = toDouble(xmax); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
        FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&nbPoints_x)), sizeof(int));    
	break;
      case LINE_AXI :
	// line joining A=(xmin,ymin,zmin) and B = (xmax, ymax, zmax)
	tmp = toDouble(xmin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(xmax*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(ymin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(ymax*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(zmin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
        tmp = toDouble(zmax*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	nb = nbPoints_x; FileStream.write(reinterpret_cast<char*>(&nb), sizeof(int));
	break;
      case POINT_AXI :
	tmp = toDouble(xmin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(ymin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(zmin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	break;
      case THREE_PLANES_AXI :
	// for 2-D regular grid, xmin, xmax, nbx, ymin, ymax, nby
	tmp = toDouble(xmin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(xmax*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&nbPoints_x)),sizeof(int));
	tmp = toDouble(ymin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp),sizeof(real));
	tmp = toDouble(ymax*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp),sizeof(real));
	FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&nbPoints_y)),sizeof(int));	
	tmp = toDouble(zmin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp),sizeof(real));
	tmp = toDouble(zmax*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp),sizeof(real));
	FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&nbPoints_z)),sizeof(int));	
	tmp = toDouble(center(0)*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp),sizeof(real));
	tmp = toDouble(center(1)*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp),sizeof(real));
	tmp = toDouble(center(2)*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp),sizeof(real));
	break;
      case CIRCLE_AXI :
	tmp = toDouble(center(0)*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(center(1)*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(center(2)*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(xmin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(ymin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(zmin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(xmax*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(ymax*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&nbPoints_x)),sizeof(int));
	break;
      case POINTS_FILE_AXI :
	{
	  nb = GlobalPoints3D.GetM();
	  FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&nb)),sizeof(int));
	  Vector<real> Pts3D(3*nb);
	  for (int i = 0; i < nb; i++)
	    {
	      Pts3D(3*i) = toDouble(GlobalPoints3D(i)(0)*z0_adim); 
	      Pts3D(3*i+1) = toDouble(GlobalPoints3D(i)(1)*z0_adim);
	      Pts3D(3*i+2) = toDouble(GlobalPoints3D(i)(2)*z0_adim);
	    }
	  FileStream.write(reinterpret_cast<char*>(const_cast<real*>(Pts3D.GetData())),
			   3*nb*sizeof(real));
	}
	break;
      }

  }
  
    
  //! displays details of class GridInterpolationFull<Dimension1>
  ostream& operator<<(ostream& out, const GridInterpolationFull<Dimension1>& g)
  {
    out << "Number of points = " << g.GetNbPoints() << '\n';
    return out;
  }
#endif
  
  
  /*********************************
   * GridInterpolation<Dimension2> *
   *********************************/
  
  
#ifdef MONTJOIE_WITH_TWO_DIM
  
  //! initialization in order to localize each element on a regular grid
  /*!
    A regular grid is split in Nx x Ny cells
    where Nx = nb_subdiv_grid_x and Ny = nb_subdiv_grid_y
    on each element i of the mesh, we retrieve all the intersecting
    cells of the regular grid with this element
    These cells are stored in the array ListeBoxGrid_Element(i)
   */
  void GridInterpolation<Dimension2>::
  InitInterpolationGrid(const Mesh<Dimension2>& mesh)
  {
    // s will be list of vertices, box the straight polygon bounding element
    // enveloppe the rectangle bounding element
    VectR2 s, box(4), res(2);
    TinyVector<R2, 2> enveloppe;
    Matrix2_2 dfj;
    SetPoints<Dimension2> PointsElem;
    
    if (nb_subdiv_grid_x <= 0)
      {
	// if no default grid is set, we create a default grid
	xmin = mesh.GetXmin(); ymin = mesh.GetYmin();
	xmax = mesh.GetXmax(); ymax = mesh.GetYmax();
	
	Real_wp deltax = xmax-xmin, deltay = ymax-ymin,
	  delta_glob = abs(deltax)+abs(deltay);
	
	// grid step along x, y
	Real_wp step_x(0), step_y(0);
	
	int ndiv = toInteger(sqrt(Real_wp(mesh.GetNbElt())));
	int nbx = 2*max(toInteger(ceil(deltax/delta_glob*ndiv)), 1);
	int nby = 2*max(toInteger(ceil(deltay/delta_glob*ndiv)), 1);
	if (nbx >= 1)
	  step_x = deltax/nbx;
	if (nby >= 1)
	  step_y = deltay/nby;
	
	nb_subdiv_grid_x = nbx; nb_subdiv_grid_y = nby;
	step_subdiv_grid_x = step_x; step_subdiv_grid_y = step_y;
      }

    Real_wp step_x = step_subdiv_grid_x, step_y = step_subdiv_grid_y;
    int nbx = nb_subdiv_grid_x, nby = nb_subdiv_grid_y;

    // filling the main array ListeBoxGrid_Element
    ListeBoxGrid_Element.Reallocate(mesh.GetNbElt());
    // loop on each element
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
	// we get vertices of element i
	mesh.GetVerticesElement(i, s);
	// computation of nodal points on the real mesh
	mesh.FjElemNodal(s, PointsElem, mesh, i);
	
	// bounding rectangle in enveloppe and bounding polygon in box
	mesh.GetBoundingBox(s, mesh.GetSafetyCoef(i), box, enveloppe);
	
	// interval along y, where the interpolation points of the grid can be
	int iymin = 0, iymax = nby-1;
	if (abs(step_y) > 0)
	  iymin = max(toInteger(floor((enveloppe(0)(1)-ymin)/step_y)), 0);
	
	if (abs(step_y) > 0)
	  iymax = min(toInteger(ceil((enveloppe(1)(1)-ymin)/step_y)), nby-1);
	
	// loop on y
	for (int iy = iymin; iy <= iymax; iy++)
	  {
	    int ixmin = 0, ixmax = nbx-1;
	    // interval along x, where the interpolation points of the grid can be
	    if (abs(step_x) > 0)
	      ixmin = max(toInteger(floor((enveloppe(0)(0)-xmin)/step_x)), 0);
	    
	    if (abs(step_x) > 0)
	      ixmax = min(toInteger(ceil((enveloppe(1)(0)-xmin)/step_x)),nbx-1);
	    
	    // loop on x
	    for (int ix = ixmin; ix <= ixmax; ix++)
	      {
		int num_box = iy*nbx+ix;
		ListeBoxGrid_Element(i).PushBack(num_box);
	      }
	  }
      }
  }
  
  
  //! pre-localization of points on the regular grid
  /*!
    \param[out] NumBoxGrid_Point grid cell where is each point
    \param[out] ListPoints_Grid for each grid cell, list of points inside
    \param[out] NbPoints_Grid for each grid cell, number of points inside
    \param[in] compute_list if true, ListPoints_Grid and NbPoints_Grid are computed
   */
  void GridInterpolation<Dimension2>::
  PreLocalizePoints(IVect& NumBoxGrid_Point, Vector<IVect>& ListPoints_Grid,
		    IVect& NbPoints_Grid, bool compute_list)
  {
    int N = GlobalCoord.GetM();
    VectR2& Points2D = GlobalCoord;

    // points are localized on a regular grid
    int nb_subdiv_grid = nb_subdiv_grid_x*nb_subdiv_grid_y;
    NumBoxGrid_Point.Reallocate(N); NumBoxGrid_Point.Fill(-1);
    if (compute_list)
      {
	ListPoints_Grid.Reallocate(nb_subdiv_grid);
	NbPoints_Grid.Reallocate(nb_subdiv_grid); NbPoints_Grid.Zero();
      }
    
    for (int i = 0; i < Points2D.GetM(); i++)
      {
	Real_wp x = Points2D(i)(0), y = Points2D(i)(1);
	int ix = 0, iy = 0;
	ix = min(toInteger(floor((x-xmin)/step_subdiv_grid_x)), nb_subdiv_grid_x-1);
	iy = min(toInteger(floor((y-ymin)/step_subdiv_grid_y)), nb_subdiv_grid_y-1);
	      
	ix = max(0,ix); iy = max(0,iy);	
	NumBoxGrid_Point(i) = iy*nb_subdiv_grid_x + ix;
	if (compute_list)
	  NbPoints_Grid(NumBoxGrid_Point(i))++;
	
      }
    
    if (compute_list)
      for (int i = 0; i < nb_subdiv_grid; i++)
	{
	  ListPoints_Grid(i).Reallocate(NbPoints_Grid(i));
	  ListPoints_Grid(i).Fill(-1);
	}
    
    if (compute_list)
      NbPoints_Grid.Zero();
    
    if (compute_list)
      for (int i = 0; i < Points2D.GetM(); i++)
	{
	  int num_box = NumBoxGrid_Point(i);
	  if (num_box >= 0)
	    ListPoints_Grid(num_box)(NbPoints_Grid(num_box)++) = i;
	}
  }
  
  
  //! localization of points on the grid, using pre-computed localization of elements
  /*!
    \param[in] mesh mesh where the points GlobalCoord are localized
    It is supposed that a call of InitInterpolationGrid has been done before 
    in order to have a localization of elements on a regular grid
    For each point i of the array GlobalCoord, this method
    will compute 
    - the element where this point lies : ElementInterp(i)
    - the local coordinates of the point inside the element
    \f$ (\hat{x},\hat{y}) \f$ : CoorInterp(i)
    - if asked, transformation DF_i^{-1} at this point : dfjm1(i)
   */
  void GridInterpolation<Dimension2>::
  LocalizePoints(const Mesh<Dimension2>& mesh)
  {
    LocalizePoints(mesh, this->GlobalCoord, this->ElementInterp, this->CoorInterp, 
		   this->Epsilon_InsideElement, this->point_treated, this->dfjm1, this->dfjm1_store);
  }

  
  //! Localization of points2D in the mesh
  void GridInterpolation<Dimension2>::
  LocalizePoints(const Mesh<Dimension2>& mesh, const VectR2& Points2D,
		 IVect& ElementInt, VectR2& CoorInt, VectReal_wp& EpsilonInt,
		 Vector<bool>& point_treat, Vector<Matrix2_2>& mat_dfjm1, bool store_df)
  {
    int N = Points2D.GetM();
    
    // we initialize the array, which will contain number of elements
    // where each interpolation is.
    ElementInt.Resize(N);
    
    // initialization of local coordinates also
    CoorInt.Resize(N);
    EpsilonInt.Resize(N);
    EpsilonInt.Fill(-1);
    if (store_df)
      mat_dfjm1.Resize(N);
        
    // resizing array point_treated
    int Nold = point_treat.GetM();
    point_treat.Resize(N);
    for (int i = Nold; i < N; i++)
      point_treat(i) = false;
    
    // for points to localize, initialization of ElementInterp and CoorInterp
    for (int i = 0; i < N; i++)
      if (!point_treat(i))
	{
	  ElementInt(i) = -1;
	  CoorInt(i).Zero();
	}
    
    // points are localized on the regular grid
    // NumBoxGrid_Point(i) is the number of the cell of the regular grid 
    // where the point i lies
    // NbPoints_Grid(i) is the number of points which belong to the cell i
    // ListPoints_Grid(i) is the list of points which belong to the cell i
    int nb_subdiv_grid = this->nb_subdiv_grid_x*this->nb_subdiv_grid_y;
    IVect NumBoxGrid_Point(N); NumBoxGrid_Point.Fill(-1);
    Vector<IVect> ListPoints_Grid(nb_subdiv_grid);
    IVect NbPoints_Grid(nb_subdiv_grid); NbPoints_Grid.Zero();
    for (int i = 0; i < Points2D.GetM(); i++)
      if (!point_treat(i))
	{
	  Real_wp x = Points2D(i)(0), y = Points2D(i)(1);
	  int ix = 0, iy = 0;
	  if ((x > xmax + R2::threshold)||(x < xmin - R2::threshold)
	      ||(y > ymax + R2::threshold)||(y < ymin - R2::threshold))
	    {
	      NumBoxGrid_Point(i) = -1;
	    }
	  else
	    {
	      if (abs(this->step_subdiv_grid_x) > R2::threshold)
		ix = min(toInteger(floor((x-xmin)/this->step_subdiv_grid_x)),
			 this->nb_subdiv_grid_x-1);
	      
	      if (abs(this->step_subdiv_grid_y) > R2::threshold)
		iy = min(toInteger(floor((y-ymin)/this->step_subdiv_grid_y)),
			 this->nb_subdiv_grid_y-1);
	      
	      ix = max(0,ix); iy = max(0,iy);
	      NumBoxGrid_Point(i) = iy*this->nb_subdiv_grid_x + ix;
	      NbPoints_Grid(NumBoxGrid_Point(i))++;
	    }
	  
	  point_treat(i) = true;
	}
    
    for (int i = 0; i < nb_subdiv_grid; i++)
      {
	ListPoints_Grid(i).Reallocate(NbPoints_Grid(i));
	ListPoints_Grid(i).Fill(-1);
      }
    
    NbPoints_Grid.Zero();
    for (int i = 0; i < Points2D.GetM(); i++)
      {
	int num_box = NumBoxGrid_Point(i);
	if (num_box >= 0)
	  ListPoints_Grid(num_box)(NbPoints_Grid(num_box)++) = i;
      }


    R2 pointA, point_loc;
    Matrix2_2 dfj;

    int old_percent = 0, new_percent = 0;
    VectR2 box(4);
    TinyVector<R2, 2> enveloppe;
    int nb_points_tested = 0, nb_points_inverted = 0;
    // loop on elements
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
	// displays a progress bar
	new_percent = toInteger(round(Real_wp(i)/(mesh.GetNbElt()-1)*80));
        if (mesh.print_level >= 2)
          for (int percent = old_percent; percent < new_percent; percent++)
            { cout<<"#"; cout.flush(); }
	
	old_percent = new_percent;
	
	// first, we look for any point near the element
	bool near_point = false;
	for (int j = 0; j < ListeBoxGrid_Element(i).GetM(); j++)
	  if (NbPoints_Grid(ListeBoxGrid_Element(i)(j)) > 0)
	    near_point = true;	
	
	if (near_point)
	  {
	    // if there is a point to localize near the element i
	    // we try to localize ...
	    FjInverseProblem<Dimension2> inverseFj(mesh, i);
	    VectR2& s = inverseFj.GetVertices();
	    SetPoints<Dimension2>& PointsElem = inverseFj.GetSetPoints();
	    
	    // bounding rectangle in enveloppe and bounding polygon in box
	    mesh.GetBoundingBox(s, mesh.GetSafetyCoef(i), box, enveloppe);
	    
	    // loop on each box containing element i
	    for (int j = 0; j < ListeBoxGrid_Element(i).GetM(); j++)
	      {
		int num_box = ListeBoxGrid_Element(i)(j);
		for (int n = 0; n < NbPoints_Grid(num_box); n++)
		  {
		    nb_points_tested++;
		    // number of the interpolation point
		    int num_point = ListPoints_Grid(num_box)(n);
		    pointA = Points2D(num_point);		    
		    
		    // point inside the enveloppe ?
                    if (PointInsideBoundingBox(pointA, enveloppe))
                      {
                        nb_points_inverted++;
                        // Is this interpolation point inside the element i ?
                        // point_loc = F_i^{-1} (pointA) is computed 
                        bool test_inside = inverseFj.Solve(pointA, point_loc);
                        Real_wp distance_boundary_elt
                          = mesh.GetDistanceToBoundary(point_loc, i);
                        
                        if (distance_boundary_elt < 0)
                          mesh.ProjectPointOnBoundary(point_loc, i);
                        
                        if (distance_boundary_elt < EpsilonInt(num_point))
                          test_inside = false;
                        
                        if (test_inside)
                          {
                            // Yes, it is
                            // we update ElementInterp and CoorInterp
                            ElementInt(num_point) = i;
                            CoorInt(num_point) = point_loc;
                            EpsilonInt(num_point) = distance_boundary_elt;
                            if (store_df)
                              {
                                // and DF_i^{-1} if necessary
                                mesh.DFj(s, PointsElem, point_loc, dfj, mesh, i);
                                GetInverse(dfj, mat_dfjm1(num_point));
                              }
                          }
                      }		    
		  }
	      }
	  }
      }

    if (mesh.print_level >= 2)
      {
        cout << endl;
        cout << "Number of tested points : " << nb_points_tested << endl;
        cout << "Number of points for which F_i has been inverted : " << nb_points_inverted << endl;
        cout << "Average number of iterations for non-linear solver : "
             << Real_wp(FjInverseProblem<Dimension2>::nb_iter_all) / nb_points_inverted << endl;
      }    
  }
  
  
  //! displays details of class GridInterpolation<Dimension2>
  template<class Dimension>
  ostream& operator<<(ostream& out, const GridInterpolation<Dimension2>& g)
  {
    out<<static_cast<const GridInterpolation_Base<Dimension2>& >(g);
    return out;
  }
  
  
  /*************************************
   * GridInterpolationFull<Dimension2> *
   *************************************/
  
  
  //! reading of a line of the data file
  /*!
    \param[in] description_field keyword of the line
    \param[out] parameters values of the line
    \param[in] nb_param number of values
    this function fills correct parameters
    depending on the matching line of the data file
   */
  void GridInterpolationFull<Dimension2>::
  SetInputData(const string& description_field,
	       const VectString& parameters)
  {
    int nb = 0;
    if (!description_field.compare("SismoPlane"))
      {
	if (parameters.GetM() <= 2)
	  {
	    cout << "In SetInputData of GridInterpolation" << endl;
	    cout << "SismoPlane needs more parameters, for instance :" << endl;
	    cout << "SismoPlane = AUTO Nx Ny" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	type_output_file = PLANE;
	if (!parameters(nb).compare("AUTO"))
	  {
	    grid_interval_to_be_computed = true;
	    nb++;
	  }
	else
	  {
	    if (parameters.GetM() <= 5)
	      {
		cout << "In SetInputData of GridInterpolation" << endl;
		cout << "SismoPlane needs more parameters, for instance :" << endl;
		cout << "SismoPlane = xmin xmax ymin ymax Nx Ny" << endl;
		cout << "Current parameters are : " << endl << parameters << endl;
		abort();
	      }
	    
	    xmin = to_num<Real_wp>(parameters(nb++));
	    xmax = to_num<Real_wp>(parameters(nb++));
	    ymin = to_num<Real_wp>(parameters(nb++));
	    ymax = to_num<Real_wp>(parameters(nb++));
	  }
	
	nbPoints_x = to_num<int>(parameters(nb++));
	nbPoints_y = to_num<int>(parameters(nb++));
      }
    else if (!description_field.compare("SismoPlaneAxi"))
      {
	if (parameters.GetM() <= 10)
	  {
	    cout << "In SetInputData of GridInterpolation" << endl;
	    cout << "SismoPlaneAxi needs more parameters, for instance :" << endl;
	    cout << "SismoPlaneAxi = xmin ymin zmin xmax ymax zmax x0 y0 z0 Nx Ny" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	type_output_file = PLANE_AXI;
	xmin = to_num<Real_wp>(parameters(nb++));
	ymin = to_num<Real_wp>(parameters(nb++));
	zmin = to_num<Real_wp>(parameters(nb++));
	xmax = to_num<Real_wp>(parameters(nb++));
	ymax = to_num<Real_wp>(parameters(nb++));
	zmax = to_num<Real_wp>(parameters(nb++));
	center(0) = to_num<Real_wp>(parameters(nb++));
	center(1) = to_num<Real_wp>(parameters(nb++));
	center(2) = to_num<Real_wp>(parameters(nb++));
	nbPoints_x = to_num<int>(parameters(nb++));
	nbPoints_y = to_num<int>(parameters(nb++));
	nbPoints_z = 1;
      }
    else if (!description_field.compare("SismoGrille3D"))
      {
	if (parameters.GetM() <= 3)
	  {
	    cout << "In SetInputData of GridInterpolation" << endl;
	    cout << "SismoGrille3D needs more parameters, for instance :" << endl;
	    cout << "SismoGrille3D = AUTO Nx Ny Nz" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	type_output_file = VOLUME_AXI;
	if (!parameters(nb).compare("AUTO"))
	  {
	    grid_interval_to_be_computed = true;
	    nb++;
	  }
	else
	  {
	    if (parameters.GetM() <= 8)
	      {
		cout << "In SetInputData of GridInterpolation" << endl;
		cout << "SismoGrille3D needs more parameters, for instance :" << endl;
		cout << "SismoGrille3D = xmin xmax ymin ymax zmin zmax Nx Ny Nz" << endl;
		cout << "Current parameters are : " << endl << parameters << endl;
		abort();
	      }

	    xmin = to_num<Real_wp>(parameters(nb++));
	    xmax = to_num<Real_wp>(parameters(nb++));
	    ymin = to_num<Real_wp>(parameters(nb++));
	    ymax = to_num<Real_wp>(parameters(nb++));
	    zmin = to_num<Real_wp>(parameters(nb++));
	    zmax = to_num<Real_wp>(parameters(nb++));
	  }
        
	nbPoints_x = to_num<int>(parameters(nb++));
	nbPoints_y = to_num<int>(parameters(nb++));
	nbPoints_z = to_num<int>(parameters(nb++));
      }
    else if (!description_field.compare("SismoGrille"))
      {
	if (parameters.GetM() <= 6)
	  {
	    cout << "In SetInputData of GridInterpolation" << endl;
	    cout << "SismoGrille needs more parameters, for instance :" << endl;
	    cout << "SismoGrille = x0 y0 z0 AUTO Nx Ny Nz" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	// 3-D display for axisymmetric case
	type_output_file = THREE_PLANES_AXI;
	center(0) = to_num<Real_wp>(parameters(nb++));
	center(1) = to_num<Real_wp>(parameters(nb++));
	center(2) = to_num<Real_wp>(parameters(nb++));
	if (!parameters(nb).compare("AUTO"))
	  {
	    grid_interval_to_be_computed = true;
	    nb++;
	  }
	else
	  {
	    if (parameters.GetM() <= 11)
	      {
		cout << "In SetInputData of GridInterpolation" << endl;
		cout << "SismoGrille needs more parameters, for instance :" << endl;
		cout << "SismoGrille = x0 y0 z0 xmin xmax ymin ymax zmin zmax Nx Ny Nz" << endl;
		cout << "Current parameters are : " << endl << parameters << endl;
		abort();
	      }

	    xmin = to_num<Real_wp>(parameters(nb++));
	    xmax = to_num<Real_wp>(parameters(nb++));
	    ymin = to_num<Real_wp>(parameters(nb++));
	    ymax = to_num<Real_wp>(parameters(nb++));
	    zmin = to_num<Real_wp>(parameters(nb++));
	    zmax = to_num<Real_wp>(parameters(nb++));
	  }
	
	nbPoints_x = to_num<int>(parameters(nb++));
	nbPoints_y = to_num<int>(parameters(nb++));
	nbPoints_z = to_num<int>(parameters(nb++));
      }
    else if (!description_field.compare("SismoLine"))
      {
	if (parameters.GetM() <= 4)
	  {
	    cout << "In SetInputData of GridInterpolation" << endl;
	    cout << "SismoLine needs more parameters, for instance :" << endl;
	    cout << "SismoLine = xmin xmax ymin ymax Nx" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
	type_output_file = LINE;
	xmin = to_num<Real_wp>(parameters(nb++));
	xmax = to_num<Real_wp>(parameters(nb++));
	ymin = to_num<Real_wp>(parameters(nb++));
	ymax = to_num<Real_wp>(parameters(nb++));
	nbPoints_x = to_num<int>(parameters(nb++));
	nbPoints_y = 1;
      }
    else if (!description_field.compare("SismoLineAxi"))
      {
	if (parameters.GetM() <= 6)
	  {
	    cout << "In SetInputData of GridInterpolation" << endl;
	    cout << "SismoLineAxi needs more parameters, for instance :" << endl;
	    cout << "SismoLineAxi = xmin xmax ymin ymax zmin zmax N" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	type_output_file = LINE_AXI;
	xmin = to_num<Real_wp>(parameters(nb++));
	xmax = to_num<Real_wp>(parameters(nb++));
	ymin = to_num<Real_wp>(parameters(nb++));
	ymax = to_num<Real_wp>(parameters(nb++));
	zmin = to_num<Real_wp>(parameters(nb++));
	zmax = to_num<Real_wp>(parameters(nb++));
	nbPoints_x = to_num<int>(parameters(nb++));
	nbPoints_y = 1; nbPoints_z = 1;
      }
    else if (!description_field.compare("SismoCircle"))
      {
	if (parameters.GetM() <= 3)
	  {
	    cout << "In SetInputData of GridInterpolation" << endl;
	    cout << "SismoCircle needs more parameters, for instance :" << endl;
	    cout << "SismoCircle = x0 y0 r Nx" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
        
	type_output_file = CIRCLE;
	center(0) = to_num<Real_wp>(parameters(nb++));
	center(1) = to_num<Real_wp>(parameters(nb++));
	Real_wp radius = to_num<Real_wp>(parameters(nb++));
	nbPoints_x = to_num<int>(parameters(nb++));
	nbPoints_y = 1;
	xmin = center(0) - radius; xmax = center(0) + radius;
	ymin = center(1) - radius; ymax = center(1) + radius;
      }
    else if (!description_field.compare("SismoCircleAxi"))
      {
	if (parameters.GetM() <= 8)
	  {
	    cout << "In SetInputData of GridInterpolation" << endl;
	    cout << "SismoCircleAxi needs more parameters, for instance :" << endl;
	    cout << "SismoCircleAxi = x0 y0 z0 nx ny nz rx ry N" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	type_output_file = CIRCLE_AXI;
	center(0) = to_num<Real_wp>(parameters(nb++));
	center(1) = to_num<Real_wp>(parameters(nb++));
	center(2) = to_num<Real_wp>(parameters(nb++));
        // normale to the plane where the circle lies
	xmin = to_num<Real_wp>(parameters(nb++));
	ymin = to_num<Real_wp>(parameters(nb++));
	zmin = to_num<Real_wp>(parameters(nb++));
	// radius of the circle
	xmax = to_num<Real_wp>(parameters(nb++));
	ymax = to_num<Real_wp>(parameters(nb++));	
	nbPoints_x = to_num<int>(parameters(nb++));
	nbPoints_y = 1;
      }
    else if (!description_field.compare("SismoPoint"))
      {
	if (parameters.GetM() <= 1)
	  {
	    cout << "In SetInputData of GridInterpolation" << endl;
	    cout << "SismoPoint needs more parameters, for instance :" << endl;
	    cout << "SismoPoint = x0 y0" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	type_output_file = POINT;
	xmin = to_num<Real_wp>(parameters(nb++));
	ymin = to_num<Real_wp>(parameters(nb++));
	xmax = xmin; ymax = ymin;
      }
    else if (!description_field.compare("SismoPointAxi"))
      {
	if (parameters.GetM() <= 2)
	  {
	    cout << "In SetInputData of GridInterpolation" << endl;
	    cout << "SismoPointAxi needs more parameters, for instance :" << endl;
	    cout << "SismoPointAxi = x0 y0 z0" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }
	
	type_output_file = POINT_AXI;
	xmin = to_num<Real_wp>(parameters(nb++));
	ymin = to_num<Real_wp>(parameters(nb++));
	zmin = to_num<Real_wp>(parameters(nb++));
	xmax = xmin; ymax = ymin; zmax = zmin;
	nbPoints_x = 1; nbPoints_y = 1; nbPoints_z = 1;
      }
    else if (!description_field.compare("SismoPointsFileAxi"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of GridInterpolation" << endl;
	    cout << "SismoPointsFileAxi needs more parameters, for instance :" << endl;
	    cout << "SismoPointsFileAxi = " << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	type_output_file = POINTS_FILE_AXI;
	PointsInputFile = parameters(0);
      }
    else if (!description_field.compare("SismoPointsFile"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of GridInterpolation" << endl;
	    cout << "SismoPointsFile needs more parameters, for instance :" << endl;
	    cout << "SismoPointsFile = " << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	type_output_file = POINTS_FILE ;
	PointsInputFile = parameters(0);
      }
  }
  
  
  //! returns type of regular grid (plane, line, circle, etc)
  int GridInterpolationFull<Dimension2>::GetType(const string& description)
  {
    if (!description.compare("FileOutputPlane"))
      return PLANE;
    else if (!description.compare("FileOutputPlaneAxi"))
      return PLANE_AXI;
    else if (!description.compare("FileOutputLine"))
      return LINE;
    else if (!description.compare("FileOutputLineAxi"))
      return LINE_AXI;
    else if (!description.compare("FileOutputPoint"))
      return POINT;
    else if (!description.compare("FileOutputPointAxi"))
      return POINT_AXI;
    else if (!description.compare("FileOutputGrille"))
      return THREE_PLANES_AXI;
    else if (!description.compare("FileOutputGrille3D"))
      return VOLUME_AXI;
    else if (!description.compare("FileOutputCircle"))
      return CIRCLE;
    else if (!description.compare("FileOutputCircleAxi"))
      return CIRCLE_AXI;
    else if (!description.compare("FileOutputPointsFile"))
      return POINTS_FILE;
    else if (!description.compare("FileOutputPointsFileAxi"))
      return POINTS_FILE_AXI;
    
    return -1;
  }
  
  
  //! returns the dimension of the output grid (3 for axisymmetric outputs)
  int GridInterpolationFull<Dimension2>::GetDimension(int type)
  {
    switch (type)
      {
      case PLANE_AXI :
      case LINE_AXI :
      case POINT_AXI :
      case THREE_PLANES_AXI :
      case VOLUME_AXI :
      case CIRCLE_AXI :
      case POINTS_FILE_AXI :
        return 3;
      }
    
    return 2;
  }
  
  
  //! initializes grid with regular points on a rectangle
  void GridInterpolationFull<Dimension2>::
  SetPlaneOutput(const Real_wp& xmin_, const Real_wp& xmax_,
		 const Real_wp& ymin_, const Real_wp& ymax_,
		 int nbx, int nby)
  {
    this->type_output_file = PLANE;
    this->xmin = xmin_; this->xmax = xmax_;
    this->ymin = ymin_; this->ymax = ymax_;
    this->nbPoints_x = nbx;
    this->nbPoints_y = nby;
    this->nbPoints_z = 1;
  }

  
  //! initializes grid with regular points on a 3-D plane
  void GridInterpolationFull<Dimension2>::
  SetPlaneAxiOutput(const R3& ptO, const R3& ptA, const R3& ptB,
		    int nbx, int nby)
  {
    this->type_output_file = PLANE_AXI;
    this->xmin = ptO(0); this->ymin = ptO(1); this->zmin = ptO(2);
    this->xmax = ptA(0); this->ymax = ptA(2); this->zmax = ptA(2);
    this->center = ptB;
    this->nbPoints_x = nbx;
    this->nbPoints_y = nby;
    this->nbPoints_z = 1;
  }
  
  
  //! initialization of the grid with regular points on a line
  void GridInterpolationFull<Dimension2>::
  SetLineOutput(const Real_wp& xmin_, const Real_wp& xmax_,
		const Real_wp& ymin_, const Real_wp& ymax_, int nbx)
  {
    this->type_output_file = LINE;
    this->xmin = xmin_; this->xmax = xmax_;
    this->ymin = ymin_; this->ymax = ymax_;
    this->nbPoints_x = nbx;    
    this->nbPoints_y = 1;
    this->nbPoints_z = 1;
  }


  //! initialization of the grid with regular points on a 3-D line
  void GridInterpolationFull<Dimension2>::
  SetLineAxiOutput(const R3& ptA, const R3& ptB, int nbx)
  {
    this->type_output_file = LINE_AXI;
    this->xmin = ptA(0); this->xmax = ptB(0);
    this->ymin = ptA(1); this->ymax = ptB(1);
    this->zmin = ptA(2); this->zmax = ptB(2);
    this->nbPoints_x = nbx;    
    this->nbPoints_y = 1;
    this->nbPoints_z = 1;
  }
  
  
  //! initialization of the grid with a single 2-D point
  void GridInterpolationFull<Dimension2>::
  SetPointOutput(const R2& ptA)
  {
    this->type_output_file = POINT;
    this->xmin = ptA(0); this->ymin = ptA(1);
    this->nbPoints_x = 1;
    this->nbPoints_y = 1;
    this->nbPoints_z = 1;
  }


  //! initialization of the grid with a single 3-D point
  void GridInterpolationFull<Dimension2>::
  SetPointAxiOutput(const R3& ptA)
  {
    this->type_output_file = POINT_AXI;
    this->xmin = ptA(0); this->ymin = ptA(1); this->zmin = ptA(2);
    this->nbPoints_x = 1;
    this->nbPoints_y = 1;
    this->nbPoints_z = 1;
  }
  

  //! initialization of the grid with regular points on a circle
  void GridInterpolationFull<Dimension2>::
  SetCircleOutput(const R2& ptA, const Real_wp& radius, int nbx)
  {
    this->type_output_file = CIRCLE;
    this->xmin = ptA(0) - radius; this->ymin = ptA(1) - radius;
    this->xmax = ptA(0) + radius; this->ymax = ptA(1) + radius;
    this->nbPoints_x = nbx;
    this->nbPoints_y = 1;
    this->nbPoints_z = 1;
  }


  //! initialization of the grid with regular points on a 3-D circle
  void GridInterpolationFull<Dimension2>::
  SetCircleAxiOutput(const R3& ptA, const R3& normale,
		     const Real_wp& rx, const Real_wp& ry, int nbx)
  {
    this->type_output_file = CIRCLE_AXI;
    this->center = ptA;
    this->xmin = normale(0); this->ymin = normale(1); this->zmin = normale(2);
    this->xmax = rx; this->ymax = ry;
    this->nbPoints_x = nbx;
    this->nbPoints_y = 1;
    this->nbPoints_z = 1;
  }

  
  //! initialization of the grid with points stored in a file
  void GridInterpolationFull<Dimension2>::SetPointsFileOutput(const string& name)
  {
    this->type_output_file = POINTS_FILE;
    PointsInputFile = name;
  }
  
  
  //! initialization of the grid with 3-D points stored in a file
  void GridInterpolationFull<Dimension2>::SetPointsFileAxiOutput(const string& name)
  {
    this->type_output_file = POINTS_FILE_AXI;
    PointsInputFile = name;
  }
  
  
  //! initialization of the grid with regular points on three planes
  void GridInterpolationFull<Dimension2>::
  SetThreePlanesAxiOutput(const R3& ptO, const R3& ptMin, const R3& ptMax,
			  int nbx_, int nby_, int nbz_)
  {
    this->type_output_file = THREE_PLANES_AXI;
    this->center = ptO;
    this->xmin = ptMin(0); this->ymin = ptMin(1); this->zmin = ptMin(2);
    this->xmax = ptMax(0); this->ymax = ptMax(1); this->zmax = ptMax(2);
    nbPoints_x = nbx_; nbPoints_y = nby_; nbPoints_z = nbz_;
  }
  
  
  //! initialization of the grid with regular points on a parallelepiped
  void GridInterpolationFull<Dimension2>::
  SetVolumeAxiOutput(const R3& ptMin, const R3& ptMax,
		     int nbx_, int nby_, int nbz_)
  {
    this->type_output_file = VOLUME_AXI;
    this->xmin = ptMin(0); this->ymin = ptMin(1); this->zmin = ptMin(2);
    this->xmax = ptMax(0); this->ymax = ptMax(1); this->zmax = ptMax(2);
    nbPoints_x = nbx_; nbPoints_y = nby_; nbPoints_z = nbz_;
  }


  //! computation of 2-D points (and theta-coordinate for axisymmetric outputs)
  void GridInterpolationFull<Dimension2>::
  GenerateGridPoints(VectR2& Points2D, VectReal_wp& TetaPoints)
  {
    GridInterpolation<Dimension2> grid; IVect Index, NumPartie;
    grid.SetXmin(-1e300); grid.SetXmax(1e300);
    grid.SetYmin(-1e300); grid.SetYmax(1e300);
    grid.SetZmin(-1e300); grid.SetZmax(1e300);

    VectR2 PointsInit;
    GenerateGridPoints(grid, Index, PointsInit, TetaPoints, NumPartie);

    Points2D.Reallocate(PointsInit.GetM());
    for (int i = 0; i < Points2D.GetM(); i++)
      Points2D(Index(i)) = PointsInit(i);
  }

  
  //! computation of 2-D points (and theta-coordinate for axisymmetric outputs)
  void GridInterpolationFull<Dimension2>::
  GenerateGridPoints(const GridInterpolation<Dimension2>& grid, IVect& Index,
		     VectR2& Points2D, VectReal_wp& TetaPoints, IVect& NumPartie)
  {    
    R3 point; R2 pt2D; Real_wp theta;
    int n;
    if (this->type_output_file == PLANE)
      {
	// output on a 2-D cartesian grid
	int nbPoints_grid = nbPoints_x*nbPoints_y;
	Real_wp step_x(0), step_y(0);
	if (nbPoints_x > 1)
	  step_x = (xmax-xmin)/(nbPoints_x-1);
	if (nbPoints_y > 1)
	  step_y = (ymax-ymin)/(nbPoints_y-1);
	
	// counting the number of points
	int num = 0;
	for (int i = 0; i < nbPoints_x; i++)
	  for (int j = 0; j < nbPoints_y; j++)
	    {
	      point(0) = xmin + i*step_x;
	      point(1) = ymin + j*step_y;
	      if (TranslatePoint(point, false, grid, pt2D, theta, n))
		num++;
	    }
	
	nb_global_grid_points = nbPoints_grid;
	Points2D.Reallocate(num);
	Index.Reallocate(num);
	NumPartie.Reallocate(num);
	TetaPoints.Reallocate(num);
	num = 0;
	for (int i = 0; i < nbPoints_x; i++)
	  for (int j = 0; j < nbPoints_y; j++)
	    {
	      int nb = j*nbPoints_x + i;
	      point(0) = xmin + i*step_x;
	      point(1) = ymin + j*step_y;
	      if (TranslatePoint(point, false, grid, pt2D, theta, n))
		{
		  Index(num) = nb;
		  Points2D(num) = pt2D;
		  TetaPoints(num) = theta;
		  NumPartie(num) = n;
		  num++;
		}
	    }
      }
    else if (type_output_file == PLANE_AXI)
      {
	// 2-D cartesian grid on a plane
	R3 vec_u, vec_v;
	vec_u(0) = xmax - xmin; vec_u(1) = ymax - ymin; vec_u(2) = zmax - zmin;
	vec_v(0) = center(0) - xmin; vec_v(1) = center(1) - ymin;
	vec_v(2) = center(2) - zmin;
	Real_wp step_x = 0.0, step_y = 0.0;
	if (nbPoints_x > 1)
	  step_x = 1.0/(nbPoints_x-1);
	if (nbPoints_y>1)
	  step_y = 1.0/(nbPoints_y-1);
	
	int num = 0;
	for (int i = 0; i < nbPoints_x; i++)
	  for (int j = 0; j < nbPoints_y; j++)
	    {
	      point(0) = xmin + i*step_x*vec_u(0) + j*step_y*vec_v(0);
	      point(1) = ymin + i*step_x*vec_u(1) + j*step_y*vec_v(1);
	      point(2) = zmin + i*step_x*vec_u(2) + j*step_y*vec_v(2);
	      if (TranslatePoint(point, true, grid, pt2D, theta, n))
		num++;
	    }
	
	nb_global_grid_points = nbPoints_x*nbPoints_y;
	Index.Reallocate(num);
        TetaPoints.Reallocate(num);
	Points2D.Reallocate(num);
	NumPartie.Reallocate(num);
	num = 0;
	for (int i = 0; i < nbPoints_x; i++)
	  for (int j = 0; j < nbPoints_y; j++)
	    {
	      int nb = j*nbPoints_x + i;
	      point(0) = xmin + i*step_x*vec_u(0) + j*step_y*vec_v(0);
	      point(1) = ymin + i*step_x*vec_u(1) + j*step_y*vec_v(1);
	      point(2) = zmin + i*step_x*vec_u(2) + j*step_y*vec_v(2);
	      if (TranslatePoint(point, true, grid, pt2D, theta, n))
		{
		  Index(num) = nb;
		  Points2D(num) = pt2D;
		  TetaPoints(num) = theta;
		  NumPartie(num) = n;
		  num++;
		}
	    }
      }
    else if  (this->type_output_file == POINTS_FILE)
      {	
	xmin = 0; xmax = 1.0; nbPoints_x = 1;
	ymin = 0; ymax = 1.0; nbPoints_y = 1;
	
	GlobalPoints2D.ReadText(PointsInputFile);
	Mlt(1.0/z0_adim, GlobalPoints2D);
	
	int num = 0;
	for (int i = 0; i < GlobalPoints2D.GetM(); i++)
	  {
	    point.Init(GlobalPoints2D(i)(0), GlobalPoints2D(i)(1));
	    if (TranslatePoint(point, false, grid, pt2D, theta, n))
	      num++;
	  }
	
	nb_global_grid_points = GlobalPoints2D.GetM();
	nbPoints_x = nb_global_grid_points;

	Points2D.Reallocate(num);
	Index.Reallocate(num);
        TetaPoints.Reallocate(num);
	NumPartie.Reallocate(num);
	num = 0;
	for (int i = 0; i < GlobalPoints2D.GetM(); i++)
	  {
	    point.Init(GlobalPoints2D(i)(0), GlobalPoints2D(i)(1));
	    if (TranslatePoint(point, false, grid, pt2D, theta, n))
	      {
		Index(num) = i;
		Points2D(num) = pt2D;
		TetaPoints(num) = theta;
		NumPartie(num) = n;
		num++;
	      }
	  }
	
#ifdef SELDON_WITH_MPI
        int rank_proc; MPI_Comm_rank(grid.GetCommunicator(), &rank_proc);
	if (rank_proc != 0)
	  GlobalPoints2D.Clear();
#endif
      }
    else if  (this->type_output_file == POINTS_FILE_AXI)
      {	
	xmin = 0; xmax = 1.0; nbPoints_x = 1;
	ymin = 0; ymax = 1.0; nbPoints_y = 1;
	zmin = 0; zmax = 1.0; nbPoints_z = 1;
	
	// we read 3-D points
	GlobalPoints3D.ReadText(PointsInputFile);
	Mlt(1.0/z0_adim, GlobalPoints3D);
	
	int num = 0;
	for (int i = 0; i < GlobalPoints3D.GetM(); i++)
	  {
	    point = GlobalPoints3D(i);
	    if (TranslatePoint(point, true, grid, pt2D, theta, n))
	      num++;
	  }

	nb_global_grid_points = GlobalPoints3D.GetM();
	nbPoints_x = nb_global_grid_points;
	
	Points2D.Reallocate(num);
	Index.Reallocate(num);
        TetaPoints.Reallocate(num);
	NumPartie.Reallocate(num);
	num = 0;
	for (int i = 0; i < GlobalPoints3D.GetM(); i++)
	  {
	    point = GlobalPoints3D(i);
	    if (TranslatePoint(point, true, grid, pt2D, theta, n))
	      {
		Index(num) = i;
		Points2D(num) = pt2D;
		TetaPoints(num) = theta;
		NumPartie(num) = n;
		num++;
	      }
	  }

#ifdef SELDON_WITH_MPI
        int rank_proc; MPI_Comm_rank(grid.GetCommunicator(), &rank_proc);
	if (rank_proc != 0)
	  GlobalPoints3D.Clear();
#endif
      }
    else if (type_output_file == THREE_PLANES_AXI)
      {
	// For an axisymmetric computation, output on 
	// three planes x = center(0), y = center(1) and z = center(2)
	// teta is computed for each point
	Real_wp step_x(0), step_y(0), step_z(0);
	if (nbPoints_x>1)
	  step_x = (xmax-xmin)/(nbPoints_x-1);
	if (nbPoints_y>1)
	  step_y = (ymax-ymin)/(nbPoints_y-1);
	if (nbPoints_z>1)
	  step_z = (zmax-zmin)/(nbPoints_z-1);
    
	int nbPoints_grid = nbPoints_x*nbPoints_y + nbPoints_x*nbPoints_z
	  + nbPoints_y*nbPoints_z;
	
	int num = 0;
	for (int j = 0; j < nbPoints_y; j++)
	  for (int k = 0; k < nbPoints_z; k++)
	    {
	      point(0) = center(0);
	      point(1) = ymin + j*step_y;
	      point(2) = zmin + k*step_z;
	      if (TranslatePoint(point, true, grid, pt2D, theta, n))
		num++;
	    }
	
	for (int i = 0; i < nbPoints_x; i++)
	  for (int k = 0; k < nbPoints_z; k++)
	    {
	      point(0) = xmin + i*step_x;
	      point(1) = center(1);
	      point(2) = zmin + k*step_z;
	      if (TranslatePoint(point, true, grid, pt2D, theta, n))
		num++;
	    }
	
	for (int i = 0; i < nbPoints_x; i++)
	  for (int j = 0; j < nbPoints_y; j++)
	    {
	      point(0) = xmin + i*step_x;
	      point(1) = ymin + j*step_y;
	      point(2) = center(2);
	      if (TranslatePoint(point, true, grid, pt2D, theta, n))
		num++;
	    }

	nb_global_grid_points = nbPoints_grid;
	
	Points2D.Reallocate(num);
	Index.Reallocate(num);
        TetaPoints.Reallocate(num);
	NumPartie.Reallocate(num);
	num = 0;
	for (int j = 0; j < nbPoints_y; j++)
	  for (int k = 0; k < nbPoints_z; k++)
	    {
	      int nb = j*nbPoints_z + k;
	      point(0) = center(0);
	      point(1) = ymin + j*step_y;
	      point(2) = zmin + k*step_z;
	      if (TranslatePoint(point, true, grid, pt2D, theta, n))
		{
		  Index(num) = nb;
		  Points2D(num) = pt2D;
		  TetaPoints(num) = theta;
		  NumPartie(num) = n;
		  num++;
		}
	    }
	
	for (int i = 0; i < nbPoints_x; i++)
	  for (int k = 0; k < nbPoints_z; k++)
	    {
	      int nb = nbPoints_z*nbPoints_y + i*nbPoints_z + k;
	      point(0) = xmin + i*step_x;
	      point(1) = center(1);
	      point(2) = zmin + k*step_z;
	      if (TranslatePoint(point, true, grid, pt2D, theta, n))
		{
		  Index(num) = nb;
		  Points2D(num) = pt2D;
		  TetaPoints(num) = theta;
		  NumPartie(num) = n;
		  num++;
		}
	    }
	
	for (int i = 0; i < nbPoints_x; i++)
	  for (int j = 0; j < nbPoints_y; j++)
	    {
	      int nb = nbPoints_z*(nbPoints_x+nbPoints_y) + i*nbPoints_y + j;
	      point(0) = xmin + i*step_x;
	      point(1) = ymin + j*step_y;
	      point(2) = center(2);
	      if (TranslatePoint(point, true, grid, pt2D, theta, n))
		{
		  Index(num) = nb;
		  Points2D(num) = pt2D;
		  TetaPoints(num) = theta;
		  NumPartie(num) = n;
		  num++;
		}
	    }
      }
    else if (type_output_file == VOLUME_AXI)
      {
        Real_wp step_x(0), step_y(0), step_z(0);
        if (nbPoints_x > 1)
          step_x = (xmax-xmin)/(nbPoints_x-1);
        if (nbPoints_y > 1)
          step_y = (ymax-ymin)/(nbPoints_y-1);
        if (nbPoints_z > 1)
          step_z = (zmax-zmin)/(nbPoints_z-1);

	// 3-D regular cartesian grid
	int nbPoints_grid = nbPoints_x*nbPoints_y*nbPoints_z;
	int num = 0;
	for (int i = 0; i < nbPoints_x; i++)
	  for (int j = 0; j < nbPoints_y; j++)
	    for (int k = 0; k < nbPoints_z; k++)
	      {
		point(0) = xmin + i*step_x;
		point(1) = ymin + j*step_y;
		point(2) = zmin + k*step_z;
		if (TranslatePoint(point, true, grid, pt2D, theta, n))
		  num++;
	      }
	
	nb_global_grid_points = nbPoints_grid;
	
	Points2D.Reallocate(num);
	Index.Reallocate(num);
        TetaPoints.Reallocate(num);
	NumPartie.Reallocate(num);
	num = 0;
	
	for (int i = 0; i < nbPoints_x; i++)
	  for (int j = 0; j < nbPoints_y; j++)
	    for (int k = 0; k < nbPoints_z; k++)
	      {
		int nb = k*nbPoints_x*nbPoints_y + j*nbPoints_x + i;
		point(0) = xmin + i*step_x;
		point(1) = ymin + j*step_y;
		point(2) = zmin + k*step_z;
		if (TranslatePoint(point, true, grid, pt2D, theta, n))
		  {
		    Index(num) = nb;
		    Points2D(num) = pt2D;
		    TetaPoints(num) = theta;
		    NumPartie(num) = n;
		    num++;
		  }
	      }
      }
    else if (this->type_output_file == LINE)
      {
	// output on a 1-D line of extremities (xmin,ymin) and (xmax,ymax)
	Real_wp step_x(0), step_y(0);
	if (nbPoints_x > 1)
	  {
	    step_x = (xmax-xmin)/(nbPoints_x-1);
	    step_y = (ymax-ymin)/(nbPoints_x-1);
	  }
	
	int num = 0;
	for (int i = 0; i < nbPoints_x; i++)
	  {
	    point(0) = xmin + i*step_x;
	    point(1) = ymin + i*step_y;
	    if (TranslatePoint(point, false, grid, pt2D, theta, n))
	      num++;
	  }

	nb_global_grid_points = nbPoints_x;
	
	Points2D.Reallocate(num);
	Index.Reallocate(num);
        TetaPoints.Reallocate(num);
	NumPartie.Reallocate(num);
	num = 0;

	for (int i = 0; i < nbPoints_x; i++)
	  {
	    point(0) = xmin + i*step_x;
	    point(1) = ymin + i*step_y;
	    if (TranslatePoint(point, false, grid, pt2D, theta, n))
	      {
		Index(num) = i;
		Points2D(num) = pt2D;
		TetaPoints(num) = theta;
		NumPartie(num) = n;
		num++;
	      }
	  }
      }
    else if (type_output_file == LINE_AXI)
      {
	// 1-D grid on a line
	Real_wp step_x = 0.0, step_y = 0.0, step_z = 0.0;
	if (nbPoints_x > 1)
	  {
	    step_x = (xmax-xmin)/(nbPoints_x-1);
	    step_y = (ymax-ymin)/(nbPoints_x-1);
	    step_z = (zmax-zmin)/(nbPoints_x-1);
	  }
	
	int num = 0;
	for (int i = 0; i < nbPoints_x; i++)
	  {
	    point(0) = xmin + i*step_x;
	    point(1) = ymin + i*step_y;
	    point(2) = zmin + i*step_z;
	    if (TranslatePoint(point, true, grid, pt2D, theta, n))
	      num++;
	  }

	nb_global_grid_points = nbPoints_x;
	
	Points2D.Reallocate(num);
	Index.Reallocate(num);
        TetaPoints.Reallocate(num);
	NumPartie.Reallocate(num);
	num = 0;

	for (int i = 0; i < nbPoints_x; i++)
	  {
	    point(0) = xmin + i*step_x;
	    point(1) = ymin + i*step_y;
	    point(2) = zmin + i*step_z;
	    if (TranslatePoint(point, true, grid, pt2D, theta, n))
	      {
		Index(num) = i;
		Points2D(num) = pt2D;
		TetaPoints(num) = theta;
		NumPartie(num) = n;
		num++;
	      }
	  }
      }
    else if (this->type_output_file == CIRCLE)
      {
	// output on a 1-D circle 
	Real_wp step_teta(0), radius(0);
	if (nbPoints_x >= 1)
	  step_teta = (2.0*pi_wp)/nbPoints_x;
	
	center(0) = (xmin+xmax)/2.0;
	center(1) = (ymin+ymax)/2.0;
	radius = (xmax-xmin)/2.0; 
	
	int num = 0;
	for (int i = 0; i < nbPoints_x; i++)
	  {
	    point(0) = center(0) + radius*cos(step_teta*i);
	    point(1) = center(1) + radius*sin(step_teta*i);
	    if (TranslatePoint(point, false, grid, pt2D, theta, n))
	      num++;
	  }

	nb_global_grid_points = nbPoints_x;
	
	Points2D.Reallocate(num);
	Index.Reallocate(num);
        TetaPoints.Reallocate(num);
	NumPartie.Reallocate(num);
	num = 0;

	for (int i = 0; i < nbPoints_x; i++)
	  {
	    point(0) = center(0) + radius*cos(step_teta*i);
	    point(1) = center(1) + radius*sin(step_teta*i);
	    if (TranslatePoint(point, false, grid, pt2D, theta, n))
	      {
		Index(num) = i;
		Points2D(num) = pt2D;
		TetaPoints(num) = theta;
		NumPartie(num) = n;
		num++;
	      }
	  }
      }
    else if (this->type_output_file == CIRCLE_AXI)
      {
	Real_wp step_teta = 2.0*pi_wp/nbPoints_x;
	R3 normale, vec_u, vec_v; Real_wp radius_x = xmax, radius_y = ymax;
	normale(0) = xmin; normale(1) = ymin; normale(2) = zmin;
	GetVectorPlane(normale, vec_u, vec_v);
	//DISP(vec_u); DISP(vec_v);
	
	int num = 0;
	for (int i = 0; i < nbPoints_x; i++)
	  {
	    Real_wp teta = step_teta*i;
	    point = center; Add(radius_x*cos(teta), vec_u, point);
	    Add(radius_y*sin(teta), vec_v, point);
	    if (TranslatePoint(point, true, grid, pt2D, theta, n))
	      num++;
	  }
	
	nb_global_grid_points = nbPoints_x;
	
	Points2D.Reallocate(num);
	Index.Reallocate(num);
        TetaPoints.Reallocate(num);
	NumPartie.Reallocate(num);
	num = 0;

	for (int i = 0; i < nbPoints_x; i++)
	  {
	    Real_wp teta = step_teta*i;
	    point = center; Add(radius_x*cos(teta), vec_u, point);
	    Add(radius_y*sin(teta), vec_v, point);
	    if (TranslatePoint(point, true, grid, pt2D, theta, n))
	      {
		Index(num) = i;
		Points2D(num) = pt2D;
		TetaPoints(num) = theta;
		NumPartie(num) = n;
		num++;
	      }
	  }
      }
    else if (this->type_output_file == POINT)
      {
	// output on a single point
	point(0) = xmin; point(1) = ymin;
	nb_global_grid_points = 1;
	if (TranslatePoint(point, false, grid, pt2D, theta, n))
	  {
	    Points2D.Reallocate(1);
	    Index.Reallocate(1);
	    TetaPoints.Reallocate(1);
	    NumPartie.Reallocate(1);
	    
	    Index(0) = 0;
	    Points2D(0) = pt2D;
	    TetaPoints(0) = theta;
	    NumPartie(0) = n;
	  }
	else
	  {
	    Index.Clear();
	    Points2D.Clear();
	    TetaPoints.Clear();
	    NumPartie.Clear();
	  }
       }
    else if (type_output_file == POINT_AXI)
      {
	// output on single point
	point.Init(xmin, ymin, zmin);
	nb_global_grid_points = 1;
	if (TranslatePoint(point, true, grid, pt2D, theta, n))
	  {
	    Points2D.Reallocate(1);
	    Index.Reallocate(1);
	    TetaPoints.Reallocate(1);
	    NumPartie.Reallocate(1);
	    
	    Index(0) = 0;
	    Points2D(0) = pt2D;
	    TetaPoints(0) = theta;
	    NumPartie(0) = n;
	  }
	else
	  {
	    Index.Clear();
	    Points2D.Clear();
	    TetaPoints.Clear();
	    NumPartie.Clear();
	  }
      }        
  }


  //! translates a point to the physical domain and returns true if 
  //! the translated point is inside the grid
  bool GridInterpolationFull<Dimension2>::
  TranslatePoint(const R3& point, bool axi, const GridInterpolation<Dimension2>& grid,
		 R2& pt2D, Real_wp& theta, int& n)
  {    
    Real_wp r, phi;
    theta = Real_wp(0); n = 0;
    if (axi)
      {
	CartesianToPolar(point(0), point(1), r, theta);
	pt2D.Init(r, point(2));
      }
    else
      pt2D.Init(point(0), point(1));
    
    if (grid.GetNbCyclicSections() > 1)
      {
	// cyclic domain
	CartesianToPolar(pt2D(0), pt2D(1), r, phi);
	
	if (phi < 0)
	  phi += 2.0*pi_wp;
            
	int nb_sections = grid.GetNbCyclicSections();
	Real_wp nreal = phi*Real_wp(nb_sections)/(2.0*pi_wp);
	
	if (abs(round(nreal) - nreal) < 1e-10)
	  nreal = round(nreal);
	
	n = toInteger(floor(nreal));
	
	Real_wp phi_bis = phi - ( n * (2.0*pi_wp ) / Real_wp(nb_sections)) ;
	pt2D.Init(r * cos(phi_bis), r * sin(phi_bis));
      }
    else
      {
	// periodicity with respect to x and y ?
        int nx, ny, nz;
        grid.GetNbCartesianSections(nx, ny, nz);
        if (nx < 1)
          nx = 1;
        
        if (ny < 1)
          ny = 1;

        if ((nx != 1) || (ny != 1))
          {
            // Delta x and Delta y for the initial cell
            Real_wp dx = this->xmax0 - this->xmin0;
            Real_wp dy = this->ymax0 - this->ymin0;
            // position of the central cell
            int posx = (nx+1)/2;
            int posy = (ny+1)/2;
            // extremum values of the supercell
            Real_wp x_max_ = this->xmin0 + posx*dx;
            Real_wp y_max_ = this->ymin0 + posy*dy;
            Real_wp x_min_ = x_max_ - nx*dx;
            Real_wp y_min_ = y_max_ - ny*dy;

	    Real_wp x = pt2D(0);
	    Real_wp y = pt2D(1);
	    // we detect on which cell the point is located
	    Real_wp ix_ = (x - x_min_)/dx;
	    Real_wp iy_ = (y - y_min_)/dy;
	    int ix = toInteger(round(ix_));
	    int iy = toInteger(round(iy_));
	    if (abs(ix_ - ix) > 1e-12)
	      ix = toInteger(floor(ix_));
	    
	    if (abs(iy_ - iy) > 1e-12)
	      iy = toInteger(floor(iy_));

	    if ((ix < 0) || (ix >= nx))
	      return false;

	    if ((iy < 0) || (iy >= ny))
	      return false;
	    
	    // then we translate the point to the central cell
	    pt2D(0) = this->xmin0 + (x - x_min_ - dx*ix);
	    pt2D(1) = this->ymin0 + (y - y_min_ - dy*iy);
	    
	    n = ix*ny + iy;
	  }
      }
    
    // we test if the point is inside the bounding box of the grid
    if ((pt2D(0) >= grid.GetXmin()-R2::threshold) && (pt2D(0) <= grid.GetXmax()+R2::threshold) 
	&& (pt2D(1) >= grid.GetYmin()-R2::threshold) && (pt2D(1) <= grid.GetYmax()+R2::threshold))
      return true;
    
    return false;
  }
  
  
  //! computation of global points if it is a regular grid, line, circle
  /*!
    \param[in] grid object used to localize points on the mesh
    this methods computes global points in the case of regular grid,
    line divided in equal subdivisions, and appends these points to the object grid.
  */
  void GridInterpolationFull<Dimension2>::InitGrid(GridInterpolation<Dimension2>& grid)
  {
    this->var_interp = &grid;
    if (this->grid_interval_to_be_computed)
      {
	xmin = xmin0; ymin = ymin0; zmin = zmin0;
	xmax = xmax0; ymax = ymax0; zmax = zmax0;
	
	if ( (this->type_output_file == THREE_PLANES_AXI)
             || (this->type_output_file == VOLUME_AXI) )
	  {
	    xmin = -xmax;
	    zmin = ymin; zmax = ymax;
	    ymin = xmin; ymax = xmax;
	  }
        
        xmin += 5.0*(xmax-xmin)*epsilon_machine;
        xmax -= 5.0*(xmax-xmin)*epsilon_machine;
        ymin += 5.0*(ymax-ymin)*epsilon_machine;
        ymax -= 5.0*(ymax-ymin)*epsilon_machine;
        zmin += 5.0*(zmax-zmin)*epsilon_machine;
        zmax -= 5.0*(zmax-zmin)*epsilon_machine;
      }

    VectR2 Points2D; VectReal_wp TetaPoints; IVect NumPartie;
    GenerateGridPoints(grid, IndexPoints, Points2D, TetaPoints, NumPartie);
    
    grid.Append(Points2D, TetaPoints, NumPartie, list_points);
  }
  
  
  //! writes the header of display grid in binary format
  void GridInterpolationFull<Dimension2>::Write(ostream& FileStream, bool double_prec) const
  {
    if (double_prec)
      {
	double a;
	WriteBinary(FileStream, a);
      }
    else
      {
	float a;
	WriteBinary(FileStream, a);
      }
  }
  
  
  //! writes the header of display grid in binary format
  template<class real>
  void GridInterpolationFull<Dimension2>::WriteBinary(ostream& FileStream, real& a) const
  {
    real tmp; int nb;
    switch (type_output_file)
      {
      case PLANE :
	// for 2-D regular grid, xmin, xmax, nbx, ymin, ymax, nby
	tmp = toDouble(xmin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(xmax*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&nbPoints_x)),sizeof(int));
	tmp = toDouble(ymin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp),sizeof(real));
	tmp = toDouble(ymax*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp),sizeof(real));
	FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&nbPoints_y)),sizeof(int));	
	break;
      case PLANE_AXI :
	// plane OAB with O = (xmin, ymin, zmin)
	// A = (xmax, ymax, zmax), B = center
        tmp = toDouble(xmin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(xmax*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(ymin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(ymax*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(zmin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(zmax*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));

	tmp = toDouble(center(0)*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(center(1)*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(center(2)*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	nb = nbPoints_x; FileStream.write(reinterpret_cast<char*>(&nb), sizeof(int));
	nb = nbPoints_y; FileStream.write(reinterpret_cast<char*>(&nb), sizeof(int));
	break;
      case VOLUME_AXI :
	tmp = toDouble(xmin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(xmax*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	nb = nbPoints_x; FileStream.write(reinterpret_cast<char*>(&nb), sizeof(int));
	tmp = toDouble(ymin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(ymax*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	nb = nbPoints_y; FileStream.write(reinterpret_cast<char*>(&nb), sizeof(int));
	tmp = toDouble(zmin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(zmax*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	nb = nbPoints_z; FileStream.write(reinterpret_cast<char*>(&nb), sizeof(int));
	break;
      case LINE :
	// for a line : xmin, ymin, xmax, ymax and nbx
	tmp = toDouble(xmin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(xmax*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(ymin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp),sizeof(real));
	tmp = toDouble(ymax*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp),sizeof(real));
	FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&nbPoints_x)),sizeof(int));
	break;
      case LINE_AXI :
	// line joining A=(xmin,ymin,zmin) and B = (xmax, ymax, zmax)
	tmp = toDouble(xmin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(xmax*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(ymin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(ymax*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(zmin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(zmax*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	nb = nbPoints_x; FileStream.write(reinterpret_cast<char*>(&nb), sizeof(int));
	break;
      case POINT :
	tmp = toDouble(xmin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(ymin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	break;
      case POINT_AXI :
	tmp = toDouble(xmin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(ymin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(zmin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	break;
      case CIRCLE :
	tmp = toDouble(0.5*(xmin+xmax)*z0_adim);
        FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(0.5*(ymin+ymax)*z0_adim);
        FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(0.5*(xmax-xmin)*z0_adim);
        FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&nbPoints_x)),sizeof(int));
	break;
      case THREE_PLANES_AXI :
	// for 2-D regular grid, xmin, xmax, nbx, ymin, ymax, nby
	tmp = toDouble(xmin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(xmax*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&nbPoints_x)),sizeof(int));
	tmp = toDouble(ymin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp),sizeof(real));
	tmp = toDouble(ymax*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp),sizeof(real));
	FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&nbPoints_y)),sizeof(int));	
	tmp = toDouble(zmin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp),sizeof(real));
	tmp = toDouble(zmax*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp),sizeof(real));
	FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&nbPoints_z)),sizeof(int));	
	tmp = toDouble(center(0)*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp),sizeof(real));
	tmp = toDouble(center(1)*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp),sizeof(real));
	tmp = toDouble(center(2)*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp),sizeof(real));
	break;
      case CIRCLE_AXI :
	tmp = toDouble(center(0)*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(center(1)*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(center(2)*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(xmin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(ymin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(zmin*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(xmax*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(ymax*z0_adim); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&nbPoints_x)),sizeof(int));
	break;
      case POINTS_FILE :
	{
	  nb = GlobalPoints2D.GetM();
	  FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&nb)),sizeof(int));
	  Vector<real> Pts2D(2*nb);
	  for (int i = 0; i < nb; i++)
	    {
	      Pts2D(2*i) = toDouble(GlobalPoints2D(i)(0)*z0_adim); 
	      Pts2D(2*i+1) = toDouble(GlobalPoints2D(i)(1)*z0_adim);
	    }
	  FileStream.write(reinterpret_cast<char*>(const_cast<real*>(Pts2D.GetData())),
			   2*nb*sizeof(real));
	}
	break;
      case POINTS_FILE_AXI :
	{
	  nb = GlobalPoints3D.GetM();
	  FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&nb)),sizeof(int));
	  Vector<real> Pts3D(3*nb);
	  for (int i = 0; i < nb; i++)
	    {
	      Pts3D(3*i) = toDouble(GlobalPoints3D(i)(0)*z0_adim); 
	      Pts3D(3*i+1) = toDouble(GlobalPoints3D(i)(1)*z0_adim);
	      Pts3D(3*i+2) = toDouble(GlobalPoints3D(i)(2)*z0_adim);
	    }
	  FileStream.write(reinterpret_cast<char*>(const_cast<real*>(Pts3D.GetData())),
			   3*nb*sizeof(real));
	}
	break;
      }
  }
  
  
  //! writes headers of display grid in text format
  void GridInterpolationFull<Dimension2>
  ::WriteText(ostream& FileStream, int type_precision) const
  {
    switch (type_precision)
      {
      case OutputTypeEnum::SINGLE_PRECISION : FileStream.precision(7); break;
      case OutputTypeEnum::DOUBLE_PRECISION : FileStream.precision(15); break;
      case OutputTypeEnum::QUADRUPLE_PRECISION : FileStream.precision(32); break;
      case OutputTypeEnum::MULTIPLE_PRECISION :
        {
          int n = toInteger(-log10(epsilon_machine));
          FileStream.precision(n);
        }
        break;
      }
    
    switch (this->type_output_file)
      {
      case PLANE :
	FileStream << xmin*z0_adim << endl << xmax*z0_adim << endl << nbPoints_x << endl;
	FileStream << ymin*z0_adim << endl << ymax*z0_adim << endl << nbPoints_y << endl;        
	break;
      case PLANE_AXI :
	// plane OAB with O = (xmin, ymin, zmin)
	// A = (xmax, ymax, zmax), B = center
	FileStream << xmin*z0_adim << endl << xmax*z0_adim << endl;
	FileStream << ymin*z0_adim << endl << ymax*z0_adim << endl;
	FileStream << zmin*z0_adim << endl << zmax*z0_adim << endl;
	FileStream << center(0)*z0_adim << endl << center(1)*z0_adim << endl
                   << center(2)*z0_adim << endl;
	FileStream << nbPoints_x << endl << nbPoints_y << endl;
	break;
      case LINE :
	FileStream << xmin*z0_adim << endl << xmax*z0_adim << endl;
	FileStream << ymin*z0_adim << endl << ymax*z0_adim << endl << nbPoints_x << endl;
	break;
      case LINE_AXI :
	// line joining A=(xmin,ymin,zmin) and B = (xmax, ymax, zmax)
	FileStream << xmin*z0_adim << endl << xmax*z0_adim << endl;
	FileStream << ymin*z0_adim << endl << ymax*z0_adim << endl;
	FileStream << zmin*z0_adim << endl << zmax*z0_adim << endl;
	FileStream << nbPoints_x << endl;
	break;
      case POINT :
	FileStream << xmin*z0_adim << endl << ymin*z0_adim << endl;
	break;
      case POINT_AXI :
	FileStream << xmin*z0_adim << endl;
	FileStream << ymin*z0_adim << endl;
	FileStream << zmin*z0_adim << endl;
	break;
      case THREE_PLANES_AXI :
	FileStream << xmin*z0_adim << endl << xmax*z0_adim << endl << nbPoints_x << endl;
	FileStream << ymin*z0_adim << endl << ymax*z0_adim << endl << nbPoints_y << endl;
	FileStream << zmin*z0_adim << endl << zmax*z0_adim << endl << nbPoints_z << endl;
	FileStream << center(0)*z0_adim << endl << center(1)*z0_adim 
                   << endl << center(2)*z0_adim << endl;
	break;
      case CIRCLE :
	FileStream << 0.5*(xmin+xmax)*z0_adim << endl << 0.5*(ymin+ymax)*z0_adim << endl
		   << 0.5*(xmax-xmin)*z0_adim << endl << nbPoints_x << endl;
	break;
      case CIRCLE_AXI :
	FileStream << center(0)*z0_adim << endl << center(1)*z0_adim << endl
		   << center(2)*z0_adim << endl << xmin*z0_adim << endl << ymin*z0_adim << endl
		   << zmin*z0_adim << endl << xmax*z0_adim << endl 
                   << ymax*z0_adim << endl << nbPoints_x << endl;
	break;
      case POINTS_FILE :
	FileStream << GlobalPoints2D.GetM() << endl;
	for (int i = 0; i < GlobalPoints2D.GetM(); i++)
	  {
            R2 pt = GlobalPoints2D(i); pt *= z0_adim;
	    PrintNoBrackets(FileStream, pt);
	    FileStream << '\n';
	  }	
	break;
      case VOLUME_AXI :
	FileStream << xmin*z0_adim << endl << xmax*z0_adim << endl << nbPoints_x << endl;
	FileStream << ymin*z0_adim << endl << ymax*z0_adim << endl << nbPoints_y << endl;
	FileStream << zmin*z0_adim << endl << zmax*z0_adim << endl << nbPoints_z << endl;
	break;
      case POINTS_FILE_AXI :
	FileStream << GlobalPoints3D.GetM() << endl;
	for (int i = 0; i < GlobalPoints3D.GetM(); i++)
	  {
            R3 pt = GlobalPoints3D(i); pt *= z0_adim;
	    PrintNoBrackets(FileStream, pt);
	    FileStream << '\n';
	  }
	break;
      }
    
  }
  
  
  //! writes headers of the grid in Vtk format
  void GridInterpolationFull<Dimension2>
  ::WriteVtk(ostream& FileStream, bool double_prec, bool ascii) const
  {    
    // main header of vtk file
    FileStream << "# vtk DataFile Version 3.0\n";
    FileStream << "Montjoie output\n";
    
    if (!ascii)
      FileStream << "BINARY\n";
    else
      FileStream << "ASCII\n";
    
    if (double_prec)
      FileStream.precision(15);
    else
      FileStream.precision(7);
    
    switch (this->type_output_file)
      {
      case PLANE :
        {
          Real_wp dx(1), dy(1), dz(1);
          if (nbPoints_x > 1)
            dx = (xmax-xmin)/(nbPoints_x-1);
          
          if (nbPoints_y > 1)
            dy = (ymax-ymin)/(nbPoints_y-1);
          
          FileStream << "DATASET STRUCTURED_POINTS\n";
          FileStream << "DIMENSIONS " << nbPoints_x << " " << nbPoints_y << " 1\n";
          FileStream << "ORIGIN " << xmin*z0_adim << " " << ymin*z0_adim << " " << zmin <<'\n';
          FileStream << "SPACING " << dx*z0_adim << " " << dy*z0_adim 
                     << " " << dz*z0_adim <<'\n';         
	  FileStream << "POINT_DATA " << nbPoints_x*nbPoints_y << endl; 
        }
        break;
      case PLANE_AXI :
        {
	  R3 ptO(xmin, ymin, zmin), OA(xmax, ymax, zmax), OB = center, pt;
	  OA -= ptO; OB -= ptO;
	  
	  Real_wp lambda, mu, dx = 1.0, dy = 1.0;
	  if (nbPoints_x > 1)
	    dx = 1.0/(nbPoints_x-1);
	  
	  if (nbPoints_y > 1)
	    dy = 1.0/(nbPoints_y-1);
	  
	  FileStream << "DATASET STRUCTURED_GRID\n";
	  FileStream << "DIMENSIONS " << nbPoints_x << " " << nbPoints_y << " 1\n";
	  if (double_prec)
	    FileStream << "POINTS " << nbPoints_x*nbPoints_y << " double\n";
	  else
	    FileStream << "POINTS " << nbPoints_x*nbPoints_y << " float\n";
	  
	  VectReal_wp ptV(3);
	  for (int j = 0; j < nbPoints_y; j++)
	    for (int i = 0; i < nbPoints_x; i++)
	      {
		// computing position of point (i, j)
		lambda = i*dx;
		mu = j*dy;
		pt = ptO + lambda*OA + mu*OB;
		if (ascii)
		  FileStream << pt(0)*z0_adim << ' ' << pt(1)*z0_adim << ' ' 
                             << pt(2)*z0_adim << '\n';
		else
		  {
		    ptV(0) = pt(0)*z0_adim; ptV(1) = pt(1)*z0_adim; ptV(2) = pt(2)*z0_adim;
		    WriteBinaryDoubleOrFloat(ptV, FileStream, double_prec, false, true);
		  }
	      }
	  
	  if (!ascii)
	    FileStream << '\n';
	  
	  FileStream << "POINT_DATA " << nbPoints_x*nbPoints_y << endl; 
        }
        break;
      case VOLUME_AXI :
        {
          Real_wp dx(1), dy(1), dz(1);
          if (nbPoints_x > 1)
            dx = (xmax-xmin)/(nbPoints_x-1);
          
          if (nbPoints_y > 1)
            dy = (ymax-ymin)/(nbPoints_y-1);

          if (nbPoints_z > 1)
            dz = (zmax-zmin)/(nbPoints_z-1);
          
          FileStream << "DATASET STRUCTURED_POINTS\n";
          FileStream << "DIMENSIONS " << nbPoints_x << " " << nbPoints_y 
                     << " " << nbPoints_z <<'\n';
          FileStream << "ORIGIN " << xmin*z0_adim << " " << ymin*z0_adim 
                     << " " << zmin*z0_adim <<'\n';
          FileStream << "SPACING " << dx*z0_adim << " " << dy*z0_adim << " " << dz*z0_adim <<'\n';
	  FileStream << "POINT_DATA " << nbPoints_x*nbPoints_y*nbPoints_z << endl; 
        }
        break;
      case LINE :
	{
	  R2 ptA(xmin, ymin), ptB(xmax, ymax), pt;
	  
	  Real_wp lambda, dx = 1.0;
	  if (nbPoints_x > 1)
	    dx = 1.0/(nbPoints_x-1);
	  
	  FileStream << "DATASET STRUCTURED_GRID\n";
	  FileStream << "DIMENSIONS " << nbPoints_x << " 1 1\n";
	  if (double_prec)
	    FileStream << "POINTS " << nbPoints_x << " double\n";
	  else
	    FileStream << "POINTS " << nbPoints_x << " float\n";
	  
	  VectReal_wp ptV(3);
	  for (int i = 0; i < nbPoints_x; i++)
	    {
	      // computing position of point i
	      lambda = i*dx;
	      pt = (1.0-lambda)*ptA + lambda*ptB;
	      if (ascii)
		FileStream << pt(0)*z0_adim << ' ' << pt(1)*z0_adim << " 0\n";
	      else
		{
		  ptV(0) = pt(0)*z0_adim; ptV(1) = pt(1)*z0_adim; ptV(2) = 0.0;
		  WriteBinaryDoubleOrFloat(ptV, FileStream, double_prec, false, true);
		}
	    }
	  
	  if (!ascii)
	    FileStream << '\n';
	  
	  FileStream << "POINT_DATA " << nbPoints_x << endl; 
	}
	break;
      case LINE_AXI :
	{
	  R3 ptA(xmin, ymin, zmin), ptB(xmax, ymax, zmax), pt;
	  
	  Real_wp lambda, dx = 1.0;
	  if (nbPoints_x > 1)
	    dx = 1.0/(nbPoints_x-1);
	  
	  FileStream << "DATASET STRUCTURED_GRID\n";
	  FileStream << "DIMENSIONS " << nbPoints_x << " 1 1\n";
	  if (double_prec)
	    FileStream << "POINTS " << nbPoints_x << " double\n";
	  else
	    FileStream << "POINTS " << nbPoints_x << " float\n";
	  
	  VectReal_wp ptV(3);
	  for (int i = 0; i < nbPoints_x; i++)
	    {
	      // computing position of point i
	      lambda = i*dx;
	      pt = (1.0-lambda)*ptA + lambda*ptB;
	      if (ascii)
		FileStream << pt(0)*z0_adim << ' ' << pt(1)*z0_adim 
                           << ' ' << pt(2)*z0_adim << '\n';
	      else
		{
		  ptV(0) = pt(0)*z0_adim; ptV(1) = pt(1)*z0_adim; ptV(2) = pt(2)*z0_adim;
		  WriteBinaryDoubleOrFloat(ptV, FileStream, double_prec, false, true);
		}
	    }
	  
	  if (!ascii)
	    FileStream << '\n';
	  
	  FileStream << "POINT_DATA " << nbPoints_x << endl; 
	}
	break;
      case POINT :
	{
	  FileStream << "DATASET STRUCTURED_GRID\n";
	  FileStream << "DIMENSIONS 1 1 1\n";
	  if (double_prec)
	    FileStream << "POINTS 1 double\n";
	  else
	    FileStream << "POINTS 1 float\n";
	  
	  VectReal_wp ptV(3);
	  if (ascii)
	    FileStream << xmin*z0_adim << ' ' << ymin*z0_adim << " 0\n";
	  else
	    {
	      ptV(0) = xmin*z0_adim; ptV(1) = ymin*z0_adim; ptV(2) = 0.0;
	      WriteBinaryDoubleOrFloat(ptV, FileStream, double_prec, false, true);
	    }
	  
	  if (!ascii)
	    FileStream << '\n';
	  
	  FileStream << "POINT_DATA 1" << endl; 
	}
	break;
      case POINT_AXI :
	{
	  FileStream << "DATASET STRUCTURED_GRID\n";
	  FileStream << "DIMENSIONS 1 1 1\n";
	  if (double_prec)
	    FileStream << "POINTS 1 double\n";
	  else
	    FileStream << "POINTS 1 float\n";
	  
	  VectReal_wp ptV(3);
	  if (ascii)
	    FileStream << xmin*z0_adim << ' ' << ymin*z0_adim << ' ' << zmin*z0_adim << '\n';
	  else
	    {
	      ptV(0) = xmin*z0_adim; ptV(1) = ymin*z0_adim; ptV(2) = zmin*z0_adim;
	      WriteBinaryDoubleOrFloat(ptV, FileStream, double_prec, false, true);
	    }
	  
	  if (!ascii)
	    FileStream << '\n';
	  
	  FileStream << "POINT_DATA 1" << endl; 
	}
	break;
      case THREE_PLANES_AXI :
	{
	  FileStream << "DATASET UNSTRUCTURED_GRID\n";
	  int nbPoints = nbPoints_x*(nbPoints_y+nbPoints_z) + nbPoints_y*nbPoints_z;
	  FileStream << "POINTS " << nbPoints;
	  
	  if (double_prec)
	    FileStream << " double\n";
	  else
	    FileStream << " float\n";
	  
	  Real_wp x, dx(0), y, dy(0), z, dz(0);
	  if (nbPoints_x > 1)
	    dx = (xmax - xmin)/(nbPoints_x-1);
	  if (nbPoints_y > 1)
	    dy = (ymax - ymin)/(nbPoints_y-1);
	  if (nbPoints_z > 1)
	    dz = (zmax - zmin)/(nbPoints_z-1);
	  
	  VectReal_wp ptV(3);
	  for (int j = 0; j < nbPoints_y; j++)
	    for (int k = 0; k < nbPoints_z; k++)
	      {
		x = center(0);
		y = ymin + j*dy;
		z = zmin + k*dz;
		if (ascii)
		  FileStream << x*z0_adim << ' ' << y*z0_adim << ' ' << z*z0_adim << '\n';
		else
		  {
		    ptV(0) = x*z0_adim; ptV(1) = y*z0_adim; ptV(2) = z*z0_adim;
		    WriteBinaryDoubleOrFloat(ptV, FileStream, double_prec, false, true);
		  }				
	      }

	  for (int i = 0; i < nbPoints_x; i++)
	    for (int k = 0; k < nbPoints_z; k++)
	      {
		x = xmin + i*dx;
		y = center(1);
		z = zmin + k*dz;
		if (ascii)
		  FileStream << x*z0_adim << ' ' << y*z0_adim << ' ' << z*z0_adim << '\n';
		else
		  {
		    ptV(0) = x*z0_adim; ptV(1) = y*z0_adim; ptV(2) = z*z0_adim;
		    WriteBinaryDoubleOrFloat(ptV, FileStream, double_prec, false, true);
		  }				
	      }

	  for (int i = 0; i < nbPoints_x; i++)
	    for (int j = 0; j < nbPoints_y; j++)
	      {
		x = xmin + i*dx;
		y = ymin + j*dy;
		z = center(2);
		if (ascii)
		  FileStream << x*z0_adim << ' ' << y*z0_adim << ' ' << z*z0_adim << '\n';
		else
		  {
		    ptV(0) = x*z0_adim; ptV(1) = y*z0_adim; ptV(2) = z*z0_adim;
		    WriteBinaryDoubleOrFloat(ptV, FileStream, double_prec, false, true);
		  }				
	      }
	  
	  if (!ascii)
	    FileStream << '\n';
	  
	  int nb_elt = (nbPoints_x-1)*(nbPoints_y+nbPoints_z-2)
	    + (nbPoints_y-1)*(nbPoints_z-1);
	  
	  FileStream << "CELLS " << nb_elt << " " << 5*nb_elt << '\n';
	  int n1 = 0, n2 = nbPoints_z, n3 = nbPoints_z+1, n4 = 1; 
	  IVect num(5);
	  for (int j = 0; j < nbPoints_y-1; j++)
	    for (int k = 0; k < nbPoints_z-1; k++)
	      {
		if (ascii)
		  FileStream << "4 " << n1 << ' ' << n2 << ' ' << n3 << ' ' << n4 << '\n';
		else
		  {
		    num(0) = swapEndian(4);
		    num(1) = swapEndian(n1); num(2) = swapEndian(n2);
		    num(3) = swapEndian(n3); num(4) = swapEndian(n4);
		    FileStream.write(reinterpret_cast<char*>(num.GetData()), 5*sizeof(int));
		  }
		
		n1++; n2++; n3++; n4++;
		if (k == nbPoints_z-2)
		  {
		    n1++; n2++; n3++; n4++;
		  }
	      }
	  
	  n1 = nbPoints_y*nbPoints_z;
	  n2 = nbPoints_y*nbPoints_z + nbPoints_z;
	  n3 = nbPoints_y*nbPoints_z + nbPoints_z+1;
	  n4 = nbPoints_y*nbPoints_z + 1; 
	  for (int j = 0; j < nbPoints_x-1; j++)
	    for (int k = 0; k < nbPoints_z-1; k++)
	      {
		if (ascii)
		  FileStream << "4 " << n1 << ' ' << n2 << ' ' << n3 << ' ' << n4 << '\n';
		else
		  {
		    num(0) = swapEndian(4);
		    num(1) = swapEndian(n1); num(2) = swapEndian(n2);
		    num(3) = swapEndian(n3); num(4) = swapEndian(n4);
		    FileStream.write(reinterpret_cast<char*>(num.GetData()), 5*sizeof(int));
		  }

		n1++; n2++; n3++; n4++;
		if (k == nbPoints_z-2)
		  {
		    n1++; n2++; n3++; n4++;
		  }
	      }

	  n1 = (nbPoints_x+nbPoints_y)*nbPoints_z;
	  n2 = (nbPoints_x+nbPoints_y)*nbPoints_z + nbPoints_y;
	  n3 = (nbPoints_x+nbPoints_y)*nbPoints_z + nbPoints_y+1;
	  n4 = (nbPoints_x+nbPoints_y)*nbPoints_z + 1; 
	  for (int j = 0; j < nbPoints_x-1; j++)
	    for (int k = 0; k < nbPoints_y-1; k++)
	      {
		if (ascii)
		  FileStream << "4 " << n1 << ' ' << n2 << ' ' << n3 << ' ' << n4 << '\n';
		else
		  {
		    num(0) = swapEndian(4);
		    num(1) = swapEndian(n1); num(2) = swapEndian(n2);
		    num(3) = swapEndian(n3); num(4) = swapEndian(n4);
		    FileStream.write(reinterpret_cast<char*>(num.GetData()), 5*sizeof(int));
		  }

		n1++; n2++; n3++; n4++;
		if (k == nbPoints_y-2)
		  {
		    n1++; n2++; n3++; n4++;
		  }
	      }
	  
	  if (!ascii)
	    FileStream << '\n';
	  
	  FileStream << "CELL_TYPES " << nb_elt << '\n';
	  for (int i = 0; i < nb_elt; i++)
	    {
	      if (ascii)
		FileStream << "8\n";
	      else
		{
		  num(0) = swapEndian(int(9));
		  FileStream.write(reinterpret_cast<char*>(num.GetData()), sizeof(int));
		}
	    }

	  if (!ascii)
	    FileStream << '\n';
	  
	  FileStream << "POINT_DATA " << nbPoints << endl; 
	}
	break;
      case CIRCLE :
	{
	  R2 ptA(0.5*(xmin+xmax), 0.5*(ymin+ymax));
	  Real_wp radius = 0.5*(xmax-xmin), x, y;
	  
	  Real_wp teta, dteta = 1.0;
	  if (nbPoints_x > 1)
	    dteta = 2.0*pi_wp/nbPoints_x;
	  
	  FileStream << "DATASET STRUCTURED_GRID\n";
	  FileStream << "DIMENSIONS " << nbPoints_x << " 1 1\n";
	  if (double_prec)
	    FileStream << "POINTS " << nbPoints_x << " double\n";
	  else
	    FileStream << "POINTS " << nbPoints_x << " float\n";
	  
	  VectReal_wp ptV(3);
	  for (int i = 0; i < nbPoints_x; i++)
	    {
	      // computing position of point i
	      teta = i*dteta;
	      x = ptA(0) + radius*cos(teta);
	      y = ptA(1) + radius*sin(teta);
	      if (ascii)
		FileStream << x*z0_adim << ' ' << y*z0_adim << " 0\n";
	      else
		{
		  ptV(0) = x*z0_adim; ptV(1) = y*z0_adim; ptV(2) = 0.0;
		  WriteBinaryDoubleOrFloat(ptV, FileStream, double_prec, false, true);
		}
	    }
	  
	  if (!ascii)
	    FileStream << '\n';
	  
	  FileStream << "POINT_DATA " << nbPoints_x << endl; 
	}
	break;
      case CIRCLE_AXI :
	{
	  R3 normale, pt3D, vec_u, vec_v;
	  Real_wp radius_x = xmax, radius_y = ymax;
	  normale(0) = xmin; normale(1) = ymin; normale(2) = zmin;
	  GetVectorPlane(normale, vec_u, vec_v);
	  
	  Real_wp teta, dteta = 1.0;
	  if (nbPoints_x > 1)
	    dteta = 2.0*pi_wp/nbPoints_x;
	  
	  FileStream << "DATASET STRUCTURED_GRID\n";
	  FileStream << "DIMENSIONS " << nbPoints_x << " 1 1\n";
	  if (double_prec)
	    FileStream << "POINTS " << nbPoints_x << " double\n";
	  else
	    FileStream << "POINTS " << nbPoints_x << " float\n";
	  
	  VectReal_wp ptV(3);
	  for (int i = 0; i < nbPoints_x; i++)
	    {
	      // computing position of point i
	      teta = i*dteta;
	      pt3D = center;
	      Add(radius_x*cos(teta), vec_u, pt3D);
	      Add(radius_y*sin(teta), vec_v, pt3D);
	      if (ascii)
		FileStream << pt3D(0)*z0_adim << ' ' << pt3D(1)*z0_adim 
                           << ' ' << pt3D(2)*z0_adim << '\n';
	      else
		{
		  ptV(0) = pt3D(0)*z0_adim; ptV(1) = pt3D(1)*z0_adim; ptV(2) = pt3D(2)*z0_adim;
		  WriteBinaryDoubleOrFloat(ptV, FileStream, double_prec, false, true);
		}
	    }
	  
	  if (!ascii)
	    FileStream << '\n';
	  
	  FileStream << "POINT_DATA " << nbPoints_x << endl; 
	}
	break;
      case POINTS_FILE :
	{
	  FileStream << "DATASET STRUCTURED_GRID\n";
	  FileStream << "DIMENSIONS " << nbPoints_x << " 1 1\n";
	  if (double_prec)
	    FileStream << "POINTS " << nbPoints_x << " double\n";
	  else
	    FileStream << "POINTS " << nbPoints_x << " float\n";
	  
	  VectReal_wp ptV(3); Real_wp x, y;
	  for (int i = 0; i < nbPoints_x; i++)
	    {
	      x = GlobalPoints2D(i)(0)*z0_adim;
	      y = GlobalPoints2D(i)(1)*z0_adim;
	      if (ascii)
		FileStream << x << ' ' << y << " 0\n";
	      else
		{
		  ptV(0) = x; ptV(1) = y; ptV(2) = 0.0;
		  WriteBinaryDoubleOrFloat(ptV, FileStream, double_prec, false, true);
		}
	    }
	  
	  if (!ascii)
	    FileStream << '\n';
	  
	  FileStream << "POINT_DATA " << nbPoints_x << endl; 	  
	}
	break;
      case POINTS_FILE_AXI :
        {
	  FileStream << "DATASET STRUCTURED_GRID\n";
	  FileStream << "DIMENSIONS " << nbPoints_x << " 1 1\n";
	  if (double_prec)
	    FileStream << "POINTS " << nbPoints_x << " double\n";
	  else
	    FileStream << "POINTS " << nbPoints_x << " float\n";
	  
	  VectReal_wp ptV(3); Real_wp x, y, z;
	  for (int i = 0; i < nbPoints_x; i++)
	    {
	      x = GlobalPoints3D(i)(0)*z0_adim;
	      y = GlobalPoints3D(i)(1)*z0_adim;
	      z = GlobalPoints3D(i)(2)*z0_adim;
	      if (ascii)
		FileStream << x << ' ' << y << ' ' << z << '\n';
	      else
		{
		  ptV(0) = x; ptV(1) = y; ptV(2) = z;
		  WriteBinaryDoubleOrFloat(ptV, FileStream, double_prec, false, true);
		}
	    }
	  
	  if (!ascii)
	    FileStream << '\n';
	  
	  FileStream << "POINT_DATA " << nbPoints_x << endl; 	  
	}
	break;
      }
    
  }
  
  
  //! reads the display grid from a binary file
  void GridInterpolationFull<Dimension2>::Read(istream& FileStream, bool double_prec)
  {
    if (double_prec)
      {
	double a;
	ReadBinary(FileStream, a);
      }
    else
      {
	float a;
	ReadBinary(FileStream, a);
      }
  }
  

  //! reads the display grid from a binary file
  template<class real>
  void GridInterpolationFull<Dimension2>::ReadBinary(istream& FileStream, real& f)
  {
    real tmp; int nb;
    nbPoints_x = 1; nbPoints_y = 1; nbPoints_z = 1;
    switch (type_output_file)
      {
      case PLANE :
	// for 2-D regular grid, xmin, xmax, nbx, ymin, ymax, nby
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); xmin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); xmax = tmp;
	FileStream.read(reinterpret_cast<char*>(&nb),sizeof(int)); nbPoints_x = nb;
	FileStream.read(reinterpret_cast<char*>(&tmp),sizeof(real)); ymin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp),sizeof(real)); ymax = tmp;
	FileStream.read(reinterpret_cast<char*>(&nb),sizeof(int)); nbPoints_y = nb;
	break;
      case PLANE_AXI :
	// plane OAB with O = (xmin, ymin, zmin)
	// A = (xmax, ymax, zmax), B = center
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); xmin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); xmax = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); ymin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); ymax = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); zmin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); zmax = tmp;

	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); center(0) = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); center(1) = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); center(2) = tmp;
	FileStream.read(reinterpret_cast<char*>(&nb), sizeof(int)); nbPoints_x = nb;
	FileStream.read(reinterpret_cast<char*>(&nb), sizeof(int)); nbPoints_y = nb;
	break;
      case VOLUME_AXI :
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); xmin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); xmax = tmp;
	FileStream.read(reinterpret_cast<char*>(&nb), sizeof(int)); nbPoints_x = nb;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); ymin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); ymax = tmp;
	FileStream.read(reinterpret_cast<char*>(&nb), sizeof(int)); nbPoints_y = nb;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); zmin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); zmax = tmp;
	FileStream.read(reinterpret_cast<char*>(&nb), sizeof(int)); nbPoints_z = nb;
	break;
      case LINE :
	// for a line : xmin, ymin, xmax, ymax and nbx
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); xmin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); xmax = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp),sizeof(real)); ymin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp),sizeof(real)); ymax = tmp;
	FileStream.read(reinterpret_cast<char*>(&nb),sizeof(int)); nbPoints_x = nb;
	break;
      case LINE_AXI :
	// line joining A=(xmin,ymin,zmin) and B = (xmax, ymax, zmax)
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); xmin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); xmax = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); ymin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); ymax = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); zmin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); zmax = tmp;
	FileStream.read(reinterpret_cast<char*>(&nb), sizeof(int)); nbPoints_x = nb;
	break;
      case POINT :
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); xmin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); ymin = tmp;
	break;
      case POINT_AXI :
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); xmin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); ymin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); zmin = tmp;
	break;
      case THREE_PLANES_AXI :
	// for 2-D regular grid, xmin, xmax, nbx, ymin, ymax, nby
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); xmin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); xmax = tmp;
	FileStream.read(reinterpret_cast<char*>(&nb),sizeof(int)); nbPoints_x = nb;
	FileStream.read(reinterpret_cast<char*>(&tmp),sizeof(real)); ymin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp),sizeof(real)); ymax = tmp;
	FileStream.read(reinterpret_cast<char*>(&nb),sizeof(int)); nbPoints_y = nb;
	FileStream.read(reinterpret_cast<char*>(&tmp),sizeof(real)); zmin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp),sizeof(real)); zmax = tmp;
	FileStream.read(reinterpret_cast<char*>(&nb),sizeof(int)); nbPoints_z = nb;
	FileStream.read(reinterpret_cast<char*>(&tmp),sizeof(real)); center(0) = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp),sizeof(real)); center(1) = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp),sizeof(real)); center(2) = tmp;
	break;
      case CIRCLE :
	{
	  Real_wp x0, y0, radius;
	  FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); x0 = tmp;
	  FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); y0 = tmp;
	  FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); radius = tmp;
	  FileStream.read(reinterpret_cast<char*>(&nb),sizeof(int)); nbPoints_x = nb;
	  xmin = x0 - radius; ymin = y0 - radius; 
	  xmax = x0 + radius; ymax = y0 + radius; 
	}
	break;
      case CIRCLE_AXI :
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); center(0) = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); center(1) = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); center(2) = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); xmin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); ymin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); zmin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); xmax = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); ymax = tmp;
	FileStream.read(reinterpret_cast<char*>(&nb),sizeof(int)); nbPoints_x = nb;
	break;
      case POINTS_FILE :
	{
	  FileStream.read(reinterpret_cast<char*>(&nb),sizeof(int));
	  Vector<real> Pts2D(2*nb);
	  FileStream.read(reinterpret_cast<char*>(Pts2D.GetData()),2*nb*sizeof(real));
	  nbPoints_x = nb;
	  GlobalPoints2D.Reallocate(nb);
	  for (int i = 0; i < nb; i++)
	    GlobalPoints2D(i).Init(Pts2D(2*i), Pts2D(2*i+1));
	  
	}
	break;
      case POINTS_FILE_AXI :
	{
	  FileStream.read(reinterpret_cast<char*>(&nb),sizeof(int));
	  Vector<real> Pts3D(3*nb);
	  FileStream.read(reinterpret_cast<char*>(Pts3D.GetData()), 3*nb*sizeof(real));
	  nbPoints_x = nb;
	  GlobalPoints3D.Reallocate(nb);
	  for (int i = 0; i < nb; i++)
	    GlobalPoints3D(i).Init(Pts3D(3*i), Pts3D(3*i+1), Pts3D(3*i+2));
	}
	break;
      }
  }
    
  
  //! Reads grid parameters from a text file
  void GridInterpolationFull<Dimension2>::ReadText(istream& FileStream)
  {
    int nb;
    switch (this->type_output_file)
      {
      case PLANE :
	FileStream >> xmin >> xmax >> nbPoints_x;
	FileStream >> ymin >> ymax >> nbPoints_y;
	break;
      case PLANE_AXI :
	// plane OAB with O = (xmin, ymin, zmin)
	// A = (xmax, ymax, zmax), B = center
	FileStream >> xmin >> xmax >> ymin >> ymax >> zmin >> zmax;
	FileStream >> center(0) >> center(1) >> center(2) >> nbPoints_x >> nbPoints_y;
	break;
      case VOLUME_AXI :
	FileStream >> xmin >> xmax >> nbPoints_x;
	FileStream >> ymin >> ymax >> nbPoints_y;
	FileStream >> zmin >> zmax >> nbPoints_z;
	break;
      case LINE :
	FileStream >> xmin >> ymin;
	FileStream >> xmax >> ymax >> nb;
	break;
      case LINE_AXI :
	// line joining A=(xmin,ymin,zmin) and B = (xmax, ymax, zmax)
	FileStream >> xmin >> xmax >> ymin >> ymax >> zmin >> zmax;
	FileStream >> nbPoints_x;
	break;
      case POINT :
	FileStream >> xmin >> ymin;
	break;
      case POINT_AXI :
	FileStream >> xmin >> ymin >> zmin;
	break;
      case THREE_PLANES_AXI :
	FileStream >> xmin >> xmax >> nbPoints_x;
	FileStream >> ymin >> ymax >> nbPoints_y;
	FileStream >> zmin >> zmax >> nbPoints_z;
	FileStream >> center(0) >> center(1) >> center(2);
	break;
      case CIRCLE :
	{
	  Real_wp x0, y0, radius;
	  FileStream >> x0 >> y0 >> radius >> nb;
	  xmin = x0 - radius; ymin = y0 - radius;
	  xmax = x0 + radius; ymax = y0 + radius;
	  nbPoints_x = nb;
	}
	break;
      case CIRCLE_AXI :
	FileStream >> center(0) >> center(1) >> center(2) 
		   >> xmin >> ymin >> zmin >> xmax >> ymax >> nbPoints_x;
	break;
      case POINTS_FILE :
	FileStream >> nb;
	nbPoints_x = nb;
	GlobalPoints2D.Reallocate(nb);
	for (int i = 0; i < nb; i++)
	  FileStream >> GlobalPoints2D(i)(0) >> GlobalPoints2D(i)(1);
	
	break;
      case POINTS_FILE_AXI :
	FileStream >> nb;
	nbPoints_x = nb;
	GlobalPoints3D.Reallocate(nb);
	for (int i = 0; i < nb; i++)
	  FileStream >> GlobalPoints3D(i)(0) >> GlobalPoints3D(i)(1) >> GlobalPoints3D(i)(2);
	
	break;
      }
    
  }

  
  //! reads parameters of the grid from a vtk file
  void GridInterpolationFull<Dimension2>::ReadVtk(istream& FileStream,
						  bool& ascii, int& N)
  {
    string header, title, type_data, geometry;
    getline(FileStream, header);
    getline(FileStream, title);
    getline(FileStream, type_data);
    if (type_data == "ASCII")
      ascii = true;
    else if (type_data == "BINARY")
      ascii = false;
    else
      {
	cout << "Expected ASCII/BINARY" << endl;
	abort();
      }
    
    Vector<string> param;
    getline(FileStream, geometry);
    StringTokenize(geometry, param, " \t");
    if (param.GetM() <= 1)
      {
	cout << "Incorrect dataset" << endl;
	abort();
      }
    
    if (param(1) == "STRUCTURED_POINTS")
      {
	// regular points
	// line DIMENSIONS nx ny nz
	FileStream >> title >> nbPoints_x >> nbPoints_y >> nbPoints_z;
	
	// line ORIGIN x0 y0 z0
	FileStream >> title >> xmin >> ymin >> zmin;
	
	// line SPACING dx dy dz
	Real_wp dx, dy, dz;
	FileStream >> title >> dx >> dy >> dz;	
	
	xmax = xmin + dx*(nbPoints_x-1);
	ymax = ymin + dy*(nbPoints_y-1);
	zmax = zmin + dz*(nbPoints_z-1);
	if (nbPoints_z == 1)
	  type_output_file = PLANE;	    
	else
	  type_output_file = VOLUME_AXI;
      }
    else if (param(1) == "STRUCTURED_GRID")
      {
	// regular points
	// line DIMENSIONS nx ny nz
	FileStream >> title >> nbPoints_x >> nbPoints_y >> nbPoints_z;
	
	// line POINTS nb typeData
	FileStream >> title >> N >> type_data;
	bool double_prec = false;
	if (type_data == "double")
	  double_prec = true;
	
	bool store_points = false;
	if (nbPoints_z == 1)
	  {
	    if (nbPoints_y == 1)
	      {		
		if (nbPoints_x == 1)
		  type_output_file = POINT_AXI;
		else
		  {
		    type_output_file = LINE_AXI;
		    store_points = true;
		    GlobalPoints3D.Reallocate(nbPoints_x);
		  }
	      }
	    else
	      type_output_file = PLANE_AXI;
	  }
	
	bool line_on_plane = true;
	Real_wp x, y, z;
	VectReal_wp ptV(3);
	if (!ascii)
	  getline(FileStream, title);
	
	for (int i = 0; i < N; i++)
	  {
	    if (ascii)
	      FileStream >> x >> y >> z;
	    else
	      {
		ReadBinaryDoubleOrFloat(ptV, FileStream, double_prec, false, true);
		x = ptV(0); y = ptV(1); z = ptV(2);
	      }
	    
	    if (z != 0.0)
	      line_on_plane = false;
	    
	    if (store_points)
	      GlobalPoints3D(i).Init(x, y, z);
	    
	    if (i == 0)
	      {
		xmin = x; ymin = y; zmin = z;
	      }
	    else if (i == nbPoints_x-1)
	      {
		xmax = x; ymax = y; zmax = z;
	      }
	    else if (i == nbPoints_x*(nbPoints_y-1))
	      {
		center(0) = x; center(1) = y; center(2) = z;
	      }
	  }
	
	if (line_on_plane)
	  {
	    if (type_output_file == LINE_AXI)
	      {
		R2 vec_u(xmax - xmin, ymax-ymin), vec_v;
		if (Norm2(vec_u) < 1e-6)
		  type_output_file = POINTS_FILE;
		else
		  {
		    // checking if points are on a line		
		    type_output_file = LINE;
		    for (int i = 1; i < nbPoints_x-1; i++)
		      {
			vec_v.Init(GlobalPoints3D(i)(0) - xmin, GlobalPoints3D(i)(1) - ymin);
			Real_wp scal = TimesProd(vec_u, vec_v);
			if (abs(scal) > 1e-6)
			  type_output_file = POINTS_FILE;
		      }
		  }
		
		if (type_output_file == POINTS_FILE)
		  {
		    GlobalPoints2D.Reallocate(nbPoints_x);
		    for (int i = 0; i < nbPoints_x; i++)
		      GlobalPoints2D(i).Init(GlobalPoints3D(i)(0), GlobalPoints3D(i)(1));
		  }

		GlobalPoints3D.Clear();
	      }
	    
	    if (type_output_file == POINT_AXI)
	      type_output_file = POINT;
	  }
	else
	  {
	    if (type_output_file == LINE_AXI)
	      {
		R3 vec_u(xmax - xmin, ymax-ymin, zmax-zmin), vec_v, vec_w;
		if (Norm2(vec_v) < 1e-6)
		  type_output_file = POINTS_FILE_AXI;
		else
		  {
		    for (int i = 1; i < nbPoints_x-1; i++)
		      {
			vec_v.Init(GlobalPoints3D(i)(0) - xmin,
				   GlobalPoints3D(i)(1) - ymin,
				   GlobalPoints3D(i)(2) - zmin);
			
			TimesProd(vec_u, vec_v, vec_w);
			if (Norm2(vec_w) > 1e-6)
			  type_output_file = POINTS_FILE_AXI;
		      }
		  }
		
		if (type_output_file != POINTS_FILE_AXI)
		  GlobalPoints3D.Clear();		
	      }
	  }
      }
    else if (param(1) == "UNSTRUCTURED_GRID")
      {
	// output on three planes
	type_output_file = THREE_PLANES_AXI;
	
	// line POINTS nbPoints
	FileStream >> header >> N >> type_data;
	
	if (!ascii)
	  getline(FileStream, title);
	
	bool double_prec = false;
	if (type_data == "double")
	  double_prec = true;
	
	// points are read
	VectReal_wp ptV(3); Real_wp x, y, z, x0=0., y0=0.;
	bool nz_found = false, nyz_found = false;
	for (int i = 0; i < N; i++)
	  {
	    if (ascii)
	      FileStream >> x >> y >> z;
	    else
	      {
		ReadBinaryDoubleOrFloat(ptV, FileStream, double_prec, false, true);
		x = ptV(0); y = ptV(1); z = ptV(2);
	      }
	    
	    if (i > 0)
	      {
		if (!nz_found)
		  if (y != y0)
		    {
		      nz_found = true;
		      nbPoints_z = i;
		    }
		
		if (!nyz_found)
		  if (x != x0)
		    {
		      center(1) = y;
		      nyz_found = true;
		      nbPoints_y = i/nbPoints_z;
		    }
		
		if (x < xmin)
		  xmin = x;

		if (y < ymin)
		  ymin = y;
		
		if (z < zmin)
		  zmin = z;

		if (x > xmax)
		  xmax = x;

		if (y > ymax)
		  ymax = y;

		if (z > zmax)
		  zmax = z;
	      }
	    else
	      {
		center(0) = x;
		xmin = x; ymin = y; zmin = z;
		xmax = x; ymax = y; zmax = z;
	      }
	    
	    x0 = x; y0 = y;
	  }
	
	center(2) = z;
	nbPoints_x = (N-nbPoints_y*nbPoints_z)/(nbPoints_y+nbPoints_z);
	
	// line CELLS nbCells
	FileStream >> header >> N;

	if (!ascii)
	  getline(FileStream, title);	
	
	IVect num(5);
	int nb_vert, n1, n2, n3, n4;
	for (int i = 0; i < N; i++)
	  {
	    if (ascii)
	      FileStream >> nb_vert >> n1 >> n2 >> n3 >> n4;
	    else
	      FileStream.read(reinterpret_cast<char*>(num.GetData()), 5*sizeof(int));
	  }
	
	// line CELL_TYPES nbCells
	FileStream >> header >> N;

	if (!ascii)
	  getline(FileStream, title);	
	
	for (int i = 0; i < N; i++)
	  {
	    if (ascii)
	      FileStream >> nb_vert >> n1 >> n2 >> n3 >> n4;
	    else
	      FileStream.read(reinterpret_cast<char*>(&nb_vert), sizeof(int));
	  }
      }
    
    // line POINT_DATA nbPoints
    FileStream >> title >> N;
    if (title != "POINT_DATA")
      {
	cout << "CELL_DATA not supported" << endl;
	abort();
      }
  }
  
    
  //! displays details of class GridInterpolationFull<Dimension2>
  template<class Dimension>
  ostream& operator<<(ostream& out, const GridInterpolationFull<Dimension2>& g)
  {
    out<<static_cast<const GridInterpolationFull_Base<Dimension2>& >(g);
    return out;
  }

#endif
  
#ifdef MONTJOIE_WITH_THREE_DIM
  
  /*********************************
   * GridInterpolation<Dimension3> *
   *********************************/
  
  
  //! initialization in order to localize each element on a predefined grid
  void GridInterpolation<Dimension3>::InitInterpolationGrid(const Mesh<Dimension3>& mesh)
  {
    // s will be list of vertices, box the straight polygon bounding element
    // enveloppe the rectangle bounding element
    VectR3 s, box(8), res(8); Matrix3_3 dfj;
    TinyVector<R3, 2> enveloppe;
    SetPoints<Dimension3> PointsElem;
    
    xmin = mesh.GetXmin(); ymin = mesh.GetYmin();
    xmax = mesh.GetXmax(); ymax = mesh.GetYmax();
    zmin = mesh.GetZmin(); zmax = mesh.GetZmax();
    
    // grid step along x, y, z
    Real_wp step_x(0), step_y(0), step_z(0);
    Real_wp deltax = xmax-xmin, deltay = ymax-ymin, deltaz = zmax - zmin,
      delta_glob = abs(deltax)+abs(deltay)+abs(deltaz);
    int ndiv = int(std::pow(double(mesh.GetNbElt()), 1.0/3));
    int nbx = 2*max(toInteger(ceil(deltax/delta_glob*ndiv)), 1);
    int nby = 2*max(toInteger(ceil(deltay/delta_glob*ndiv)), 1);
    int nbz = 2*max(toInteger(ceil(deltaz/delta_glob*ndiv)), 1);
    if (nbx >= 1)
      step_x = deltax/nbx;
    if (nby >= 1)
      step_y = deltay/nby;
    if (nbz >= 1)
      step_z = deltaz/nbz;
    
    nb_subdiv_grid_x = nbx; nb_subdiv_grid_y = nby; nb_subdiv_grid_z = nbz;
    step_subdiv_grid_x = step_x; step_subdiv_grid_y = step_y;
    step_subdiv_grid_z = step_z;
    
    ListeBoxGrid_Element.Reallocate(mesh.GetNbElt());
    // loop on each element
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
	// we get vertices of element i
	mesh.GetVerticesElement(i, s);
	
	// computation of nodal points on the real mesh
	mesh.FjElemNodal(s, PointsElem, mesh, i);
	
	// bounding rectangle in enveloppe and bounding polygon in box
	mesh.GetBoundingBox(s, mesh.GetSafetyCoef(i), box, enveloppe);
	
	// interval along y, where the interpolation points of the grid can be
	int iymin = 0, iymax = nby-1, izmin = 0, izmax = nbz-1,
	  ixmin = 0, ixmax = nbx-1;
	
	// interval along x, y or z, where the interpolation points of the grid can be
	if (abs(step_y) > 0)
	  iymin = max(toInteger(floor((enveloppe(0)(1)-ymin)/step_y)), 0);
	
	if (abs(step_y) > 0)
	  iymax = min(toInteger(ceil((enveloppe(1)(1)-ymin)/step_y)), nby-1);
	
	
	if (abs(step_x) > 0)
	  ixmin = max(toInteger(floor((enveloppe(0)(0)-xmin)/step_x)), 0);
	
	if (abs(step_x) > 0)
	  ixmax = min(toInteger(ceil((enveloppe(1)(0)-xmin)/step_x)),nbx-1);
	
	// interval along x, where the interpolation points of the grid can be
	if (abs(step_z) > 0)
	  izmin = max(toInteger(floor((enveloppe(0)(2)-zmin)/step_z)), 0);
	
	if (abs(step_z) > 0)
	  izmax = min(toInteger(ceil((enveloppe(1)(2)-zmin)/step_z)),nbz-1);
	
	// loop on boxes
	for (int ix = ixmin; ix <= ixmax; ix++)
	  for (int iy = iymin; iy <= iymax; iy++)
	    for (int iz = izmin; iz <= izmax; iz++)
	      {
		int num_box = iz*nby*nbx + iy*nbx + ix;
		ListeBoxGrid_Element(i).PushBack(num_box);
	      }
      }
  }
  
  
  //! pre-localization of points on the regular grid
  /*!
    \param[out] NumBoxGrid_Point grid cell where is each point
    \param[out] ListPoints_Grid for each grid cell, list of points inside
    \param[out] NbPoints_Grid for each grid cell, number of points inside
    \param[in] compute_list if true, ListPoints_Grid and NbPoints_Grid are computed
   */
  void GridInterpolation<Dimension3>::
  PreLocalizePoints(IVect& NumBoxGrid_Point, Vector<IVect>& ListPoints_Grid,
		    IVect& NbPoints_Grid, bool compute_list)
  {
    int N = GlobalCoord.GetM();
    VectR3& Points3D = GlobalCoord;
    
    // points are localized on a regular grid
    int nb_subdiv_grid = nb_subdiv_grid_x*nb_subdiv_grid_y*nb_subdiv_grid_z;
    NumBoxGrid_Point.Reallocate(N); NumBoxGrid_Point.Fill(-1);
    if (compute_list)
      {
	ListPoints_Grid.Reallocate(nb_subdiv_grid);
	NbPoints_Grid.Reallocate(nb_subdiv_grid); NbPoints_Grid.Fill(0);
      }
    
    for (int i = 0; i < Points3D.GetM(); i++)
      {
	Real_wp x = Points3D(i)(0), y = Points3D(i)(1), z = Points3D(i)(2);
	int ix = 0, iy = 0, iz = 0;
	ix = min(toInteger(floor((x-xmin)/step_subdiv_grid_x)), nb_subdiv_grid_x-1);
	iy = min(toInteger(floor((y-ymin)/step_subdiv_grid_y)), nb_subdiv_grid_y-1);
	iz = min(toInteger(floor((z-zmin)/step_subdiv_grid_z)), nb_subdiv_grid_z-1);
	      
	ix = max(0,ix); iy = max(0,iy); iz = max(0,iz);
	NumBoxGrid_Point(i)
	  = iz*nb_subdiv_grid_x*nb_subdiv_grid_y + iy*nb_subdiv_grid_x + ix;
	
	if (compute_list)
	  NbPoints_Grid(NumBoxGrid_Point(i))++;
      }
    
    if (compute_list)
      for (int i = 0; i < nb_subdiv_grid; i++)
	{
	  ListPoints_Grid(i).Reallocate(NbPoints_Grid(i));
	  ListPoints_Grid(i).Fill(-1);
	}
    
    if (compute_list)
      NbPoints_Grid.Fill(0);
    
    if (compute_list)
      for (int i = 0; i < Points3D.GetM(); i++)
	{
	  int num_box = NumBoxGrid_Point(i);
	  if (num_box >= 0)
	    ListPoints_Grid(num_box)(NbPoints_Grid(num_box)++) = i;
	}
  }
  
  
  //! localization of points on the grid, using pre-computed localization of elements
  void GridInterpolation<Dimension3>::LocalizePoints(const Mesh<Dimension3>& mesh)
  {
    LocalizePoints(mesh, this->GlobalCoord, this->ElementInterp, this->CoorInterp, 
		   this->Epsilon_InsideElement, this->point_treated, this->dfjm1, this->dfjm1_store);
  }

  
  //! Localization of points2D in the mesh
  void GridInterpolation<Dimension3>::
  LocalizePoints(const Mesh<Dimension3>& mesh, const VectR3& Points3D,
		 IVect& ElementInt, VectR3& CoorInt, VectReal_wp& EpsilonInt,
		 Vector<bool>& point_treat, Vector<Matrix3_3>& mat_dfjm1, bool store_df)
  {
    int N = Points3D.GetM();
    
    // we initialize the array, which will contain number of elements
    // where each interpolation is.
    ElementInt.Reallocate(N);
    // initialization of local coordinates also
    CoorInt.Reallocate(N);
    ElementInt.Fill(-1);
    CoorInt.Fill(R3(0,0,0));
    
    // resizing array point_treated
    int Nold = point_treat.GetM();
    point_treat.Resize(N);
    for (int i = Nold; i < N; i++)
      point_treat(i) = false;

    EpsilonInt.Reallocate(N);
    EpsilonInt.Fill(-1);    
    if (store_df)
      mat_dfjm1.Reallocate(N);
    
    // points are localized on a regular grid
    int nb_subdiv_grid = this->nb_subdiv_grid_x*this->nb_subdiv_grid_y*this->nb_subdiv_grid_z;
    IVect NumBoxGrid_Point(N); NumBoxGrid_Point.Fill(-1);
    Vector<IVect> ListPoints_Grid(nb_subdiv_grid);
    IVect NbPoints_Grid(nb_subdiv_grid); NbPoints_Grid.Fill(0);
    for (int i = 0; i < Points3D.GetM(); i++)
      if (!point_treat(i))
	{
	  Real_wp x = Points3D(i)(0), y = Points3D(i)(1), z = Points3D(i)(2);
	  int ix = 0, iy = 0, iz = 0;
	  if ((x > xmax+R_N::threshold)||(x < xmin-R_N::threshold)
	      ||(y > ymax+R_N::threshold)||(y < ymin-R_N::threshold)
	      ||(z > zmax+R_N::threshold)||(z < zmin-R_N::threshold))
	    {
	      NumBoxGrid_Point(i) = -1;
	    }
	  else
	    {
	      if (abs(this->step_subdiv_grid_x) > R_N::threshold)
		ix = min(toInteger(floor((x-xmin)/this->step_subdiv_grid_x)),
			 this->nb_subdiv_grid_x-1);
	      
	      if (abs(this->step_subdiv_grid_y) > R_N::threshold)
		iy = min(toInteger(floor((y-ymin)/this->step_subdiv_grid_y)),
			 this->nb_subdiv_grid_y-1);
	      
	      if (abs(this->step_subdiv_grid_z) > R_N::threshold)
		iz = min(toInteger(floor((z-zmin)/this->step_subdiv_grid_z)),
			 this->nb_subdiv_grid_z-1);
	      
	      ix = max(0, ix); iy = max(0, iy); iz = max(0, iz);
	      NumBoxGrid_Point(i) 
		= iz*this->nb_subdiv_grid_x*this->nb_subdiv_grid_y + iy*this->nb_subdiv_grid_x + ix;
	      
	      NbPoints_Grid(NumBoxGrid_Point(i))++;
	    }
	  
	  point_treat(i) = true;
	}
    
    for (int i = 0; i < nb_subdiv_grid; i++)
      {
	ListPoints_Grid(i).Reallocate(NbPoints_Grid(i));
	ListPoints_Grid(i).Fill(-1);
      }
    
    NbPoints_Grid.Fill(0);
    for (int i = 0; i < Points3D.GetM(); i++)
      {
	int num_box = NumBoxGrid_Point(i);
	if (num_box >= 0)
	  ListPoints_Grid(num_box)(NbPoints_Grid(num_box)++) = i;
      }
    
    R3 pointA, point_loc;
    int old_percent = 0, new_percent = 0;
    TinyVector<R3, 2> enveloppe; Matrix3_3 dfj;
    
    int nb_points_tested = 0, nb_points_inverted = 0;    
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
	// displays a progress bar
	new_percent = toInteger(round(Real_wp(i)/(mesh.GetNbElt()-1)*80));
        if (mesh.print_level >= 2)
          for (int percent = old_percent; percent < new_percent; percent++)
            { cout<<"#"; cout.flush(); }
	
	old_percent = new_percent;
	
	// first, we look for any point near of the element
	bool near_point = false;
	for (int j = 0; j < ListeBoxGrid_Element(i).GetM(); j++)
	  if (NbPoints_Grid(ListeBoxGrid_Element(i)(j)) > 0)
	    near_point = true;
	
	if (near_point)
	  {
	    FjInverseProblem<Dimension3> inverseFj(mesh, i);
	    VectR3& s = inverseFj.GetVertices();
	    SetPoints<Dimension3>& PointsElem = inverseFj.GetSetPoints();
	    
	    // we use nodal points to compute the enveloppe
	    mesh.GetBoundingBox(i, s, PointsElem, enveloppe);
	    
	    // loop on each box containing element i
	    for (int j = 0; j < ListeBoxGrid_Element(i).GetM(); j++)
	      {
		int num_box = ListeBoxGrid_Element(i)(j);
		for (int n = 0; n < NbPoints_Grid(num_box); n++)
		  {
		    nb_points_tested++; 
		    // number of the interpolation point
		    int num_point = ListPoints_Grid(num_box)(n);
		    pointA = Points3D(num_point);
		    
		    // point inside the enveloppe ?
		    if (PointInsideBoundingBox(pointA, enveloppe))
		      {
			
			nb_points_inverted++;
			// Is this interpolation point inside the element i ?
			// point_loc = F_i^{-1} (pointA) is computed 
			bool test_inside = inverseFj.Solve(pointA, point_loc);
			Real_wp distance_boundary_elt
			  = mesh.GetDistanceToBoundary(point_loc, i);
			
			if (distance_boundary_elt < 0)
			  mesh.ProjectPointOnBoundary(point_loc, i);
			
			if (distance_boundary_elt < EpsilonInt(num_point))
			  test_inside = false;
			
			if (test_inside)
			  {
			    // Yes, it is
			    // we update ElementInterp and CoorInterp
			    ElementInt(num_point) = i;
			    CoorInt(num_point) = point_loc;
			    EpsilonInt(num_point) = distance_boundary_elt;
			    if (store_df)
			      {
				mesh.DFj(s, PointsElem, point_loc, dfj, mesh, i);
				GetInverse(dfj, mat_dfjm1(num_point));
			      }
			  }
		      }
		  }
	      }
	  }
      }
    
    if (mesh.print_level >= 2)
      {
        cout<<endl;
        cout << "Number of tested points : " << nb_points_tested << endl;
        cout << "Number of points for which F_i has been inverted : " << nb_points_inverted << endl;
        cout << "Average number of iterations for non-linear solver : "
             << Real_wp(FjInverseProblem<Dimension3>::nb_iter_all) / nb_points_inverted << endl;
      }
  }

  
  //! displays details of class GridInterpolation<Dimension3>
  template<class Dimension>
  ostream& operator<<(ostream& out, const GridInterpolation<Dimension3>& g)
  {
    out<<static_cast<const GridInterpolation_Base<Dimension3>& >(g);
    return out;
  }

  
  /*************************************
   * GridInterpolationFull<Dimension3> *
   *************************************/
  
  
  //! reading of a line of the data file
  /*!
    \param[in] description_field keyword of the line
    \param[out] parameters values of the line
    \param[in] nb_param number of values
    this function fills correct parameters depending
    on the matching line of the data file
   */
  void GridInterpolationFull<Dimension3>::
  SetInputData(const string& description_field,
	       const VectString& parameters)
  {
    int nb = 0;
    if (!description_field.compare("SismoGrille3D"))
      {
	if (parameters.GetM() <= 3)
	  {
	    cout << "In SetInputData of GridInterpolation" << endl;
	    cout << "SismoGrille3D needs more parameters, for instance :" << endl;
	    cout << "SismoGrille3D = AUTO Nx Ny Nz" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	type_output_file = VOLUME;
	if (!parameters(nb).compare("AUTO"))
	  {
	    grid_interval_to_be_computed = true;
	    nb++;
	  }
	else
	  {
	    if (parameters.GetM() <= 8)
	      {
		cout << "In SetInputData of GridInterpolation" << endl;
		cout << "SismoGrille3D needs more parameters, for instance :" << endl;
		cout << "SismoGrille3D = xmin xmax ymin ymax zmin zmax Nx Ny Nz" << endl;
		cout << "Current parameters are : " << endl << parameters << endl;
		abort();
	      }

	    xmin = to_num<Real_wp>(parameters(nb++));
	    xmax = to_num<Real_wp>(parameters(nb++));
	    ymin = to_num<Real_wp>(parameters(nb++));
	    ymax = to_num<Real_wp>(parameters(nb++));
	    zmin = to_num<Real_wp>(parameters(nb++));
	    zmax = to_num<Real_wp>(parameters(nb++));
	  }
	
	nbPoints_x = to_num<int>(parameters(nb++));
	nbPoints_y = to_num<int>(parameters(nb++));
	nbPoints_z = to_num<int>(parameters(nb++));
      }
    else if (!description_field.compare("SismoGrille"))
      {
	if (parameters.GetM() <= 6)
	  {
	    cout << "In SetInputData of GridInterpolation" << endl;
	    cout << "SismoGrille needs more parameters, for instance :" << endl;
	    cout << "SismoGrille = x0 y0 z0 AUTO Nx Ny Nz" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	type_output_file = THREE_PLANES;
	center(0) = to_num<Real_wp>(parameters(nb++));
	center(1) = to_num<Real_wp>(parameters(nb++));
	center(2) = to_num<Real_wp>(parameters(nb++));
	if (!parameters(nb).compare("AUTO"))
	  {
	    grid_interval_to_be_computed = true;
	    nb++;
	  }
	else
	  {
	    if (parameters.GetM() <= 11)
	      {
		cout << "In SetInputData of GridInterpolation" << endl;
		cout << "SismoGrille needs more parameters, for instance :" << endl;
		cout << "SismoGrille = x0 y0 z0 xmin xmax ymin ymax zmin zmax Nx Ny Nz" << endl;
		cout << "Current parameters are : " << endl << parameters << endl;
		abort();
	      }

	    xmin = to_num<Real_wp>(parameters(nb++));
	    xmax = to_num<Real_wp>(parameters(nb++));
	    ymin = to_num<Real_wp>(parameters(nb++));
	    ymax = to_num<Real_wp>(parameters(nb++));
	    zmin = to_num<Real_wp>(parameters(nb++));
	    zmax = to_num<Real_wp>(parameters(nb++));
	  }
	
	nbPoints_x = to_num<int>(parameters(nb++));
	nbPoints_y = to_num<int>(parameters(nb++));
	nbPoints_z = to_num<int>(parameters(nb++));
      }
    else if (!description_field.compare("SismoPlane"))
      {
	if (parameters.GetM() <= 10)
	  {
	    cout << "In SetInputData of GridInterpolation" << endl;
	    cout << "SismoPlane needs more parameters, for instance :" << endl;
	    cout << "SismoPlane = xmin ymin zmin xmax ymax zmax x0 y0 z0 Nx Ny" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	type_output_file = PLANE;
	xmin = to_num<Real_wp>(parameters(nb++));
	ymin = to_num<Real_wp>(parameters(nb++));
	zmin = to_num<Real_wp>(parameters(nb++));
	xmax = to_num<Real_wp>(parameters(nb++));
	ymax = to_num<Real_wp>(parameters(nb++));
	zmax = to_num<Real_wp>(parameters(nb++));
	center(0) = to_num<Real_wp>(parameters(nb++));
	center(1) = to_num<Real_wp>(parameters(nb++));
	center(2) = to_num<Real_wp>(parameters(nb++));
	nbPoints_x = to_num<int>(parameters(nb++));
	nbPoints_y = to_num<int>(parameters(nb++));
	nbPoints_z = 1;
      }
    else if (!description_field.compare("SismoCircle"))
      {
	if (parameters.GetM() <= 8)
	  {
	    cout << "In SetInputData of GridInterpolation" << endl;
	    cout << "SismoCircle needs more parameters, for instance :" << endl;
	    cout << "SismoCircle = x0 y0 z0 nx ny nz rx ry N" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	type_output_file = CIRCLE;
	center(0) = to_num<Real_wp>(parameters(nb++));
	center(1) = to_num<Real_wp>(parameters(nb++));
	center(2) = to_num<Real_wp>(parameters(nb++));
        // normale to the plane where the circle lies
	xmin = to_num<Real_wp>(parameters(nb++));
	ymin = to_num<Real_wp>(parameters(nb++));
	zmin = to_num<Real_wp>(parameters(nb++));
	// radius of the circle
	xmax = to_num<Real_wp>(parameters(nb++));
	ymax = to_num<Real_wp>(parameters(nb++));	
	nbPoints_x = to_num<int>(parameters(nb++));
	nbPoints_y = 1;
      }
    else if (!description_field.compare("SismoLine"))
      {
	if (parameters.GetM() <= 6)
	  {
	    cout << "In SetInputData of GridInterpolation" << endl;
	    cout << "SismoLine needs more parameters, for instance :" << endl;
	    cout << "SismoLine = xmin xmax ymin ymax zmin zmax N" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	type_output_file = LINE;
	xmin = to_num<Real_wp>(parameters(nb++));
	xmax = to_num<Real_wp>(parameters(nb++));
	ymin = to_num<Real_wp>(parameters(nb++));
	ymax = to_num<Real_wp>(parameters(nb++));
	zmin = to_num<Real_wp>(parameters(nb++));
	zmax = to_num<Real_wp>(parameters(nb++));
	nbPoints_x = to_num<int>(parameters(nb++));
	nbPoints_y = 1; nbPoints_z = 1;
      }
    else if (!description_field.compare("SismoPoint"))
      {
	if (parameters.GetM() <= 2)
	  {
	    cout << "In SetInputData of GridInterpolation" << endl;
	    cout << "SismoPoint needs more parameters, for instance :" << endl;
	    cout << "SismoPoint = x0 y0 z0" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	type_output_file = POINT;
	xmin = to_num<Real_wp>(parameters(nb++));
	ymin = to_num<Real_wp>(parameters(nb++));
	zmin = to_num<Real_wp>(parameters(nb++));
	xmax = xmin; ymax = ymin; zmax = zmin;
	nbPoints_x = 1; nbPoints_y = 1; nbPoints_z = 1;
      }
    else if (!description_field.compare("SismoPointsFile"))
      {
	if (parameters.GetM() <= 0)
	  {
	    cout << "In SetInputData of GridInterpolation" << endl;
	    cout << "SismoPointsFile needs more parameters, for instance :" << endl;
	    cout << "SismoPointsFile = " << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	type_output_file = POINTS_FILE ;
	PointsInputFile = parameters(0);
      }
  }
  
  //! returns type of regular grid (3-d grid, plane, line, point, etc)
  int GridInterpolationFull<Dimension3>::GetType(const string& description)
  {
    if (!description.compare("FileOutputGrille3D"))
      return VOLUME;
    else if (!description.compare("FileOutputGrille"))
      return THREE_PLANES;
    else if (!description.compare("FileOutputPlane"))
      return PLANE;
    else if (!description.compare("FileOutputLine"))
      return LINE;
    else if (!description.compare("FileOutputPoint"))
      return POINT;
    else if (!description.compare("FileOutputPointsFile"))
      return POINTS_FILE;
    else if (!description.compare("FileOutputCircle"))
      return CIRCLE;
    
    return -1;
  }

  
  //! returns the dimension of the output grid
  int GridInterpolationFull<Dimension3>::GetDimension(int type)
  {
    return 3;
  }
  
  
  //! initializes grid with regular points on a 3-D plane
  void GridInterpolationFull<Dimension3>::
  SetPlaneOutput(const R3& ptO, const R3& ptA, const R3& ptB,
		    int nbx, int nby)
  {
    this->type_output_file = PLANE;
    this->xmin = ptO(0); this->ymin = ptO(1); this->zmin = ptO(2);
    this->xmax = ptA(0); this->ymax = ptA(2); this->zmax = ptA(2);
    this->center = ptB;
    this->nbPoints_x = nbx;
    this->nbPoints_y = nby;
    this->nbPoints_z = 1;
  }
  
  
  //! initialization of the grid with regular points on a 3-D line
  void GridInterpolationFull<Dimension3>::
  SetLineOutput(const R3& ptA, const R3& ptB, int nbx)
  {
    this->type_output_file = LINE;
    this->xmin = ptA(0); this->xmax = ptB(0);
    this->ymin = ptA(1); this->ymax = ptB(1);
    this->zmin = ptA(2); this->zmax = ptB(2);
    this->nbPoints_x = nbx;    
    this->nbPoints_y = 1;
    this->nbPoints_z = 1;
  }
  
  
  //! initialization of the grid with a single 3-D point
  void GridInterpolationFull<Dimension3>::
  SetPointOutput(const R3& ptA)
  {
    this->type_output_file = POINT;
    this->xmin = ptA(0); this->ymin = ptA(1); this->zmin = ptA(2);
    this->nbPoints_x = 1;
    this->nbPoints_y = 1;
    this->nbPoints_z = 1;
  }
  

  //! initialization of the grid with regular points on a 3-D circle
  void GridInterpolationFull<Dimension3>::
  SetCircleOutput(const R3& ptA, const R3& normale,
		  const Real_wp& rx, const Real_wp& ry, int nbx)
  {
    this->type_output_file = CIRCLE;
    this->center = ptA;
    this->xmin = normale(0); this->ymin = normale(1); this->zmin = normale(2);
    this->xmax = rx; this->ymax = ry;
    this->nbPoints_x = nbx;
    this->nbPoints_y = 1;
    this->nbPoints_z = 1;
  }

  
  //! initialization of the grid with 3-D points stored in a file
  void GridInterpolationFull<Dimension3>::SetPointsFileOutput(const string& name)
  {
    this->type_output_file = POINTS_FILE;
    PointsInputFile = name;
  }
  
  
  //! initialization of the grid with regular points on three planes
  void GridInterpolationFull<Dimension3>::
  SetThreePlanesOutput(const R3& ptO, const R3& ptMin, const R3& ptMax,
		       int nbx_, int nby_, int nbz_)
  {
    this->type_output_file = THREE_PLANES;
    this->center = ptO;
    this->xmin = ptMin(0); this->ymin = ptMin(1); this->zmin = ptMin(2);
    this->xmax = ptMax(0); this->ymax = ptMax(1); this->zmax = ptMax(2);
    nbPoints_x = nbx_; nbPoints_y = nby_; nbPoints_z = nbz_;
  }
  
  
  //! initialization of the grid with regular points on a parallelepiped
  void GridInterpolationFull<Dimension3>::
  SetVolumeOutput(const R3& ptMin, const R3& ptMax,
		  int nbx_, int nby_, int nbz_)
  {
    this->type_output_file = VOLUME;
    this->xmin = ptMin(0); this->ymin = ptMin(1); this->zmin = ptMin(2);
    this->xmax = ptMax(0); this->ymax = ptMax(1); this->zmax = ptMax(2);
    nbPoints_x = nbx_; nbPoints_y = nby_; nbPoints_z = nbz_;
  }


  //! computation of 3-D points
  void GridInterpolationFull<Dimension3>
  ::GenerateGridPoints(VectR3& Points3D)
  {
    GridInterpolation<Dimension3> grid; IVect Index, NumPartie;
    grid.SetXmin(-1e300); grid.SetXmax(1e300);
    grid.SetYmin(-1e300); grid.SetYmax(1e300);
    grid.SetZmin(-1e300); grid.SetZmax(1e300);

    VectR3 PointsInit;
    GenerateGridPoints(grid, Index, PointsInit, NumPartie);

    Points3D.Reallocate(PointsInit.GetM());
    for (int i = 0; i < Points3D.GetM(); i++)
      Points3D(Index(i)) = PointsInit(i);
  }
  

  //! computation of 3-D points
  void GridInterpolationFull<Dimension3>
  ::GenerateGridPoints(const GridInterpolation<Dimension3>& grid, IVect& Index,
		       VectR3& Points3D, IVect& NumPartie)
  {
    R3 point, pt3D; int n;
    Real_wp step_x(0), step_y(0), step_z(0);
    if (nbPoints_x > 1)
      step_x = (xmax-xmin)/(nbPoints_x-1);
    if (nbPoints_y > 1)
      step_y = (ymax-ymin)/(nbPoints_y-1);
    if (nbPoints_z > 1)
      step_z = (zmax-zmin)/(nbPoints_z-1);
    
    if (type_output_file == VOLUME)
      {
	// 3-D regular cartesian grid
	int nbPoints_grid = nbPoints_x*nbPoints_y*nbPoints_z;
	int num = 0;
	for (int i = 0; i < nbPoints_x; i++)
	  for (int j = 0; j < nbPoints_y; j++)
	    for (int k = 0; k < nbPoints_z; k++)
	      {
		point(0) = xmin + i*step_x;
		point(1) = ymin + j*step_y;
		point(2) = zmin + k*step_z;
		if (TranslatePoint(point, grid, pt3D, n))
		  num++;
	      }

	nb_global_grid_points = nbPoints_grid;

	Points3D.Reallocate(num);
	Index.Reallocate(num);
	NumPartie.Reallocate(num);
	num = 0;
	for (int i = 0; i < nbPoints_x; i++)
	  for (int j = 0; j < nbPoints_y; j++)
	    for (int k = 0; k < nbPoints_z; k++)
	      {
		int nb = k*nbPoints_x*nbPoints_y + j*nbPoints_x + i;
		point(0) = xmin + i*step_x;
		point(1) = ymin + j*step_y;
		point(2) = zmin + k*step_z;
		if (TranslatePoint(point, grid, pt3D, n))
		  {
		    Index(num) = nb;
		    Points3D(num) = pt3D;
		    NumPartie(num) = n;
		    num++;
		  }
	      }	
      }
    else if (type_output_file == THREE_PLANES)
      {
	// 2-D cartesian grid on three planes
	// x = center(0), y = center(1) and z = center(2)
	int nbPoints_grid = nbPoints_x*nbPoints_y + nbPoints_x*nbPoints_z
	  + nbPoints_y*nbPoints_z;
	
	int num = 0;
	for (int j = 0; j < nbPoints_y; j++)
	  for (int k = 0; k < nbPoints_z; k++)
	    {
	      point(0) = center(0);
	      point(1) = ymin + j*step_y;
	      point(2) = zmin + k*step_z;
	      if (TranslatePoint(point, grid, pt3D, n))
		num++;
	    }
	
	for (int i = 0; i < nbPoints_x; i++)
	  for (int k = 0; k < nbPoints_z; k++)
	    {
	      point(0) = xmin + i*step_x;
	      point(1) = center(1);
	      point(2) = zmin + k*step_z;
	      if (TranslatePoint(point, grid, pt3D, n))
		num++;
	    }
	
	for (int i = 0; i < nbPoints_x; i++)
	  for (int j = 0; j < nbPoints_y; j++)
	    {
	      point(0) = xmin + i*step_x;
	      point(1) = ymin + j*step_y;
	      point(2) = center(2);
	      if (TranslatePoint(point, grid, pt3D, n))
		num++;
	    }
	
	nb_global_grid_points = nbPoints_grid;

	Points3D.Reallocate(num);
	Index.Reallocate(num);
	NumPartie.Reallocate(num);
	num = 0;

	for (int j = 0; j < nbPoints_y; j++)
	  for (int k = 0; k < nbPoints_z; k++)
	    {
	      int nb = j*nbPoints_z + k;
	      point(0) = center(0);
	      point(1) = ymin + j*step_y;
	      point(2) = zmin + k*step_z;
	      if (TranslatePoint(point, grid, pt3D, n))
		{
		  Index(num) = nb;
		  Points3D(num) = pt3D;
		  NumPartie(num) = n;
		  num++;
		}
	    }
	
	for (int i = 0; i < nbPoints_x; i++)
	  for (int k = 0; k < nbPoints_z; k++)
	    {
	      int nb = nbPoints_z*nbPoints_y + i*nbPoints_z + k;
	      point(0) = xmin + i*step_x;
	      point(1) = center(1);
	      point(2) = zmin + k*step_z;
	      if (TranslatePoint(point, grid, pt3D, n))
		{
		  Index(num) = nb;
		  Points3D(num) = pt3D;
		  NumPartie(num) = n;
		  num++;
		}
	    }
	
	for (int i = 0; i < nbPoints_x; i++)
	  for (int j = 0; j < nbPoints_y; j++)
	    {
	      int nb = nbPoints_z*(nbPoints_x+nbPoints_y) + i*nbPoints_y + j;
	      point(0) = xmin + i*step_x;
	      point(1) = ymin + j*step_y;
	      point(2) = center(2);
	      if (TranslatePoint(point, grid, pt3D, n))
		{
		  Index(num) = nb;
		  Points3D(num) = pt3D;
		  NumPartie(num) = n;
		  num++;
		}
	    }
      }
    else if (type_output_file == PLANE)
      {
	// 2-D cartesian grid on a plane
	R3 vec_u, vec_v;
	vec_u(0) = xmax - xmin; vec_u(1) = ymax - ymin; vec_u(2) = zmax - zmin;
	vec_v(0) = center(0) - xmin; vec_v(1) = center(1) - ymin;
	vec_v(2) = center(2) - zmin;
	step_x = 0.0; step_y = 0.0;
	if (nbPoints_x > 1)
	  step_x = 1.0/(nbPoints_x-1);
	if (nbPoints_y > 1)
	  step_y = 1.0/(nbPoints_y-1);
	
	int num = 0;
	for (int i = 0; i < nbPoints_x; i++)
	  for (int j = 0; j < nbPoints_y; j++)
	    {
	      point(0) = xmin + i*step_x*vec_u(0) + j*step_y*vec_v(0);
	      point(1) = ymin + i*step_x*vec_u(1) + j*step_y*vec_v(1);
	      point(2) = zmin + i*step_x*vec_u(2) + j*step_y*vec_v(2);
	      if (TranslatePoint(point, grid, pt3D, n))
		num++;
	    }
	
		
	nb_global_grid_points = nbPoints_x * nbPoints_y;

	Points3D.Reallocate(num);
	Index.Reallocate(num);
	NumPartie.Reallocate(num);
	num = 0;

	for (int i = 0; i < nbPoints_x; i++)
	  for (int j = 0; j < nbPoints_y; j++)
	    {
	      int nb = j*nbPoints_x + i;
	      point(0) = xmin + i*step_x*vec_u(0) + j*step_y*vec_v(0);
	      point(1) = ymin + i*step_x*vec_u(1) + j*step_y*vec_v(1);
	      point(2) = zmin + i*step_x*vec_u(2) + j*step_y*vec_v(2);
	      if (TranslatePoint(point, grid, pt3D, n))
		{
		  Index(num) = nb;
		  Points3D(num) = pt3D;
		  NumPartie(num) = n;
		  num++;
		}
	    }
      }
    else if (type_output_file == CIRCLE)
      {
	Real_wp step_teta = 2.0*pi_wp/nbPoints_x;
	R3 normale, vec_u, vec_v; Real_wp radius_x = xmax, radius_y = ymax;
	normale(0) = xmin; normale(1) = ymin; normale(2) = zmin;
	GetVectorPlane(normale, vec_u, vec_v);

	int num = 0;
	for (int i = 0; i < nbPoints_x; i++)
	  {
	    Real_wp teta = step_teta*i;
	    point = center; Add(radius_x*cos(teta), vec_u, point);
	    Add(radius_y*sin(teta), vec_v, point);
	    if (TranslatePoint(point, grid, pt3D, n))
	      num++;
	  }
        
	nb_global_grid_points = nbPoints_x;

	Points3D.Reallocate(num);
	Index.Reallocate(num);
	NumPartie.Reallocate(num);
	num = 0;

	for (int i = 0; i < nbPoints_x; i++)
	  {
	    Real_wp teta = step_teta*i;
	    point = center; Add(radius_x*cos(teta), vec_u, point);
	    Add(radius_y*sin(teta), vec_v, point);
	    if (TranslatePoint(point, grid, pt3D, n))
	      {
		Index(num) = i;
		Points3D(num) = pt3D;
		NumPartie(num) = n;
		num++;
	      }
	  }
      }
    else if (type_output_file == LINE)
      {
	// 1-D grid on a line
	step_x = 0.0; step_y = 0.0; step_z = 0.0;
	if (nbPoints_x > 1)
	  {
	    step_x = (xmax-xmin)/(nbPoints_x-1);
	    step_y = (ymax-ymin)/(nbPoints_x-1);
	    step_z = (zmax-zmin)/(nbPoints_x-1);
	  }
	
	int num = 0;
	for (int i = 0; i < nbPoints_x; i++)
	  {
	    point(0) = xmin + i*step_x;
	    point(1) = ymin + i*step_y;
	    point(2) = zmin + i*step_z;
	    if (TranslatePoint(point, grid, pt3D, n))
	      num++;
	  }
	
	nb_global_grid_points = nbPoints_x;

	Points3D.Reallocate(num);
	Index.Reallocate(num);
	NumPartie.Reallocate(num);
	num = 0;

	for (int i = 0; i < nbPoints_x; i++)
	  {
	    point(0) = xmin + i*step_x;
	    point(1) = ymin + i*step_y;
	    point(2) = zmin + i*step_z;
	    if (TranslatePoint(point, grid, pt3D, n))
	      {
		Index(num) = i;
		Points3D(num) = pt3D;
		NumPartie(num) = n;
		num++;
	      }
	  }
      }
    else if (type_output_file == POINT)
      {
	// output on single point
	point.Init(xmin, ymin, zmin);
	nb_global_grid_points = 1;
	if (TranslatePoint(point, grid, pt3D, n))
	  {
	    Points3D.Reallocate(1);
	    Index.Reallocate(1);
	    NumPartie.Reallocate(1);
	    
	    Index(0) = 0;
	    Points3D(0) = pt3D;
	    NumPartie(0) = n;
	  }
	else
	  {
	    Index.Clear();
	    Points3D.Clear();
	    NumPartie.Clear();
	  }
      }
    else if  (this->type_output_file == POINTS_FILE)
      {
	xmin = 0; xmax = 1.0; nbPoints_x = 10;
	ymin = 0; ymax = 1.0; nbPoints_y = 10;
	
	GlobalPoints3D.ReadText(PointsInputFile);
	Mlt(1.0/z0_adim, GlobalPoints3D);
	
	int num = 0;
	for (int i = 0; i < GlobalPoints3D.GetM(); i++)
	  {
	    point = GlobalPoints3D(i);
	    if (TranslatePoint(point, grid, pt3D, n))
	      num++;
	  }
	
	nb_global_grid_points = GlobalPoints3D.GetM();
	nbPoints_x = nb_global_grid_points;
	
	Points3D.Reallocate(num);
	Index.Reallocate(num);
	NumPartie.Reallocate(num);
	num = 0;
	for (int i = 0; i < GlobalPoints3D.GetM(); i++)
	  {
	    point = GlobalPoints3D(i);
	    if (TranslatePoint(point, grid, pt3D, n))
	      {
		Index(num) = i;
		Points3D(num) = pt3D;
		NumPartie(num) = n;
		num++;
	      }
	  }
	
#ifdef SELDON_WITH_MPI
        int rank_proc; MPI_Comm_rank(grid.GetCommunicator(), &rank_proc);
	if (rank_proc != 0)
	  GlobalPoints3D.Clear();
#endif
      }
  }
  
  
  //! computation of global points if it is a regular grid, plane, line, etc
  /*!
    \param[in] grid object used to localize points on the mesh
    this methods computes global points in the case of regular grid, plane, 
    line divided in equal subdivisions, and appends these points to the object grid.
   */
  void GridInterpolationFull<Dimension3>::
  InitGrid(GridInterpolation<Dimension3>& grid)
  {
    this->var_interp = &grid;    
    
    if (this->grid_interval_to_be_computed)
      {
	xmin = xmin0; ymin = ymin0; zmin = zmin0;
	xmax = xmax0; ymax = ymax0; zmax = zmax0;

        xmin += 5.0*(xmax-xmin)*epsilon_machine;
        xmax -= 5.0*(xmax-xmin)*epsilon_machine;
        ymin += 5.0*(ymax-ymin)*epsilon_machine;
        ymax -= 5.0*(ymax-ymin)*epsilon_machine;
        zmin += 5.0*(zmax-zmin)*epsilon_machine;
        zmax -= 5.0*(zmax-zmin)*epsilon_machine;
      }
    
    VectR3 Points3D; IVect NumPartie;
    GenerateGridPoints(grid, IndexPoints, Points3D, NumPartie);
    
    int N = Points3D.GetM();
    VectReal_wp false_teta(N);
    false_teta.Zero();
    
    grid.Append(Points3D, false_teta, NumPartie, list_points);
  }
  

  //! translates a point to the physical domain and returns true if 
  //! the translated point is inside the grid
  bool GridInterpolationFull<Dimension3>::
  TranslatePoint(const R3& point, const GridInterpolation<Dimension3>& grid,
	         R3& pt3D, int& num_partie)
  {    	
    pt3D = point;
    int nb_sections = grid.GetNbCyclicSections();
    if (nb_sections > 1) 
      {
        int nx, ny, nz;
        grid.GetNbCartesianSections(nx, ny, nz);
        if (nz < 1)
          nz = 1;
	
	// cyclic domain
	Real_wp theta = 0, thetabis , r;
        Real_wp dz = this->zmax0 - this->zmin0;
        int posz = (nz+1)/2;
        Real_wp z_max_ = this->zmin0 + posz*dz;
        Real_wp z_min_ = z_max_ - nz*dz;
	int n;

	CartesianToPolar(point(0), point(1), r, theta);
	
	if (theta < 0)
	  theta += 2.0*pi_wp;
	
	Real_wp z = point(2);            	    
	Real_wp nreal = theta*Real_wp(nb_sections)/(2.0*pi_wp);
	if (abs(round(nreal) - nreal) < 1e-10)
	  nreal = round(nreal);
	
	n = toInteger(floor(nreal));
        
	thetabis = theta - ( n * (2.0*pi_wp ) / Real_wp(nb_sections)) ;
                        
	pt3D(0) = r * cos(thetabis) ;
	pt3D(1) = r * sin(thetabis) ;
	
	Real_wp iz_ = (z - z_min_)/dz;
	int iz = toInteger(round(iz_));
	if (abs(iz_ - iz) > 1e-12)
	  iz = toInteger(floor(iz_));
	
	if ((iz < 0) || (iz >= nz))
	  return false;
        
	pt3D(2) = this->zmin0 + (z - z_min_ - dz*iz);
	num_partie = nz*n + iz;
      }
    else
      {
        // periodicity with respect to x and y ?
        int nx, ny, nz;
        grid.GetNbCartesianSections(nx, ny, nz);
        if (nx < 1)
          nx = 1;
        
        if (ny < 1)
          ny = 1;

        if (nz < 1)
          nz = 1;
        
        if ((nx != 1) || (ny != 1) || (nz != 1))
          {
            // Delta x and Delta y for the initial cell
            Real_wp dx = this->xmax0 - this->xmin0;
            Real_wp dy = this->ymax0 - this->ymin0;
            Real_wp dz = this->zmax0 - this->zmin0;
            
            // position of the central cell
            int posx = (nx+1)/2;
            int posy = (ny+1)/2;
            int posz = (nz+1)/2;
            
            // extremum values of the supercell
            Real_wp x_max_ = this->xmin0 + posx*dx;
            Real_wp y_max_ = this->ymin0 + posy*dy;
            Real_wp z_max_ = this->zmin0 + posz*dz;
            
            Real_wp x_min_ = x_max_ - nx*dx;
            Real_wp y_min_ = y_max_ - ny*dy;
            Real_wp z_min_ = z_max_ - nz*dz;
            
	    Real_wp x = point(0);
	    Real_wp y = point(1);
	    Real_wp z = point(2);
	    
	    // we detect on which cell the point is located
	    Real_wp ix_ = (x - x_min_)/dx;
	    Real_wp iy_ = (y - y_min_)/dy;
	    Real_wp iz_ = (z - z_min_)/dz;
	    int ix = toInteger(round(ix_));
	    int iy = toInteger(round(iy_));
	    int iz = toInteger(round(iz_));
	    if (abs(ix_ - ix) > 1e-12)
	      ix = toInteger(floor(ix_));
	    
	    if (abs(iy_ - iy) > 1e-12)
	      iy = toInteger(floor(iy_));
	    
	    if (abs(iz_ - iz) > 1e-12)
	      iz = toInteger(floor(iz_));

	    if ((ix < 0) || (ix >= nx))
	      return false;

	    if ((iy < 0) || (iy >= ny))
	      return false;

	    if ((iz < 0) || (iz >= nz))
	      return false;
	    
	    // then we translate the point to the central cell
	    pt3D(0) = this->xmin0 + (x - x_min_ - dx*ix);
	    pt3D(1) = this->ymin0 + (y - y_min_ - dy*iy);
            pt3D(2) = this->zmin0 + (z - z_min_ - dz*iz);
	    
	    num_partie = nz*(ix*ny + iy) + iz;
          }
      }

    // we test if the point is inside the bounding box of the grid
    if ((pt3D(0) >= grid.GetXmin()-R3::threshold) && (pt3D(0) <= grid.GetXmax()+R3::threshold) 
	&& (pt3D(1) >= grid.GetYmin()-R3::threshold) && (pt3D(1) <= grid.GetYmax()+R3::threshold)
	&& (pt3D(2) >= grid.GetZmin()-R3::threshold) && (pt3D(2) <= grid.GetZmax()+R3::threshold))
      return true;
    
    return false;
  }

  
  //! writes parameters of the grid in binary format
  void GridInterpolationFull<Dimension3>::
  Write(ostream& FileStream, bool double_prec) const
  {
    if (double_prec)
      {
	double a;
	WriteBinary(FileStream, a);
      }
    else
      {
	float a;
	WriteBinary(FileStream, a);
      }
  }
  
  
  //! writes parameters of the grid in binary format
  template<class real>
  void GridInterpolationFull<Dimension3>::
  WriteBinary(ostream& FileStream, real& f) const
  {
    real tmp; int nb;
    
    switch (type_output_file)
      {
      case VOLUME :
	tmp = toDouble(xmin); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(xmax); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	nb = nbPoints_x; FileStream.write(reinterpret_cast<char*>(&nb), sizeof(int));
	tmp = toDouble(ymin); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(ymax); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	nb = nbPoints_y; FileStream.write(reinterpret_cast<char*>(&nb), sizeof(int));
	tmp = toDouble(zmin); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(zmax); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	nb = nbPoints_z; FileStream.write(reinterpret_cast<char*>(&nb), sizeof(int));
	break;
      case THREE_PLANES :
	tmp = toDouble(xmin); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(xmax); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	nb = nbPoints_x; FileStream.write(reinterpret_cast<char*>(&nb), sizeof(int));
	tmp = toDouble(ymin); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(ymax); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	nb = nbPoints_y; FileStream.write(reinterpret_cast<char*>(&nb), sizeof(int));
	tmp = toDouble(zmin); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(zmax); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	nb = nbPoints_z; FileStream.write(reinterpret_cast<char*>(&nb), sizeof(int));
        
	tmp = toDouble(center(0)); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(center(1)); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(center(2)); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	break;
      case PLANE :
	// plane OAB with O = (xmin, ymin, zmin)
	// A = (xmax, ymax, zmax), B = center
	tmp = toDouble(xmin); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(xmax); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(ymin); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(ymax); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(zmin); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(zmax); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));

	tmp = toDouble(center(0)); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(center(1)); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(center(2)); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	nb = nbPoints_x; FileStream.write(reinterpret_cast<char*>(&nb), sizeof(int));
	nb = nbPoints_y; FileStream.write(reinterpret_cast<char*>(&nb), sizeof(int));
	break;
      case LINE :
	// line joining A=(xmin,ymin,zmin) and B = (xmax, ymax, zmax)
	tmp = toDouble(xmin); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(xmax); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(ymin); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(ymax); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(zmin); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(zmax); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	nb = nbPoints_x; FileStream.write(reinterpret_cast<char*>(&nb), sizeof(int));
	break;
      case POINT :
	tmp = toDouble(xmin); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(ymin); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(zmin); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	break;
      case CIRCLE :
	tmp = toDouble(center(0)); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(center(1)); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(center(2)); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(xmin); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(ymin); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(zmin); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(xmax); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	tmp = toDouble(ymax); FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	FileStream.write(reinterpret_cast<char*>(const_cast<int*>(&nbPoints_x)),sizeof(int));
	break;
      case POINTS_FILE :
	{
	  nb = GlobalPoints3D.GetM(); FileStream.write(reinterpret_cast<char*>(&nb), sizeof(int));
	  for (int i = 0; i < nb; i++)
	    {
	      tmp = toDouble(GlobalPoints3D(i)(0)); 
              FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	      tmp = toDouble(GlobalPoints3D(i)(1)); 
              FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
              tmp = toDouble(GlobalPoints3D(i)(2)); 
              FileStream.write(reinterpret_cast<char*>(&tmp), sizeof(real));
	    }
	  
	}
	break;
      }
  }
  
  
  //! writes parameters of the grid in text format
  void GridInterpolationFull<Dimension3>::
  WriteText(ostream& FileStream, int type_precision) const
  {
    switch (type_precision)
      {
      case OutputTypeEnum::SINGLE_PRECISION : FileStream.precision(7); break;
      case OutputTypeEnum::DOUBLE_PRECISION : FileStream.precision(15); break;
      case OutputTypeEnum::QUADRUPLE_PRECISION : FileStream.precision(32); break;
      case OutputTypeEnum::MULTIPLE_PRECISION :
        {
          int n = toInteger(-log10(epsilon_machine));
          FileStream.precision(n);
        }
        break;
      }
    
    switch (type_output_file)
      {
      case VOLUME :
	FileStream << xmin << endl << xmax << endl << nbPoints_x << endl;
	FileStream << ymin << endl << ymax << endl << nbPoints_y << endl;
	FileStream << zmin << endl << zmax << endl << nbPoints_z << endl;
	break;
      case THREE_PLANES :
	FileStream << xmin << endl << xmax << endl << nbPoints_x << endl;
	FileStream << ymin << endl << ymax << endl << nbPoints_y << endl;
	FileStream << zmin << endl << zmax << endl << nbPoints_z << endl;
	FileStream << center(0) << endl << center(1) << endl << center(2) << endl;
	break;
      case PLANE :
	// plane OAB with O = (xmin, ymin, zmin)
	// A = (xmax, ymax, zmax), B = center
	FileStream << xmin << endl << xmax << endl;
	FileStream << ymin << endl << ymax << endl;
	FileStream << zmin << endl << zmax << endl;
	FileStream << center(0) << endl << center(1) << endl << center(2) << endl;
	FileStream << nbPoints_x << endl << nbPoints_y << endl;
	break;
      case LINE :
	// line joining A=(xmin,ymin,zmin) and B = (xmax, ymax, zmax)
	FileStream << xmin << endl << xmax << endl;
	FileStream << ymin << endl << ymax << endl;
	FileStream << zmin << endl << zmax << endl;
	FileStream << nbPoints_x << endl;
	break;
      case POINT :
	FileStream << xmin << endl;
	FileStream << ymin << endl;
	FileStream << zmin << endl;
	break;
      case POINTS_FILE :
	{
	  FileStream << GlobalPoints3D.GetM() << endl;
	  for (int i = 0; i < GlobalPoints3D.GetM(); i++)
	    {
	      PrintNoBrackets(FileStream, GlobalPoints3D(i));
	      FileStream << '\n';
	    }	  
	}
	break;
      case CIRCLE :
	FileStream << center(0) << endl << center(1)<<endl
		   << center(2) << endl << xmin << endl << ymin << endl
		   << zmin << endl << xmax << endl << ymax << endl << nbPoints_x << endl;
	break;
      }
	
  }


  //! writes parameters of the interpolation grid in vtk format
  void GridInterpolationFull<Dimension3>::
  WriteVtk(ostream& FileStream, bool double_prec, bool ascii) const
  {
    GridInterpolationFull<Dimension2> grid;
    switch (this->type_output_file)
      {
      case VOLUME :
	grid.SetVolumeAxiOutput(R3(xmin, ymin, zmin), R3(xmax, ymax, zmax),
				nbPoints_x, nbPoints_y, nbPoints_z);
	break;
      case THREE_PLANES :
	grid.SetThreePlanesAxiOutput(center, R3(xmin, ymin, zmin), R3(xmax, ymax, zmax),
				     nbPoints_x, nbPoints_y, nbPoints_z);
	break;
      case PLANE :
	grid.SetPlaneAxiOutput(R3(xmin, ymin, zmin), R3(xmax, ymax, zmax),
			       center, nbPoints_x, nbPoints_y);
	break;
      case LINE :
	grid.SetLineAxiOutput(R3(xmin, ymin, zmin), R3(xmax, ymax, zmax), nbPoints_x);
	break;
      case POINT :
	grid.SetPointAxiOutput(R3(xmin, ymin, zmin));
	break;
      case CIRCLE :
	grid.SetCircleAxiOutput(center, R3(xmin, ymin, zmin), xmax, ymax, nbPoints_x);
	break;
      case POINTS_FILE :
	grid.SetPointsFileAxiOutput(PointsInputFile);
	grid.GlobalPoints3D = GlobalPoints3D;
	grid.SetNbPointsX(GlobalPoints3D.GetM());
	break;
      }
    
    grid.WriteVtk(FileStream, double_prec, ascii);
  }
  
  
  //! reads the display grid from a binary file
  void GridInterpolationFull<Dimension3>::Read(istream& FileStream, bool double_prec)
  {
    if (double_prec)
      {
	double a;
	ReadBinary(FileStream, a);
      }
    else
      {
	float a;
	ReadBinary(FileStream, a);
      }
  }


  //! reads the display grid from a binary file  
  template<class real>
  void GridInterpolationFull<Dimension3>::ReadBinary(istream& FileStream, real& f)
  {
    real tmp; int nb;
    nbPoints_x = 1; nbPoints_y = 1; nbPoints_z = 1;
    switch (type_output_file)
      {
      case VOLUME :
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); xmin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); xmax = tmp;
	FileStream.read(reinterpret_cast<char*>(&nb), sizeof(int)); nbPoints_x = nb;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); ymin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); ymax = tmp;
	FileStream.read(reinterpret_cast<char*>(&nb), sizeof(int)); nbPoints_y = nb;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); zmin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); zmax = tmp;
	FileStream.read(reinterpret_cast<char*>(&nb), sizeof(int)); nbPoints_z = nb;
	break;
      case THREE_PLANES :
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); xmin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); xmax = tmp;
	FileStream.read(reinterpret_cast<char*>(&nb), sizeof(int)); nbPoints_x = nb;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); ymin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); ymax = tmp;
	FileStream.read(reinterpret_cast<char*>(&nb), sizeof(int)); nbPoints_y = nb;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); zmin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); zmax = tmp;
	FileStream.read(reinterpret_cast<char*>(&nb), sizeof(int)); nbPoints_z = nb;

	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); center(0) = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); center(1) = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); center(2) = tmp;
	break;
      case PLANE :
	// plane OAB with O = (xmin, ymin, zmin)
	// A = (xmax, ymax, zmax), B = center
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); xmin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); xmax = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); ymin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); ymax = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); zmin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); zmax = tmp;

	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); center(0) = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); center(1) = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); center(2) = tmp;
	FileStream.read(reinterpret_cast<char*>(&nb), sizeof(int)); nbPoints_x = nb;
	FileStream.read(reinterpret_cast<char*>(&nb), sizeof(int)); nbPoints_y = nb;
	break;
      case LINE :
	// line joining A=(xmin,ymin,zmin) and B = (xmax, ymax, zmax)
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); xmin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); xmax = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); ymin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); ymax = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); zmin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); zmax = tmp;
	FileStream.read(reinterpret_cast<char*>(&nb), sizeof(int)); nbPoints_x = nb;
	break;
      case POINT :
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); xmin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); ymin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); zmin = tmp;
	break;
      case CIRCLE :
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); center(0) = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); center(1) = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); center(2) = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); xmin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); ymin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); zmin = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); xmax = tmp;
	FileStream.read(reinterpret_cast<char*>(&tmp), sizeof(real)); ymax = tmp;
	FileStream.read(reinterpret_cast<char*>(&nb),sizeof(int)); nbPoints_x = nb;
	break;
      case POINTS_FILE :
	{
	  FileStream.read(reinterpret_cast<char*>(&nb),sizeof(int));
	  Vector<real> Pts3D(3*nb);
	  FileStream.read(reinterpret_cast<char*>(Pts3D.GetData()), 3*nb*sizeof(real));
	  nbPoints_x = nb;
	  GlobalPoints3D.Reallocate(nb);
	  for (int i = 0; i < nb; i++)
	    GlobalPoints3D(i).Init(Pts3D(3*i), Pts3D(3*i+1), Pts3D(3*i+2));
	  
	}
	break;
      }
  }
  
  
  //! Reads grid parameters from a text file
  void GridInterpolationFull<Dimension3>::ReadText(istream& FileStream)
  {
    int nb;
    switch (type_output_file)
      {
      case VOLUME :
	FileStream >> xmin >> xmax >> nbPoints_x;
	FileStream >> ymin >> ymax >> nbPoints_y;
	FileStream >> zmin >> zmax >> nbPoints_z;
	break;
      case THREE_PLANES :
	FileStream >> xmin >> xmax >> nbPoints_x;
	FileStream >> ymin >> ymax >> nbPoints_y;
	FileStream >> zmin >> zmax >> nbPoints_z;
	FileStream >> center(0) >> center(1) >> center(2);
	break;
      case PLANE :
	// plane OAB with O = (xmin, ymin, zmin)
	// A = (xmax, ymax, zmax), B = center
	FileStream >> xmin >> xmax >> ymin >> ymax >> zmin >> zmax;
	FileStream >> center(0) >> center(1) >> center(2) >> nbPoints_x >> nbPoints_y;
	break;
      case LINE :
	// line joining A=(xmin,ymin,zmin) and B = (xmax, ymax, zmax)
	FileStream >> xmin >> xmax >> ymin >> ymax >> zmin >> zmax;
	FileStream >> nbPoints_x;
	break;
      case POINT :
	FileStream >> xmin >> ymin >> zmin;
	break;
      case CIRCLE :
	FileStream >> center(0) >> center(1) >> center(2) 
		   >> xmin >> ymin >> zmin >> xmax >> ymax >> nbPoints_x;
	break;
      case POINTS_FILE :
	{
	  FileStream >> nb;
	  nbPoints_x = nb;
	  GlobalPoints3D.Reallocate(nb);
	  for (int i = 0; i < nb; i++)
	    FileStream >> GlobalPoints3D(i)(0) >> GlobalPoints3D(i)(1) >> GlobalPoints3D(i)(2);
	}
	break;
      }

  }


  //! Reads grid parameters from a vtk file
  void GridInterpolationFull<Dimension3>::ReadVtk(istream& FileStream, bool& ascii, int& N)
  {
    GridInterpolationFull<Dimension2> grid;
    grid.ReadVtk(FileStream, ascii, N);

    xmin = grid.GetXmin(); xmax = grid.GetXmax();
    ymin = grid.GetYmin(); ymax = grid.GetYmax();
    zmin = grid.GetZmin(); zmax = grid.GetZmax();
    
    nbPoints_x = grid.GetNbPointsX();
    nbPoints_y = grid.GetNbPointsY();
    nbPoints_z = grid.GetNbPointsZ();
    
    center = grid.GetCenter();
    
    GlobalPoints3D = grid.GlobalPoints3D;
    
    switch (grid.GetOutputType())
      {
      case GridInterpolationFull<Dimension2>::VOLUME_AXI :
	this->type_output_file = VOLUME;
	break;
      case GridInterpolationFull<Dimension2>::THREE_PLANES_AXI :
	this->type_output_file = THREE_PLANES;
	break;
      case GridInterpolationFull<Dimension2>::PLANE_AXI :
	this->type_output_file = PLANE;
	break;
      case GridInterpolationFull<Dimension2>::LINE_AXI :
	this->type_output_file = LINE;
	break;
      case GridInterpolationFull<Dimension2>::POINT_AXI :
	this->type_output_file = POINT;
	break;
      case GridInterpolationFull<Dimension2>::CIRCLE_AXI :
	this->type_output_file = CIRCLE;
	break;
      case GridInterpolationFull<Dimension2>::POINTS_FILE_AXI :
	this->type_output_file = POINTS_FILE;
	break;
      }
  }
  
  
  //! displays details of class GridInterpolationFull<Dimension3>
  template<class Dimension>
  ostream& operator<<(ostream& out, const GridInterpolationFull<Dimension3>& g)
  {
    out<<static_cast<const GridInterpolationFull_Base<Dimension3>& >(g);
    return out;
  }
#endif

  
  //! Writes vtk file with the parameters of the grid and datas
  /*!
    \param[in] val values on the points of the grid that will be written in the vtk file
    \param[in] grid interpolation grid
    \param[in] file_name name of the vtk file
    \param[in] double_prec if true datas are written in double precision,
               otherwise in single precision
    \param[in] ascii if false the vtk file is written in binary format
   */
  template<class T, class Dimension>
  void WriteVtk(const Vector<T>& val, const GridInterpolationFull<Dimension>& grid,
		const string& file_name, bool double_prec, bool ascii)
  {
    ofstream file_out(file_name.data());
    if (!file_out.is_open())
      {
	cout << "Failed to open file " << file_name << endl;
	abort();
      }

    // writing header
    grid.WriteVtk(file_out, double_prec, ascii);
    
    // then datas
    WriteVtk(val, "u", file_out, double_prec, ascii);
    file_out.close();
  }
  
  
  //! reads the interpolation grid and data from a vtk file
  /*!
    \param[out] val values on the points of the grid
    \param[out] grid interpolation grid
    \param[in] file_name name of the vtk file
   */
  template<class T, class Dimension>
  void ReadVtk(Vector<T>& val, GridInterpolationFull<Dimension>& grid,
	       const string& file_name, bool double_prec, bool ascii)
  {
    ifstream file_in(file_name.data());
    if (!file_in.is_open())
      {
	cout << "Failed to open file " << file_name << endl;
	abort();
      }

    // reading header
    int N;
    grid.ReadVtk(file_in, ascii, N);
    
    // then datas
    string name; val.Reallocate(N);
    ReadVtk(val, name, file_in, ascii);
    file_in.close();
  }


  void WriteBinaryTypeData_DoubleOrFloat(const VectReal_wp& output_vector,
					 ostream& file_out, bool double_prec)
  {
    int itmp = 0;
    if (double_prec)
      itmp = 1;
    
    file_out.write(reinterpret_cast<char*>(const_cast<int*>(&itmp)), sizeof(int));
  }

  
  void WriteBinaryTypeData_DoubleOrFloat(const VectComplex_wp& output_vector,
					 ostream& file_out, bool double_prec)
  {
    int itmp = 2;
    if (double_prec)
      itmp = 3;
    
    file_out.write(reinterpret_cast<char*>(const_cast<int*>(&itmp)), sizeof(int));
  }
  
  
  void WriteTypeData_DoubleOrFloat(const VectReal_wp& output_vector,
				   ostream& file_out, int type_prec)
  {
    if (type_prec == OutputTypeEnum::DOUBLE_PRECISION)
      file_out<<"1"<<endl;
    else if (type_prec == OutputTypeEnum::QUADRUPLE_PRECISION)
      file_out<<"4"<<endl;
    else if (type_prec == OutputTypeEnum::MULTIPLE_PRECISION)
      file_out<<"5"<<endl;
    else
      file_out<<"0"<<endl;
  }
  
  
  void WriteTypeData_DoubleOrFloat(const VectComplex_wp& output_vector,
				    ostream& file_out, int type_prec)
  {
    if (type_prec == OutputTypeEnum::DOUBLE_PRECISION)
      file_out<<"3"<<endl;
    else if (type_prec == OutputTypeEnum::QUADRUPLE_PRECISION)
      file_out<<"6"<<endl;
    else if (type_prec == OutputTypeEnum::MULTIPLE_PRECISION)
      file_out<<"7"<<endl;
    else
      file_out<<"2"<<endl;
  }
  
  
  //! writes datas and interpolation in Matlab format (.dat)
  /*!
    \param[in] val values of the data on points of the interpolation grid
    \param[in] grid interpolation grid
    \param[in] file_name name of the Matlab file
    \param[in] double_prec if true, the datas will be written in double precision
    \param[in] ascii if false, the datas will be written in binary format
    The file can be read in Matlab with the function loadND
   */
  template<class T, class Dimension>
  void WriteMatlab(const Vector<T>& val,
		   GridInterpolationFull<Dimension>& grid,
		   const string& file_name, int double_prec, bool ascii)
  {
    Vector<GridInterpolationFull<Dimension> > var_grid;
    var_grid.SetData(1, &grid);
    
    WriteMatlab(val, 1, var_grid, grid.GetOutputType(),
		file_name, double_prec, ascii);
    
    var_grid.Nullify();
  }	  


  //! writes datas and interpolation in Matlab format (.dat)
  /*!
    \param[in] output_vector values of the data on points of the interpolation grids
    \param[in] nb_true the number of interpolation grids
    \param[in] grid vector of interpolation grids
    \param[in] type the type of grid to take into account
    \param[in] file_name name of the Matlab file
    \param[in] double_prec if true, the datas will be written in double precision
    \param[in] ascii if false, the datas will be written in binary format
    The file can be read in Matlab with the function loadND
   */
  template<class T, class Dimension>
  void WriteMatlab(const Vector<T>& output_vector, int nb_true,
		   const Vector<GridInterpolationFull<Dimension> >& grid, int type,
		   const string& file_name, int type_precision, bool ascii)
  {
    ofstream file_out;
    
    file_out.open(file_name.data());
    if (!file_out.is_open())
      {
	cout << "Failed to open file " << file_name << endl;
	abort();
      }
    
    if (!ascii)
      {	
        bool double_prec = false;
        if (type_precision == OutputTypeEnum::DOUBLE_PRECISION)
          double_prec = true;
        
	// binary file for Matlab/Python
        
	// first field : number of grids in the file
	file_out.write(reinterpret_cast<char*>(const_cast<int*>(&nb_true)),
		       sizeof(int));
	
	// second field : dimension
	nb_true = GridInterpolationFull<Dimension>::GetDimension(type);
	file_out.write(reinterpret_cast<char*>(const_cast<int*>(&nb_true)),
		       sizeof(int));
	
	// third field : type of data (float, double, complex<float> or complex<double>)
	WriteBinaryTypeData_DoubleOrFloat(output_vector, file_out, double_prec);
	    
	// fourth field : type of grid
	nb_true = ParamOutputClass::GetTrueType(type, Dimension::dim_N);
	file_out.write(reinterpret_cast<char*>(const_cast<int*>(&nb_true)),
		       sizeof(int));
	
	// then for each grid, we write parameters of the grid
	for (int i = 0; i < grid.GetM(); i++)
	  if (grid(i).GetOutputType() == type)
	    grid(i).Write(file_out, double_prec);
	
	// then we write all the values (for all the grids)
	WriteBinaryDoubleOrFloat(output_vector, file_out, double_prec);
      }
    else
      {
	// ascii file for Matlab/Python
	
	// first field : number of grids
	file_out << nb_true << endl;
        
	// second field : dimension
	file_out << GridInterpolationFull<Dimension>::GetDimension(type) << endl;
        
	// third field : type of data
	WriteTypeData_DoubleOrFloat(output_vector, file_out, type_precision);
        
	// fourth field : type of grid
	file_out << ParamOutputClass::GetTrueType(type, Dimension::dim_N) << endl;
        
	// then for each grid, we write parameters of the grid
	for (int i = 0; i < grid.GetM(); i++)
	  if (grid(i).GetOutputType() == type)
	    grid(i).WriteText(file_out, type_precision);
	
	// then all the values
	for (int i = 0; i < output_vector.GetM(); i++)
	  file_out << output_vector(i) << '\n';
	
      }
    
    file_out.close();
  }  

  
  //! reads datas and interpolation grid from Matlab file (.dat)
  /*!
    \param[out] val values of the data on points of the interpolation grid
    \param[out] grid interpolation grid
    \param[in] file_name name of the Matlab file
    \param[in] ascii if false, the datas are stored in binary format
   */
  template<class T, class Dimension>
  void ReadMatlab(Vector<T>& val,
		  GridInterpolationFull<Dimension>& grid,
		  const string& file_name, bool ascii)
  {
    Vector<GridInterpolationFull<Dimension> > var_grid;
    
    ReadMatlab(val, var_grid, file_name, ascii);
    grid = var_grid(0);
  }	  
  
  
  //! reads datas and interpolation grids from Matlab file (.dat)
  /*!
    \param[out] output_vector values of the data on points of the interpolation grids
    \param[out] grid interpolation grids
    \param[in] file_name name of the Matlab file
    \param[in] ascii if false, the datas are stored in binary format
   */
  template<class T, class Dimension>
  void ReadMatlab(Vector<T>& output_vector,
		  Vector<GridInterpolationFull<Dimension> >& grid,
		  const string& file_name, bool ascii)
  {
    ifstream file_in;
    
    file_in.open(file_name.data());
    if (!file_in.is_open())
      {
	cout << "Failed to open file " << file_name << endl;
	abort();
      }
    
    if (!ascii)
      {	
	// binary file for Matlab/Python

	// first field : number of grids in the file
	int nb_grid;
	file_in.read(reinterpret_cast<char*>(&nb_grid), sizeof(int));
	grid.Reallocate(nb_grid);
	
	// second field : dimension
	int dim;
	file_in.read(reinterpret_cast<char*>(&dim), sizeof(int));
	
	// third field : type of data (float, double, complex<float> or complex<double>)
	int type_data;
	file_in.read(reinterpret_cast<char*>(&type_data), sizeof(int));
	bool double_prec = false;
	if (type_data%2 == 1)
	  double_prec = true;
		    
	// fourth field : type of grid
	int type_grid;
	file_in.read(reinterpret_cast<char*>(&type_grid), sizeof(int));
			
	int type_output = ParamOutputClass::GetFalseType(type_grid, dim, Dimension::dim_N);
	// then for each grid, we read parameters of the grid
	for (int i = 0; i < grid.GetM(); i++)
	  {
	    grid(i).SetOutputType(type_output);
	    grid(i).Read(file_in, double_prec);
	  }
	
	// then we read all the values (for all the grids)
	if (type_data == 0)
	  {
	    Vector<float> val;
	    val.Read(file_in);
	    CopyVector(val, output_vector);
	  }
	else if (type_data == 1)
	  {
	    Vector<double> val;
	    val.Read(file_in);
	    CopyVector(val, output_vector);
	  }
	else if (type_data == 2)
	  {
	    Vector<complex<float> > val;
	    val.Read(file_in);
	    CopyVector(val, output_vector);
	  }
	else if (type_data == 3)
	  {
	    Vector<complex<double> > val;
	    val.Read(file_in);
	    CopyVector(val, output_vector);
	  }
      }
    else
      {
	// ascii file for Matlab/Python
	
	// first field : number of grids
	int nb_grid;
	file_in >> nb_grid;
	grid.Reallocate(nb_grid);
        
	// second field : dimension
	int dim;
	file_in >> dim;
        
	// third field : type of data
	int type_data;
	file_in >> type_data;
	
	// fourth field : type of grid
	int type_grid;
	file_in >> type_grid;
        
	// then for each grid, we read parameters of the grid
	int type_output = ParamOutputClass::GetFalseType(type_grid, dim, Dimension::dim_N);
	for (int i = 0; i < grid.GetM(); i++)
	  {
	    grid(i).SetOutputType(type_output);
	    grid(i).ReadText(file_in);
	  }
	
	// then all the values
	output_vector.ReadText(file_in);	
      }
    
    file_in.close();
  }  
  
}

#define MONTJOIE_FILE_GRID_INTERPOLATION_CXX
#endif
