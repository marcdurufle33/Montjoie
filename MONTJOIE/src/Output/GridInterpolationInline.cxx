#ifndef MONTJOIE_FILE_GRID_INTERPOLATION_INLINE_CXX

namespace Montjoie
{
  /**************************
   * GridInterpolation_Base *
   **************************/
  
  
  //! returns the number of points stored
  template<class Dimension>
  inline int GridInterpolation_Base<Dimension>::GetNbElt() const
  {
    return ElementInterp.GetM();
  }


  template<class Dimension>
  inline int GridInterpolation_Base<Dimension>::GetNbGlobalPoints() const
  {
    return GlobalCoord.GetM();
  }

  
  //! returns the element where the point i is
  template<class Dimension>
  inline int GridInterpolation_Base<Dimension>::GetElementNumber(int i) const
  {
    return ElementInterp(i);
  }
  
  
  //! sets the element where the point is
  template<class Dimension>
  inline void GridInterpolation_Base<Dimension>::SetElementNumber(int i, int n)
  {
    ElementInterp(i) = n;
  }
  
  
  //! returns local coordinates of point i on the element
  template<class Dimension>
  inline const typename Dimension::R_N& GridInterpolation_Base<Dimension>::
  GetLocalCoordinate(int i) const
  {
    return CoorInterp(i);
  }

  
  //! return global coordinates of point i
  template<class Dimension>
  inline const typename Dimension::R_N& GridInterpolation_Base<Dimension>::
  GetGlobalCoordinate(int i) const
  {
    return GlobalCoord(i);
  }

  
  //! return global coordinates of all points
  template<class Dimension>
  inline const typename Dimension::VectR_N& GridInterpolation_Base<Dimension>::
  GetGlobalCoordinate() const
  {
    return GlobalCoord;
  }


  //! return global coordinates of all points
  template<class Dimension>
  inline typename Dimension::VectR_N& GridInterpolation_Base<Dimension>::GetGlobalCoordinate()
  {
    return GlobalCoord;
  }

  
  //! sets local coordinates of point i
  template<class Dimension>
  inline void GridInterpolation_Base<Dimension>::SetLocalCoordinate(int i, const R_N& x)
  {
    CoorInterp(i) = x;
  }
  
  
  //! sets global coordinates of point i
  template<class Dimension>
  inline void GridInterpolation_Base<Dimension>::SetGlobalCoordinate(int i, const R_N& x)
  {
    GlobalCoord(i) = x;
  }
  
  
  //! returns the number of theta-parameters stored
  template<class Dimension>
  inline int GridInterpolation_Base<Dimension>::GetNbTheta() const
  {
    return TetaInterp.GetM();
  }
  

  //! returns theta-parameter of point i (for axisymmetric computations)
  template<class Dimension>
  inline const Real_wp& GridInterpolation_Base<Dimension>::GetTheta(int i) const
  {
    return TetaInterp(i);
  }


  //! returns theta-parameter of point i (for axisymmetric computations)
  template<class Dimension>
  inline const VectReal_wp& GridInterpolation_Base<Dimension>::GetTheta() const
  {
    return TetaInterp;
  }


  //! returns theta-parameter of point i (for axisymmetric computations)
  template<class Dimension>
  inline VectReal_wp& GridInterpolation_Base<Dimension>::GetTheta()
  {
    return TetaInterp;
  }

  
  //! Sets theta parameter of point i
  template<class Dimension>
  inline void GridInterpolation_Base<Dimension>::SetTheta(int i, const Real_wp& teta)
  {
    TetaInterp(i) = teta;
  }
  
  
  template<class Dimension>
  inline bool GridInterpolation_Base<Dimension>::StoreDFjm1() const
  {
    return dfjm1_store;
  }
  
  
  //! returns inverse of jacobian matrix DF_i^{-1} on point i
  template<class Dimension>
  inline const typename Dimension::MatrixN_N&
  GridInterpolation_Base<Dimension>::GetDFjm1(int i) const
  {
    return dfjm1(i);
  }
  

  template<class Dimension>
  inline void GridInterpolation_Base<Dimension>::SetDFjm1(int i, const MatrixN_N& A)
  {
    dfjm1(i) = A;
  }


  template<class Dimension>
  inline int GridInterpolation_Base<Dimension>::GetNbSectionNumber() const
  {
    return CyclicSection.GetM();
  }
  

  template<class Dimension>
  inline int GridInterpolation_Base<Dimension>::GetSectionNumber(int i) const
  {
    return CyclicSection(i);
  }

  
  template<class Dimension>
  inline void GridInterpolation_Base<Dimension>::SetSectionNumber(int i, int n)
  {
    CyclicSection(i) = n;
  }


  template<class Dimension>
  inline Real_wp GridInterpolation_Base<Dimension>::GetEpsilonInsideElement(int i) const
  {
    return Epsilon_InsideElement(i);
  }
  
  
  template<class Dimension>
  inline void GridInterpolation_Base<Dimension>::SetEpsilonInsideElement(int i, const Real_wp& eps)
  {
    Epsilon_InsideElement(i) = eps;
  }
  
  
  template<class Dimension>
  inline const Real_wp& GridInterpolation_Base<Dimension>::GetXmin() const
  {
    return xmin;
  }
  
  
  template<class Dimension>
  inline const Real_wp& GridInterpolation_Base<Dimension>::GetXmax() const
  {
    return xmax;
  }
  
  
  template<class Dimension>
  inline const Real_wp& GridInterpolation_Base<Dimension>::GetYmin() const
  {
    return ymin;
  }
  
  
  template<class Dimension>
  inline const Real_wp& GridInterpolation_Base<Dimension>::GetYmax() const
  {
    return ymax;
  }
  
  
  template<class Dimension>
  inline const Real_wp& GridInterpolation_Base<Dimension>::GetZmin() const
  {
    return zmin;
  }
  
  
  template<class Dimension>
  inline const Real_wp& GridInterpolation_Base<Dimension>::GetZmax() const
  {
    return zmax;
  }

  
  template<class Dimension>
  inline void GridInterpolation_Base<Dimension>::SetXmin(const Real_wp& x)
  {
    xmin = x;
  }
  
  
  template<class Dimension>
  inline void GridInterpolation_Base<Dimension>::SetXmax(const Real_wp& x)
  {
    xmax = x;
  }
  
  
  template<class Dimension>
  inline void GridInterpolation_Base<Dimension>::SetYmin(const Real_wp& y)
  {
    ymin = y;
  }
  
  
  template<class Dimension>
  inline void GridInterpolation_Base<Dimension>::SetYmax(const Real_wp& y)
  {
    ymax = y;
  }
  
  
  template<class Dimension>
  inline void GridInterpolation_Base<Dimension>::SetZmin(const Real_wp& z)
  {
    zmin = z;
  }
  
  
  template<class Dimension>
  inline void GridInterpolation_Base<Dimension>::SetZmax(const Real_wp& z)
  {
    zmax = z;
  }

  
  template<class Dimension>
  inline void GridInterpolation_Base<Dimension>::SetNbSubdivisions(int nx, int ny)
  {
    nb_subdiv_grid_x = nx;
    nb_subdiv_grid_y = ny;
  }
  
  
  template<class Dimension>
  inline void GridInterpolation_Base<Dimension>::SetSubdivisionStep(const Real_wp& dx, const Real_wp& dy)
  {
    step_subdiv_grid_x = dx;
    step_subdiv_grid_y = dy;
  }

  
  template<class Dimension>
  inline void GridInterpolation_Base<Dimension>::SetNbSubdivisions(int nx, int ny, int nz)
  {
    nb_subdiv_grid_x = nx; nb_subdiv_grid_y = ny; nb_subdiv_grid_z = nz; 
  }
  
  
  template<class Dimension>
  inline void GridInterpolation_Base<Dimension>
  ::SetSubdivisionStep(const Real_wp& dx, const Real_wp& dy, const Real_wp& dz)
  { 
    step_subdiv_grid_x = dx; step_subdiv_grid_y = dy; step_subdiv_grid_z = dz;
  }


#ifdef SELDON_WITH_MPI
  template<class Dimension>
  inline void GridInterpolation_Base<Dimension>::SetCommunicator(const MPI_Comm& comm)
  {
    comm_ = comm;
  }


  template<class Dimension>
  inline const MPI_Comm& GridInterpolation_Base<Dimension>::GetCommunicator() const
  {
    return comm_;
  }
#endif
  

  /******************************
   * GridInterpolationFull_Base *
   ******************************/

  //! returns true if the grid intervals have to be retrieved from the mesh
  template<class Dimension>
  inline bool GridInterpolationFull_Base<Dimension>::GridIntervalToBeComputed() const
  {
    return grid_interval_to_be_computed;
  }
  
  
  //! returns the type of grid
  template<class Dimension>
  inline int GridInterpolationFull_Base<Dimension>::GetOutputType() const
  {
    return type_output_file;
  }
  
  
  //! sets the type of grid
  template<class Dimension>
  inline void GridInterpolationFull_Base<Dimension>::SetOutputType(int i)
  {
    type_output_file = i;
  }
        
  
  //! returns the number of stored points in this grid
  template<class Dimension>
  inline int GridInterpolationFull_Base<Dimension>::GetNbPoints() const
  {
    return list_points.GetM();
  }


  //! returns the global number of stored points in this grid
  template<class Dimension>
  inline int GridInterpolationFull_Base<Dimension>::GetNbGlobalPoints() const
  {
    return nb_global_grid_points;
  }


  //! returns the number of point in GridInterpolation object
  template<class Dimension>
  inline int GridInterpolationFull_Base<Dimension>::GetPointNumber(int i) const
  {
    return list_points(i);
  }


  //! returns the point numbers in GridInterpolation object
  template<class Dimension>
  inline const IVect& GridInterpolationFull_Base<Dimension>::GetPointNumber() const
  {
    return list_points;
  }


  //! returns the point numbers for each processor
  template<class Dimension>
  inline const Vector<IVect>& GridInterpolationFull_Base<Dimension>::GetPointNumberAllProc() const
  {
    return list_points_proc;
  }
  
  
  //! returns the number of points along x
  template<class Dimension>
  inline int GridInterpolationFull_Base<Dimension>::GetNbPointsX() const
  {
    return nbPoints_x;
  }
  

  //! returns the number of points along y  
  template<class Dimension>
  inline int GridInterpolationFull_Base<Dimension>::GetNbPointsY() const
  {
    return nbPoints_y;
  }
  
  
  //! returns the number of points along z
  template<class Dimension>
  inline int GridInterpolationFull_Base<Dimension>::GetNbPointsZ() const
  {
    return nbPoints_z;
  }
  
  
  //! sets the number of points along x
  template<class Dimension>
  inline void GridInterpolationFull_Base<Dimension>::SetNbPointsX(int i)
  {
    nbPoints_x = i;
  }
  
  
  //! sets the number of points along y
  template<class Dimension>
  inline void GridInterpolationFull_Base<Dimension>::SetNbPointsY(int i)
  {
    nbPoints_y = i;
  }
  
  
  //! sets the number of points along z
  template<class Dimension>
  inline const R3& GridInterpolationFull_Base<Dimension>::GetCenter() const
  {
    return center;
  }
  
  
  //! returns the minimal value of x
  template<class Dimension>
  inline const Real_wp& GridInterpolationFull_Base<Dimension>::GetXmin() const
  {
    return xmin;
  }
  
  
  //! returns the maximal value of x
  template<class Dimension>
  inline const Real_wp& GridInterpolationFull_Base<Dimension>::GetXmax() const
  {
    return xmax;
  }
  

  //! returns the minimal value of y
  template<class Dimension>
  inline const Real_wp& GridInterpolationFull_Base<Dimension>::GetYmin() const
  {
    return ymin;
  }
  
  
  //! returns the maximal value of y
  template<class Dimension>
  inline const Real_wp& GridInterpolationFull_Base<Dimension>::GetYmax() const
  {
    return ymax;
  }
  
  
  //! returns the minimal value of z
  template<class Dimension>
  inline const Real_wp& GridInterpolationFull_Base<Dimension>::GetZmin() const
  {
    return zmin;
  }
  
  
  //! returns the maximal value of z
  template<class Dimension>
  inline const Real_wp& GridInterpolationFull_Base<Dimension>::GetZmax() const
  {
    return zmax;
  }
    
  
  //! sets the minimal value of x
  template<class Dimension>
  inline void GridInterpolationFull_Base<Dimension>::SetXmin(const Real_wp& x)
  {
    xmin = x;
  }
  
  
  //! sets the maximal value of x
  template<class Dimension>
  inline void GridInterpolationFull_Base<Dimension>::SetXmax(const Real_wp& x)
  {
    xmax = x;
  }
  

  //! sets the minimal value of y
  template<class Dimension>
  inline void GridInterpolationFull_Base<Dimension>::SetYmin(const Real_wp& x)
  {
    ymin = x;
  }
  
  
  //! sets the maximal value of y
  template<class Dimension>
  inline void GridInterpolationFull_Base<Dimension>::SetYmax(const Real_wp& x)
  {
    ymax = x;
  }
  
  
  //! sets the minimal value of z
  template<class Dimension>
  inline void GridInterpolationFull_Base<Dimension>::SetZmin(const Real_wp& x)
  {
    zmin = x;
  }
  
  
  //! sets the maximal value of z
  template<class Dimension>
  inline void GridInterpolationFull_Base<Dimension>::SetZmax(const Real_wp& x)
  {
    zmax = x;
  }


  //! sets the minimal value of x
  template<class Dimension>
  inline void GridInterpolationFull_Base<Dimension>::SetXmin0(const Real_wp& x)
  {
    xmin0 = x;
  }

  
  //! sets the maximal value of x
  template<class Dimension>
  inline void GridInterpolationFull_Base<Dimension>::SetXmax0(const Real_wp& x)
  {
    xmax0 = x;
  }
  

  //! sets the minimal value of y
  template<class Dimension>
  inline void GridInterpolationFull_Base<Dimension>::SetYmin0(const Real_wp& x)
  {
    ymin0 = x;
  }
  
  
  //! sets the maximal value of y
  template<class Dimension>
  inline void GridInterpolationFull_Base<Dimension>::SetYmax0(const Real_wp& x)
  {
    ymax0 = x;
  }
  
  
  //! sets the minimal value of z
  template<class Dimension>
  inline void GridInterpolationFull_Base<Dimension>::SetZmin0(const Real_wp& x)
  {
    zmin0 = x;
  }
  
  
  //! sets the maximal value of z
  template<class Dimension>
  inline void GridInterpolationFull_Base<Dimension>::SetZmax0(const Real_wp& x)
  {
    zmax0 = x;
  }

  
  /*********************************
   * GridInterpolation<Dimension1> *
   *********************************/
  
  
#ifdef MONTJOIE_WITH_ONE_DIM
  inline const Real_wp& GridInterpolation<Dimension1>::GetXmin()
  {
    return Points(0);
  }
  
  
  inline const Real_wp& GridInterpolation<Dimension1>::GetXmax()
  {
    return Points(Points.GetM()-1);
  }
  
  
  inline void GridInterpolation<Dimension1>::AddPoint(const Real_wp& x)
  {
    Points.PushBack(x);
  }

  
  //! returns the number of points
  inline int GridInterpolation<Dimension1>::GetNbPointsGrid() const
  {
    return Points.GetM();
  }
  
  
  //! returns the element number where the point i is
  inline int GridInterpolation<Dimension1>::GetElementNumber(int i) const
  {
    return ElemInterp(i);
  }
  
  
  //! returns the local coordinate of point i in the element
  inline Real_wp GridInterpolation<Dimension1>::GetLocalCoordinate(int i) const
  {
    return CoorInterp(i);
  }
  
  
  //! returns global coordinate of point i
  inline Real_wp GridInterpolation<Dimension1>::GetGlobalCoordinate(int i) const
  {
    return Points(i);
  }
  
  
  //! returns all the points
  inline const VectReal_wp& GridInterpolation<Dimension1>::
  GetGlobalCoordinate() const
  {
    return Points;
  }

  
  //! allocates arrays to contain n points in the interpolation grid
  inline void GridInterpolation<Dimension1>::ReallocatePoints(int N)
  {
    Points.Reallocate(N); Points.Fill(0);
    ElemInterp.Reallocate(N); ElemInterp.Fill(-1);
    CoorInterp.Reallocate(N);    
  }
  
  
  //! sets the element number associated with point i
  inline void GridInterpolation<Dimension1>::SetElementNumber(int i, int elt)
  {
    ElemInterp(i) = elt;
  }
  
  
  //! sets the local coordinate in the element
  inline void GridInterpolation<Dimension1>::SetLocalCoordinate(int i, const Real_wp& x)
  {
    CoorInterp(i) = x;
  }
  
  
  //! sets the global coordinate of point i
  inline void GridInterpolation<Dimension1>::SetGlobalCoordinate(int i, const Real_wp& x)
  {
    Points(i) = x;
  }
  
  
  //! returns the number of points contained in all the processors
  inline int GridInterpolation<Dimension1>::GetGlobalNumberPoints() const
  {
    return nb_global_points;
  }
  
  
  //! sets the number of points contained in all the processors
  inline void GridInterpolation<Dimension1>::SetGlobalNumberPoints(int n)
  {
    nb_global_points = n;
  }
  
  
  /*************************************
   * GridInterpolationFull<Dimension1> *
   *************************************/
  
  
  //! default constructor
  inline GridInterpolationFull<Dimension1>::GridInterpolationFull()
  {
  }
  
  
  //! returns the dimension of the output grid (3 for axisymmetric outputs)
  inline int GridInterpolationFull<Dimension1>::GetDimension(int type)
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
      case PLANE :
	return 2;
      }
    
    return 1;
  }

  
  //! not used
  inline void GridInterpolationFull<Dimension1>::ReadText(istream& FileStream)
  {
    abort();
  }
    
  
  //! not used
  inline void GridInterpolationFull<Dimension1>::Read(istream& FileStream, bool double_prec)
  {
    abort();
  }
  

  //! not used
  inline void GridInterpolationFull<Dimension1>::ReadVtk(istream& FileStream, bool&, int&)
  {
    abort();
  }
  
  
  //! writing header in binary
  inline void GridInterpolationFull<Dimension1>::Write(ostream& FileStream, bool double_prec) const
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

  //! writing header in ascii
  inline void GridInterpolationFull<Dimension1>
  ::WriteText(ostream& FileStream, int type_precision) const
  {
    abort();
  }
  

  //! writing header in ascii
  inline void GridInterpolationFull<Dimension1>
  ::WriteVtk(ostream& FileStream, bool double_prec, bool ascii) const
  {
    abort();
  }

  
  //! not used
  inline void GridInterpolationFull<Dimension1>::WriteGrid(const string& file_name)
  {
    abort();
  } 
    
#endif
  
  
  /*********************************
   * GridInterpolation<Dimension2> *
   *********************************/
  
  
#ifdef MONTJOIE_WITH_TWO_DIM
  //! default constructor
  inline GridInterpolation<Dimension2>::GridInterpolation()
    : GridInterpolation_Base<Dimension2>()
  {
  }
  
  
  /*************************************
   * GridInterpolationFull<Dimension2> *
   *************************************/
  
  
  //! default constructor
  inline GridInterpolationFull<Dimension2>::GridInterpolationFull() 
    : GridInterpolationFull_Base<Dimension2>()
  {
  }
  

  inline void GridInterpolationFull<Dimension2>::SetPointsOutput(const VectR2& ptA)
  {
    GlobalPoints2D = ptA;
    type_output_file = POINTS_FILE;
  } 
  
  
  //! not implemented
  inline void GridInterpolationFull<Dimension2>::WriteGrid(const string& file_name)
  {
    abort();
  }
  
#endif
  
#ifdef MONTJOIE_WITH_THREE_DIM
  /*********************************
   * GridInterpolation<Dimension3> *
   *********************************/
  
  
  //! default constructor
  inline GridInterpolation<Dimension3>::GridInterpolation()
    : GridInterpolation_Base<Dimension3>()
  {
  }
  
    
  /*************************************
   * GridInterpolationFull<Dimension3> *
   *************************************/
  
  
  //! default constructor
  inline GridInterpolationFull<Dimension3>::GridInterpolationFull()
  {
  }
  
  
  inline void GridInterpolationFull<Dimension3>::SetPointsOutput(const VectR3& ptA)
  {
    GlobalPoints3D = ptA;
    type_output_file = POINTS_FILE;
  } 

  
  //! writes the display grid in a file
  inline void GridInterpolationFull<Dimension3>::WriteGrid(const string& file_name)
  {
    abort();
  }
  
#endif
  
}

#define MONTJOIE_FILE_GRID_INTERPOLATION_INLINE_CXX
#endif
