#ifndef MONTJOIE_FILE_PARAMETER_OUTPUT_INLINE_CXX

namespace Montjoie
{
  inline int ParamOutputClass::GetComponent() const
  {
    return component;
  }


  inline int ParamOutputClass::GetNature() const
  {
    return nature;
  }

  
  inline void ParamOutputClass::SetNature(int i)
  {
    nature = i;
  }
  
  
  inline void ParamOutputClass::SetComponent(int i)
  {
    component = i;
  }

  
  inline int ParamOutputClass::IncrementSnapshot() const
  {
    int n = num_instantane;
    num_instantane++;  
    t_current = this->t_begin + num_instantane*this->deltat;
    return n;
  }
  
  
  inline int ParamOutputClass::GetSnapshotNumber() const
  {
    return num_instantane;
  }

  
  inline int ParamOutputClass::GetTrueType(int type, int dim)
  {
    if (dim == 2)
      {
        switch (type)
          {
#ifdef MONTJOIE_WITH_TWO_DIM
          case GridInterpolationFull<Dimension2>::PLANE : return 0;
          case GridInterpolationFull<Dimension2>::PLANE_AXI : return 2;
          case GridInterpolationFull<Dimension2>::LINE : return 1;
          case GridInterpolationFull<Dimension2>::LINE_AXI : return 3;
          case GridInterpolationFull<Dimension2>::POINT : return 2;
          case GridInterpolationFull<Dimension2>::POINT_AXI : return 4;
          case GridInterpolationFull<Dimension2>::CIRCLE : return 3;
          case GridInterpolationFull<Dimension2>::CIRCLE_AXI : return 6;
          case GridInterpolationFull<Dimension2>::THREE_PLANES_AXI : return 1;
          case GridInterpolationFull<Dimension2>::VOLUME_AXI : return 0;
          case GridInterpolationFull<Dimension2>::POINTS_FILE : return 6;
          case GridInterpolationFull<Dimension2>::POINTS_FILE_AXI : return 5;
#endif
          }
      }
    else if (dim == 1)
      {
        switch (type)
          {
#ifdef MONTJOIE_WITH_ONE_DIM
	  case GridInterpolationFull<Dimension1>::PLANE : return 0;
          case GridInterpolationFull<Dimension1>::PLANE_AXI : return 2;
          case GridInterpolationFull<Dimension1>::LINE : return 0;
          case GridInterpolationFull<Dimension1>::LINE_AXI : return 3;
          case GridInterpolationFull<Dimension1>::POINT : return 1;
          case GridInterpolationFull<Dimension1>::POINT_AXI : return 4;
          case GridInterpolationFull<Dimension1>::CIRCLE_AXI : return 6;
          case GridInterpolationFull<Dimension1>::THREE_PLANES_AXI : return 1;
          case GridInterpolationFull<Dimension1>::VOLUME_AXI : return 0;
          case GridInterpolationFull<Dimension1>::POINTS_FILE_AXI : return 5;
#endif
          }
      }
    
    return type;
  }


  inline int ParamOutputClass::GetFalseType(int type, int dim, int dim_N)
  {
    if (dim_N == 2)
      {
#ifdef MONTJOIE_WITH_TWO_DIM
        if (dim == 3)
	  switch (type)
	    {
	    case 0 : return GridInterpolationFull<Dimension2>::VOLUME_AXI;
	    case 1 : return GridInterpolationFull<Dimension2>::THREE_PLANES_AXI;
	    case 2 : return GridInterpolationFull<Dimension2>::PLANE_AXI;
	    case 3 : return GridInterpolationFull<Dimension2>::LINE_AXI;
	    case 4 : return GridInterpolationFull<Dimension2>::POINT_AXI;
	    case 5 : return GridInterpolationFull<Dimension2>::POINTS_FILE_AXI;
	    case 6 : return GridInterpolationFull<Dimension2>::CIRCLE_AXI;
	    }
	else if (dim == 2)
	  switch (type)
	    {
	    case 0 : return GridInterpolationFull<Dimension2>::PLANE;
	    case 1 : return GridInterpolationFull<Dimension2>::LINE;
	    case 2 : return GridInterpolationFull<Dimension2>::POINT;
	    case 3 : return GridInterpolationFull<Dimension2>::CIRCLE;
	    case 6 : return GridInterpolationFull<Dimension2>::POINTS_FILE;
	    }
#endif
      }
    else
      {
	if (dim == 2)
	  {
	    cout << "Case not handled" << endl;
	    abort();
	  }
      }
    
    return type;
  }
  
  
  inline const string& ParamOutputClass::GetTotalFieldFile() const
  {
    return name_file(1);
  }


  inline const string& ParamOutputClass::GetDiffractedFieldFile() const
  {
    return name_file(0);
  }

  
  inline int ParamOutputClass::GetNbFile() const
  {
    return name_file.GetM();
  }
  

  inline void ParamOutputClass::ResizeNbFile(int n)
  {
    name_file.Resize(n);
  }

  
  inline void ParamOutputClass::AdimTime(const Real_wp& t0)
  {
    t_begin /= t0; t_end /= t0;
    deltat /= t0; epsilon_time /= t0;
  }
  

  inline const Real_wp& ParamOutputClass::GetCurrentTime() const
  {
    return t_current;
  }
  
  
  inline void ParamOutputClass::SetThresholdTime(const Real_wp& eps)
  {
    epsilon_time = eps; 
  }

  
  inline void ParamOutputClass::SetTotalFieldFile(const string& name)
  {
    name_file(1) = name;
  }
  

  inline void ParamOutputClass::SetDiffractedFieldFile(const string& name)
  {
    name_file(0) = name;
  }
  

  inline void ParamOutputClass::SetFileName(int i, const string& name)
  {
    if (i >= name_file.GetM())
      name_file.Resize(i+1);
    
    name_file(i) = name;
  }
  
  
  //! returns file name related to the token num_file
  inline const string& ParamOutputClass::GetFileName(int num_file) const
  {
    return name_file(num_file);
  }
   

  //! returns true if this output corresponds to a single point (seismogramm)
  inline bool ParamOutputClass::IsSeismogramOutput() const
  {
    return sismo_point;
  }
  
} // namespace Montjoie

#define MONTJOIE_FILE_PARAMETER_OUTPUT_INLINE_CXX
#endif

