#ifndef MONTJOIE_FILE_MESH_BOUNDARIES_CXX

namespace Montjoie
{
  
  //! Default constructor
  template <class Dimension, class T>
  PmlRegionParameter_Base<Dimension, T>::PmlRegionParameter_Base()
  {
    pml_add_axis_x = PML_NO;
    pml_add_axis_y = PML_NO;
    pml_add_axis_z = PML_NO;
    
    thicknessPML = 0.0;
    nb_layers_pml = 0;
    radiusPML = 0.0;
  }


  //! sets the PML layers to add along each coordinate
  template <class Dimension, class T>
  void PmlRegionParameter_Base<Dimension, T>::
  SetAdditionPML(int add_x, int add_y, int add_z, int nb_layers)
  {
    pml_add_axis_x = add_x;
    pml_add_axis_y = add_y;
    pml_add_axis_z = add_z;
    nb_layers_pml = nb_layers;
  }
  
  
  //! Sets area for PML with a line of the data file
  template <class Dimension, class T>
  void PmlRegionParameter_Base<Dimension, T>::SetParameters(const Vector<string>& parameters)
  {
    if (parameters.GetM() <= 3)
      {
	cout << "In SetInputData of MeshBoundaries_Base" << endl;
	cout << "AddPML needs more parameters, for instance :" << endl;
	cout << "AddPML = YES XYZ delta nb_layers" << endl;
	cout << "Current parameters are : " << endl << parameters << endl;
	abort();
      }
    
    if (!parameters(0).compare("YES"))
      {
	pml_add_axis_x = this->PML_NO;
	pml_add_axis_y = this->PML_NO;
	pml_add_axis_z = this->PML_NO;
	
	if (parameters(1).find("X-",0) != string::npos)
	  pml_add_axis_x = this->PML_NEGATIVE_SIDE;
	else if (parameters(1).find("X+",0) != string::npos)
	  pml_add_axis_x = this->PML_POSITIVE_SIDE;
	else if (parameters(1).find("X",0) != string::npos)
	  pml_add_axis_x = this->PML_BOTH_SIDES;
	
	if (parameters(1).find("Y-",0) != string::npos)
	  pml_add_axis_y = this->PML_NEGATIVE_SIDE;
	else if (parameters(1).find("Y+",0) != string::npos)
	  pml_add_axis_y = this->PML_POSITIVE_SIDE;
	else if (parameters(1).find("Y",0) != string::npos)
	  pml_add_axis_y = this->PML_BOTH_SIDES;
	
	if (parameters(1).find("Z-",0) != string::npos)
	  pml_add_axis_z = this->PML_NEGATIVE_SIDE;
	else if (parameters(1).find("Z+",0) != string::npos)
	  pml_add_axis_z = this->PML_POSITIVE_SIDE;
	else if (parameters(1).find("Z",0) != string::npos)
	  pml_add_axis_z = this->PML_BOTH_SIDES;
	
	if (parameters(1).find("R",0) != string::npos)
	  pml_add_axis_x = this->PML_RADIAL;
	
	thicknessPML = to_num<T>(parameters(2));
	if (!parameters(3).compare("AUTO"))
	  nb_layers_pml = 0;
	else
	  to_num(parameters(3), nb_layers_pml);
      }
    else
      {
	int pos_word = -1;
	for (int k = 0; k < parameters.GetM(); k++) 
	  if (isalpha(parameters(k)[0]))
	    {
	      pos_word = k;
	      break;
	    }

	if ((pos_word == -1) || (parameters.GetM() <= pos_word+2))
	  {
	    cout << "The provided parameters are incorrect "
		 << endl << parameters << endl;
	    
	    abort();
	  }
	
	thicknessPML = to_num<T>(parameters(pos_word+1));
	if (!parameters(pos_word+2).compare("AUTO"))
	  nb_layers_pml = 0;
	else
	  to_num(parameters(pos_word+2), nb_layers_pml);
	
	ref_pml.Reallocate(pos_word);
	for (int i = 0; i < pos_word; i++)
	  ref_pml(i) = to_num<int>(parameters(i));
      }
  }

  
  //! Adds PML layers
  template <class Dimension, class T>
  void PmlRegionParameter_Base<Dimension, T>::AddPML(int num_pml, Mesh<Dimension>& mesh)
  {
    if (this->nb_layers_pml == 0)
      {
	T length_edge = this->GetLeafClass().GetAverageLengthOnBoundary(mesh);
	this->nb_layers_pml = max(1, toInteger(ceil(this->thicknessPML/length_edge)));
      }
    
    this->GetLeafClass().AddPMLElements(this->nb_layers_pml, num_pml, mesh);
  }

  
  /***********************
   * MeshBoundaries_Base *
   ***********************/

  
  //! default constructor
  template <class Dimension, class T>
  MeshBoundaries_Base<Dimension, T>::MeshBoundaries_Base()
  {
    this->ResizeNbReferences(10);
  }
  
  
  //! returns size of memory used by the object in bytes
  template<class Dimension, class T>
  size_t MeshBoundaries_Base<Dimension, T>::GetMemorySize() const
  {
    size_t taille = Type_curve.GetMemorySize() + Cond_curve.GetMemorySize() + 
      Body_curve.GetMemorySize() + Param_known.GetMemorySize();

    taille += Seldon::GetMemorySize(Param_curve);
    taille += lob_curve.GetMemorySize();
    taille += periodicity_references.GetMemorySize();
    taille += type_coord_periodicity.GetMemorySize();
    return taille;
  }
  
  
  //! copies curves from another mesh
  template<class Dimension, class T>
  void MeshBoundaries_Base<Dimension, T>::CopyCurves(const Mesh<Dimension>& mesh)
  {
    Type_curve = mesh.Type_curve;
    Cond_curve = mesh.Cond_curve;
    Body_curve = mesh.Body_curve;
    Param_known = mesh.Param_known;
    Param_curve = mesh.Param_curve;
  }

  
  //! copies parameters from another mesh
  template<class Dimension, class T>
  void MeshBoundaries_Base<Dimension, T>::CopyInputData(const Mesh<Dimension>& mesh)
  {
    Mesh<Dimension, T>& mesh_leaf = this->GetLeafClass();
    
    CopyCurves(mesh);

    pml_areas = mesh.pml_areas;
    
    periodicity_references = mesh.periodicity_references;
    type_coord_periodicity = mesh.type_coord_periodicity;
    
    mesh_leaf.print_level = mesh.print_level;
    mesh_leaf.mesh_path = mesh.mesh_path;
    mesh_leaf.irregular_mesh_double_cut = mesh.irregular_mesh_double_cut;
    mesh_leaf.required_mesh_type = mesh.required_mesh_type;
    mesh_leaf.ListeVertices_to_be_refined = mesh.ListeVertices_to_be_refined;
    mesh_leaf.LevelRefinement_Vertex = mesh.LevelRefinement_Vertex;
    mesh_leaf.RatioRefinement_Vertex = mesh.RatioRefinement_Vertex;
    
  }
  
  
  //! returns an unused reference of the mesh
  template<class Dimension, class T>
  int MeshBoundaries_Base<Dimension, T>::GetNewReference()
  {
    Mesh<Dimension, T>& mesh_leaf = this->GetLeafClass();
    int ref = 3;
    for (int i = 0; i < mesh_leaf.GetNbBoundaryRef(); i++)
      ref = max(ref, mesh_leaf.BoundaryRef(i).GetReference());
  
    if ((ref+1) > this->GetNbReferences())
      ResizeNbReferences(ref+1);
    
    return (ref+1);
  }

  
  //! returns two unused reference of the mesh
  template<class Dimension, class T>
  void MeshBoundaries_Base<Dimension, T>::GetNewReference(int& new_ref1, int& new_ref2)
  {
    Mesh<Dimension, T>& mesh_leaf = this->GetLeafClass();
    int ref = 3;
    for (int i = 0; i < mesh_leaf.GetNbBoundaryRef(); i++)
      ref = max(ref, mesh_leaf.BoundaryRef(i).GetReference());
  
    if ((ref+2) > this->GetNbReferences())
      ResizeNbReferences(ref+2);
    
    new_ref1 = ref+1;
    new_ref2 = ref+2;
  }
  
  
  //! resizing arrays in order to be ready to a curve referenced ref_max
  /*!
    the arrays Type_curve, Cond_curve and Body_curve are resized
   */
  template<class Dimension, class T>
  void MeshBoundaries_Base<Dimension, T>::ResizeNbReferences(int ref_max)
  {
    Mesh<Dimension, T>& mesh = this->GetLeafClass();
    int nold = this->GetNbReferences();
    
    // we need to resize arrays
    Type_curve.Resize(ref_max+1);
    Cond_curve.Resize(ref_max+1);
    Body_curve.Resize(ref_max+1);
    Param_curve.Resize(ref_max+1);
    mesh.reference_original_neighbor.Resize(ref_max+1);

    Param_known.Resize(ref_max+1);
    // and set to zero new values
    for (int k = (nold+1); k <= ref_max; k++)
      {
	Type_curve(k) = 0;
	Cond_curve(k) = 0;
	Body_curve(k) = 0;
        mesh.reference_original_neighbor(k) = -3;
	Param_known(k) = false;
      }
  }

  
  //! sets a reference of the mesh to a given type of curve and boundary condition
  /*!
    \param[in] ref reference of a boundary
    \param[in] type curve type of the boundary (straight, circle, spline, sphere ...)
    \param[in] cond boundary condition associated to the boundary
   */
  template<class Dimension, class T>
  void MeshBoundaries_Base<Dimension, T>::SetNewReference(int ref, int type, int cond)
  {
    if (ref <= 0)
      {
        cout << "Reference should be greater than one " << endl;
        abort();
      }
    
    Type_curve(ref) = type;
    Cond_curve(ref) = cond;
  }
  

  //! copies boundary conditions of another mesh
  template <class Dimension, class T> template<class Dim>
  void MeshBoundaries_Base<Dimension, T>::SetBoundaryCondition(const Mesh<Dim, T>& mesh)
  {
    int nb_ref = mesh.GetNbReferences();
    if (this->GetNbReferences() < nb_ref)
      this->ResizeNbReferences(nb_ref+1);
    
    for (int i = 0; i <= nb_ref; i++)
      Cond_curve(i) = mesh.GetBoundaryCondition(i);
  }
    
  
  //! sets the boundary condition of reference i
  template <class Dimension, class T>
  void MeshBoundaries_Base<Dimension, T>::SetBoundaryCondition(int i, int j)
  {
    if (i <= 0)
      {
        cout << "You can't affect a boundary condition for references lower than 1 " << endl;
        abort();
      }
    
    Cond_curve(i) = j;
  }
  
  
  //! sets the type of curves for reference ref
  template <class Dimension, class T>
  void MeshBoundaries_Base<Dimension, T>::SetCurveType(int ref, int type)
  {
    if (ref <= 0)
      {
        cout << "You can't affect a curve type for references lower than 1 " << endl;
        abort();
      }
    
    if (ref > GetNbReferences())
      this->ResizeNbReferences(ref);
    
    Type_curve(ref) = type;
    Param_known(ref) = false;
  }
  
  
  //! sets the body number associated with reference i
  template <class Dimension, class T>
  void MeshBoundaries_Base<Dimension, T>::SetBodyNumber(int i, int num)
  {
    if (i <= 0)
      {
        cout << "You can't affect a body number for references lower than 1 " << endl;
        abort();
      }

    Body_curve(i) = num;
  }
  
  
  //! sets the parameters associated with reference i
  template <class Dimension, class T>
  void MeshBoundaries_Base<Dimension, T>::SetCurveParameter(int i, const Vector<T>& param)
  {
    if (i <= 0)
      {
        cout << "You can't affect a curve type for references lower than 1 " << endl;
        abort();
      }
    
    Param_curve(i) = param;
    if (param.GetM() >= 1)
      Param_known(i) = true;
    else
      Param_known(i) = false;
  }
  
  
  //! clear periodic references
  template<class Dimension, class T>
  void MeshBoundaries_Base<Dimension, T>::ClearPeriodicCondition()
  {
    type_coord_periodicity.Clear();
    periodicity_references.Clear();
  }
  
  
  //! modification of parameters with a line of the data file
  template<class Dimension, class T>
  void MeshBoundaries_Base<Dimension,T>::
  SetInputData(const string& description_field, const VectString& parameters)
  {
    Mesh<Dimension, T>& mesh_leaf = this->GetLeafClass();    
    if (!description_field.compare("TypeBody"))
      {
	if (parameters.GetM() <= 1)
	  {
	    cout << "In SetInputData of MeshBoundaries_Base" << endl;
	    cout << "TypeBody needs more parameters, for instance :" << endl;
	    cout << "TypeBody = ref body_num" << endl;
	    cout << "Current parameters are : " << endl << parameters << endl;
	    abort();
	  }

	int num = to_num<int>(parameters(parameters.GetM()-1));

	if (parameters(0)[0] == 'B')
	  {
	    int pos = parameters(0).find("F");
	    if (pos > 0)
	      {
		int beg = to_num<int>(parameters(0).substr(1, pos-1));
		int final = to_num<int>(parameters(0).substr(pos+1, parameters(0).size()-pos-1));
		
		if (final > this->GetNbReferences())
		  ResizeNbReferences(final);

		for (int k = beg; k <= final; k++)
		  this->SetBodyNumber(k, num);
	      }
	  }
	else
	  for (int k = 0; k < parameters.GetM()-1; k++)
	    {
	      int ref = to_num<int>(parameters(k));
	      if (ref > this->GetNbReferences())
		ResizeNbReferences(ref);
	      
	      this->SetBodyNumber(ref, num);
	    }	
      }
    else if (!description_field.compare("TypeCurve"))
      {
	int position_word = -1;
	// we look for the first parameter, with begins by a letter
	// This parameter is the description of the boundary condition (DIRICHLET, NEUMANN, etc)
	for (int i = 0; i < parameters.GetM(); i++)
	  if (isalpha(parameters(i)[0]))
	    {
              position_word = i;
              break;
            }
        
	if (position_word > 0)
	  {
	    // All the parameters at the left of the description of the boundary condition,
	    // are references. The references are stored
	    IVect ref(position_word);
	    for (int i = 0; i < position_word; i++)
	      {
		
		ref(i) = to_num<int>(parameters(i));
		if (ref(i) > this->GetNbReferences())
		  ResizeNbReferences(ref(i));
	      }
	    
            // setting description of the curve
	    mesh_leaf.ReadCurveType(ref, parameters(position_word));
            
	    // if some values appear after the description of the boundary condition,
	    // we consider that these values are parameters of the curve
	    if (parameters.GetM() > (position_word+1))
	      {
		// we store these values in parameters_curve
		Vector<T> parameters_curve(parameters.GetM()-position_word-1);
		for (int i = (position_word+1); i < parameters.GetM(); i++)
		  parameters_curve(i-position_word-1) = to_num<T>(parameters(i));
		
		// the parameters of the curve are known
		for (int i = 0; i < ref.GetM(); i++)
		  {
		    this->SetCurveParameter(ref(i), parameters_curve);
		    mesh_leaf.CheckCurveParameter(ref(i));
		  }
	      }	    
	  }
        else if (position_word == 0)
          {
            if (parameters(0)[0] == 'B')
              {
                int pos = parameters(0).find("F");
                if (pos > 0)
                  {
                    int beg = to_num<int>(parameters(0).substr(1, pos-1));
                    int final = to_num<int>(parameters(0).substr(pos+1, parameters(0).size()-pos-1));
                    IVect ref(final-beg+1);
                    for (int i = beg; i <= final; i++)
                      ref(i-beg) = i;
                    
                    if (final > this->GetNbReferences())
                      ResizeNbReferences(final);
                    
                    mesh_leaf.ReadCurveType(ref, parameters(1));
                  }
              }
            else
              {
                cout << "case not handled" << endl;
                abort();
              }
          }
      }
    else if (!description_field.compare("AddPML"))
      {
	PmlRegionParameter<Dimension, T> new_pml;
	new_pml.SetParameters(parameters);
	pml_areas.PushBack(new_pml);
      }
  }

      
  //! returns true if two curved boundaries have same set of parameters
  /*!
    \param[in] ref1 reference of the first boundary
    \param[in] ref2 reference of the second boundary
   */
  template<class Dimension, class T>
  bool MeshBoundaries_Base<Dimension, T>::AreCurveEqual(int ref1, int ref2) const
  {
    if (ref1 == ref2)
      return true;
    
    if (GetCurveType(ref1) == GetCurveType(ref2))
      {
        if ((ref1 <= 0)|| (ref2 <= 0))
          return true;
        
        if (Param_curve(ref1).GetM() == Param_curve(ref2).GetM())
          {
            for (int i = 0; i < Param_curve(ref1).GetM(); i++)
              if (abs(Param_curve(ref1)(i)-Param_curve(ref2)(i)) > R_N::threshold)
                return false;
            
            return true;
          }
      }
    
    return false;
  }
  
  
  //! returns reference edge number of extremities n1 and n2
  template<class Dimension, class T>
  int MeshBoundaries_Base<Dimension, T>::FindEdgeRefWithExtremities(int n1, int n2) const
  {
    const Mesh<Dimension, T>& mesh_leaf = this->GetLeafClass();
    int ne;
    if (n1 > n2)
      {
	// permutation of numbers so that we are sure to find edge in head_minv(n1)/next_edge
	ne = n1;
	n1 = n2;
	n2 = ne;
      }
    
    ne = -1;
    for (int e = mesh_leaf.head_minv(n1); e != -1; e = mesh_leaf.next_edge(e))
      if (mesh_leaf.edges(e).numVertex(1) == n2)
	ne = e;
    
    if (ne < mesh_leaf.GetNbEdgesRef())
      return ne;
    
    return -1;
  }
  
  
  //! returns the edge number of the edge [n1, n2]
  template<class Dimension, class T>
  int MeshBoundaries_Base<Dimension, T>::FindEdgeWithExtremities(int n1, int n2) const
  {
    const Mesh<Dimension, T>& mesh_leaf = this->GetLeafClass();
    int ne;
    if (n1 > n2)
      {
	// permutation of numbers so that we are sure to find edge in head_minv(n1)/next_edge
	ne = n1;
	n1 = n2;
	n2 = ne;
      }
    
    ne = -1;
    for (int e = mesh_leaf.head_minv(n1); e != -1; e = mesh_leaf.next_edge(e))
      if (mesh_leaf.edges(e).numVertex(1) == n2)
	ne = e;
    
    return ne;
  }
  
  
  //! clears linking informations between edges and faces
  template<class Dimension, class T>
  void MeshBoundaries_Base<Dimension,T>::ClearConnectivity()
  {
  }
  
  
  //! clears predefined curves
  template<class Dimension, class T>
  void MeshBoundaries_Base<Dimension,T>::ClearCurves()
  {
    // removing all the informations about curves
    Type_curve.Fill(0);
    Param_known.Fill(false);
    for (int i = 0; i < Param_curve.GetM(); i++)
      Param_curve(i).Clear();
  }
  
    
  //! clears the mesh
  template<class Dimension, class T>
  void MeshBoundaries_Base<Dimension,T>::Clear()
  {
  }

} // namespace Montjoie

#define MONTJOIE_FILE_MESH_BOUNDARIES_CXX
#endif
