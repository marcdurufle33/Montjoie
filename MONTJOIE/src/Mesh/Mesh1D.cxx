#ifndef MONTJOIE_FILE_MESH_1D_CXX

#include "MeshBase.hxx"

namespace Montjoie
{
  
  //! default constructor
  template<class T>
  Mesh<Dimension1, T>::Mesh()
  {
    order = 1;
    nodl = 0;
    print_level = 0;
    periodic_ext = false;
    pml_add_type = PML_NO;
    nb_layers_pml = 0;
    thicknessPML = 0.0;
    ref_pml = 0;
  }
  

  //! returns the size used by the object in bytes
  template<class T>
  size_t Mesh<Dimension1, T>::GetMemorySize() const
  {
    size_t taille = sizeof(*this) - 10*sizeof(Vector<T>);
    taille += Nodle.GetMemorySize();
    taille += Vertices.GetMemorySize();
    taille += elements.GetMemorySize();
    taille += level_refinement_vertex.GetMemorySize();
    taille += position_refinement_vertex.GetMemorySize();
    taille += ratio_refinement_vertex.GetMemorySize();
    taille += nb_layers_to_add.GetMemorySize();
    taille += reference_layers_to_add.GetMemorySize();
    taille += position_layers_to_add.GetMemorySize();
    taille += GlobDofNumber_Subdomain.GetMemorySize();
    return taille;
  }

  
  //! returns the mesh size (average length of edges)
  template<class T>
  T Mesh<Dimension1, T>::GetMeshSize() const
  {
    T sum = 0;
    for (int i = 0; i < elements.GetM(); i++)
      sum += this->Vertex(i+1)- this->Vertex(i);
    
    sum /= this->elements.GetM();
    return sum;
  }
  
  
  //! numbers the mesh
  /*!
    map.nb_dof_vertex : number of dofs for each vertex
    map.nb_dof_edge : number of dofs strictely inside each edge
  */
  template<class T>
  void Mesh<Dimension1, T>::NumberMesh(bool discontinuous_mesh)
  {
    int nb_dof_edge = order+1;
    int num_dof = 0;
    int nb_edges = elements.GetM();
    if (discontinuous_mesh)
      {
	nodl = nb_edges*nb_dof_edge;   
	Nodle.Reallocate(nb_edges, nb_dof_edge);
	
	for (int i = 0; i < nb_edges; i++)
	  for (int j = 0; j <= order; j++)
	    Nodle(i, j) = num_dof++;
      }
    else
      {
	nodl = nb_edges*(nb_dof_edge-1)+1;   
	Nodle.Reallocate(nb_edges, nb_dof_edge);
	
	for (int i = 0; i < nb_edges; i++)
	  {
	    Nodle(i,0) = num_dof;
	    for (int j = 1; j < (nb_dof_edge-1); j++)
	      Nodle(i,j) = num_dof+j;
	    
	    Nodle(i, nb_dof_edge-1) = num_dof + nb_dof_edge - 1;
	    num_dof += nb_dof_edge-1;
	  }
        
        if (periodic_ext)
          {
            Nodle(nb_edges-1, nb_dof_edge-1) = 0;
            nodl--;
          }
      }
  }
  
  
  //! reads a mesh in a file
  template<class T>
  void Mesh<Dimension1, T>::Read(const string & file_mesh)
  {
    Vertices.ReadText(file_mesh);
    elements.Reallocate(Vertices.GetM()-1);
    for (int i = 0; i < Vertices.GetM()-1; i++)
      elements(i).Init(i, i+1, 1);
  }
  
  
  //! writes a mesh (only coordinates are written)
  template<class T>
  void Mesh<Dimension1, T>::Write(const string & file_mesh)
  {
    Vertices.WriteText(file_mesh);
  }
  
  
  //! subdivision of a mesh 
  /*!
    \param[out] fine_mesh result mesh
    \param[in] factor_refinement number of subdivisions for each reference
   */
  template<class T>
  void Mesh<Dimension1, T>::SubdivideMesh(Mesh<Dimension1, T>& fine_mesh,
                                          const IVect& factor_refinement) const
  {
    int nb_edges_fine = 0;
    int j, nbinterv;
    for (int i = 0; i < elements.GetM(); i++)
      {
	j = elements(i).GetReference();
	nb_edges_fine += factor_refinement(j);
      }
    fine_mesh.Vertices.Reallocate(nb_edges_fine+1);
    fine_mesh.elements.Reallocate(nb_edges_fine);
    int nb_vertices = 1;
    fine_mesh.Vertices(0) = Vertices(0);
    for (int i = 0; i < elements.GetM(); i++)
      {
	j = elements(i).GetReference();
	nbinterv = factor_refinement(j);
	T deltax = (Vertices(i+1)-Vertices(i))/nbinterv;
	T pos_x = Vertices(i);
	for (j = 1; j <= nbinterv; j++)
	  {
	    pos_x += deltax;
	    fine_mesh.Vertices(nb_vertices)=pos_x;
	    fine_mesh.elements(nb_vertices-1)
              .Init(nb_vertices-1, nb_vertices, elements(i).GetReference());
	    nb_vertices++;
	  }
      }
  }
  
  
  //! creation of a non-regular mesh
  /*!
    \param[in] Points vertices
    \param[in] ref_domain reference for each element
   */
  template<class T>
  void Mesh<Dimension1, T>::CreateNonRegularMesh(const Vector<T>& Points, const IVect& ref_domain)
  {    
    Vertices = Points;
    elements.Reallocate(Vertices.GetM() - 1);
    for (int i = 0; i < Vertices.GetM()-1; i++)
      elements(i).Init(i, i+1, ref_domain(i));
  }
  
  
  //! creation of a regular mesh
  /*!
    \param[in] nbPoints_x number of points
    \param[in] ref_domain reference for all the edges
   */
  template<class T>
  void Mesh<Dimension1, T>::CreateRegularMesh(const T& xmin, const T& xmax,
                                              int nbPoints_x, int ref_domain)
  {
    T step_x;
    if (nbPoints_x <=1 )
      step_x = xmax-xmin;
    else
      step_x = (xmax-xmin)/(nbPoints_x-1);
        
    Vertices.Reallocate(nbPoints_x);
    elements.Reallocate(nbPoints_x-1);
    for (int i = 0; i < nbPoints_x ; i++)
      {
	Vertices(i) = xmin + step_x*i;
	if (i < (nbPoints_x-1))
	  elements(i).Init(i, i+1, ref_domain);
      }
  }
  
  //! creates a layered mesh 
  template<class T>
  void Mesh<Dimension1, T>::CreateLayeredMesh(const Vector<T>& pos, const IVect& nb_cells, const IVect& ref)
  {
    T step_x;    
    int nbPoints = 1;
    for (int i = 0; i < nb_cells.GetM(); i++)
      nbPoints += nb_cells(i);
    
    Vertices.Reallocate(nbPoints);
    elements.Reallocate(nbPoints-1);
    int nb = 0;
    for (int i = 0; i < nb_cells.GetM(); i++)
      {
        step_x = pos(i+1)-pos(i);
        if (nb_cells(i) > 1)
          step_x /= nb_cells(i);
        
        for (int j = 0; j < nb_cells(i); j++)
          {
            Vertices(nb) = pos(i) + step_x*j;
            elements(nb).Init(nb, nb+1, ref(i));
            nb++;
          }
      }
    
    // last point
    Vertices(nb) = pos(nb_cells.GetM());
  }
  

  //! adds a layer with N cells until x = pos and with a reference ref
  //! these elements are added when ConstructMesh is called
  template<class T>
  void Mesh<Dimension1, T>::AddLayersInput(const T& pos, int N, int ref)
  {
    position_layers_to_add.PushBack(pos);
    nb_layers_to_add.PushBack(N);
    reference_layers_to_add.PushBack(ref);
  }

  
  //! Clears the layers to be added when ConstructMesh is called
  template<class T>
  void Mesh<Dimension1, T>::ClearLayersInput()
  {
    position_layers_to_add.Clear();
    nb_layers_to_add.Clear();
    reference_layers_to_add.Clear();
  }
  

  //! Adds a refinement close to the vertex pos, the refinement is performed when ConstructMesh is called
  /*!
    \param[in] pos vertex coordinate
    \param[in] lvl refinement level
    \param[in] ratio refinement ratio
  */
  template<class T>
  void Mesh<Dimension1, T>::AddRefinementVertex(const T& pos, int lvl, const T& ratio)
  {
    position_refinement_vertex.PushBack(pos);
    level_refinement_vertex.PushBack(lvl);
    ratio_refinement_vertex.PushBack(ratio);
  }


  //! Clears all refinements to be performed when ConstructMesh is called
  template<class T>
  void Mesh<Dimension1, T>::ClearRefinementVertex()
  {
    position_refinement_vertex.Clear();
    level_refinement_vertex.Clear();
    ratio_refinement_vertex.Clear();
  }
  
  
  //! constructs mesh from parameters of the data file
  template<class T>
  void Mesh<Dimension1, T>::ConstructMesh(const Vector<string>& parameters,
					  Real_wp& xmin_, Real_wp& xmax_)
  {
    if (parameters.GetM() <= 0)
      {
	cout << "Empty mesh, write a line FileMesh = something in the data file" << endl;
	abort();
      }
    
    int nb = 1;
    if (!parameters(0).compare("REGULAR"))
      {
	if (parameters.GetM() <= 3)
	  {
	    cout << "FileMesh requires more parameters. A good syntax is" << endl;
	    cout << "FileMesh = REGULAR nbx xmin xmax" << endl;
	    cout << "Current parameters are " << parameters << endl;
	    abort();
	  }
	
	int nbPoints_square = to_num<int>(parameters(1));
	Real_wp xmin = to_num<Real_wp>(parameters(2));
	Real_wp xmax = to_num<Real_wp>(parameters(3));
	this->CreateRegularMesh(xmin, xmax, nbPoints_square, 1);
        nb = 4;
      }
    else if (parameters(0) == "LAYERED")
      {
        if (parameters.GetM() <= 1)
	  {
	    cout << "FileMesh requires more parameters. A good syntax is" << endl;
	    cout << "FileMesh = LAYERED r0 r1 r2 .. r_N AUTO step" << endl;
	    cout << "Current parameters are " << parameters << endl;
	    abort();
	  }
	
        int position_word = -1;
        // we look for the first parameter that begins by a letter
	// This parameter is the description of the space step strategy
	for (int i = 1; i < parameters.GetM(); i++)
	  if (isalpha(parameters(i)[0]))
	    {
	      position_word = i;
	      break;
	    }

        if (position_word > 2)
          {
            // reading positions r0, r1, ..., r_N
            VectReal_wp pos(position_word-1);
            IVect nbSubdiv(position_word-2);
            nbSubdiv.Fill(1);
            for (int i = 1; i < position_word; i++)
              pos(i-1) = to_num<Real_wp>(parameters(i));
            
            for (int i = 0; i < pos.GetM()-1; i++)
              if (pos(i) > pos(i+1))
                {
                  cout << "The positions on the 1-D mesh must be increasing" << endl;
                  abort();
                }
            
            nb = position_word+1;            
            if (parameters(position_word) == "AUTO")
              {
                // computing the number of cells in each layer
                Real_wp dx = to_num<Real_wp>(parameters(nb++));
                for (int i = 0; i < nbSubdiv.GetM(); i++)
                  nbSubdiv(i) = toInteger(ceil((pos(i+1)-pos(i))/dx));                
              }
            else
              {
                // manual choice => the user provides the number of cells for each layer
                for (int i = 0; i < nbSubdiv.GetM(); i++)
                  nbSubdiv(i) = to_num<int>(parameters(nb++));
              }
            
            // creation of this mesh with references 1, 2, 3, etc
            IVect ref_domain(position_word-2);
            for (int i = 0; i < nbSubdiv.GetM(); i++)
              ref_domain(i) = i+1;
            
            this->CreateLayeredMesh(pos, nbSubdiv, ref_domain);
          }
      }
    else
      {
	string file_mesh = parameters(0);
	this->Read(file_mesh);
      }
    
    if (parameters.GetM() > nb+1)
      {
        if (parameters(nb) == "REFINED")
          {
            int r = to_num<int>(parameters(nb+1));
            Mesh<Dimension1> mesh(*this);
            IVect factor_refinement(mesh.GetNbElt());
            factor_refinement.Fill(r);
            mesh.SubdivideMesh(*this, factor_refinement);
          }
      }
        
    for (int i = 0; i < level_refinement_vertex.GetM(); i++)
      {
        // refinement around a vertex
        int lvl = level_refinement_vertex(i);
        Real_wp ratio = ratio_refinement_vertex(i);
        Real_wp pt = position_refinement_vertex(i);
        
        // we search the point
        int num = -1;
        for (int j = 0; j < this->Vertices.GetM(); j++)
          if (abs(pt - this->Vertices(j)) < 1e-10)
            {
              num = j;
              break;
            }
        
        // if found, the local refinement is performed
        if (num >= 0)
          {
            Vector<T> OldVertex(this->Vertices);
            Vector<Edge<Dimension1> > old_elt(elements);
            
            int nb_vert = 0;
            if ((num == 0) || (num == OldVertex.GetM()-1))
              nb_vert = lvl;
            else
              nb_vert = 2*lvl;
            
            Vertices.Reallocate(OldVertex.GetM() + nb_vert);
            elements.Reallocate(old_elt.GetM() + nb_vert);
            
            // left part of the mesh
            for (int i = 0; i < num; i++)
              {
                Vertices(i) = OldVertex(i);
                elements(i).Init(i, i+1, old_elt(i).GetReference());
              }
                
            int nb = num;
            // refined part on the left of point
            if (num > 0)
              {
                Real_wp a = OldVertex(num-1), b = OldVertex(num);
                Real_wp lambda = 1.0/ratio;
                for (int i = 0; i < lvl; i++)
                  {
                    Vertices(nb) = (1.0-lambda)*b + lambda*a;
                    elements(nb).Init(nb, nb+1, old_elt(num-1).GetReference());
                    nb++;
                    lambda /= ratio;
                  }
              }
            
            // target vertex 
            Vertices(nb) = OldVertex(num);
            if (num < OldVertex.GetM()-1)
              elements(nb).Init(nb, nb+1, old_elt(num).GetReference());
            
            nb++;
            // refined part on the right of point
            if (num < OldVertex.GetM()-1)
              {
                Real_wp a = OldVertex(num), b = OldVertex(num+1);
                Real_wp lambda = pow(1.0/ratio, lvl);
                for (int i = 0; i < lvl; i++)
                  {
                    Vertices(nb) = (1.0-lambda)*a + lambda*b;
                    elements(nb).Init(nb, nb+1, old_elt(num).GetReference());
                    nb++;
                    lambda *= ratio;
                  }
              }                
            
            // remaining elements on the right
            for (int i = num+1; i < OldVertex.GetM(); i++)
              {
                Vertices(nb) = OldVertex(i);
                if (i < OldVertex.GetM()-1)
                  elements(nb).Init(nb, nb+1, old_elt(i).GetReference());
                
                nb++;
              }
          }
	else
	  {
	    cout << "Point " << pt << " Not found in the mesh" << endl;
	    cout << "Refinement cannot be performed" << endl;
	    abort();
	  }
      }
    
    // additional layers
    for (int p = 0; p < nb_layers_to_add.GetM(); p++)
      {
	int nb_div = this->nb_layers_to_add(p);
	Real_wp xf = this->position_layers_to_add(p);
	int ref = this->reference_layers_to_add(p);

	int nb = Vertices.GetM();
	Real_wp xmax = Vertices(nb-1);
	Real_wp dx = (xf - xmax) / nb_div;
	
	Vertices.Resize(nb + nb_div);
	elements.Resize(elements.GetM() + nb_div);
	for (int i = 0; i < nb_div; i++)
	  {
	    if (i == nb_div-1)
	      Vertices(nb) = xf;
	    else
	      Vertices(nb) = xmax + (i+1)*dx;
	    
	    elements(nb-1).Init(nb-1, nb, ref);
	    nb++;
	  }
      }

    xmin_ = 0; xmax_ = 0;
    if (this->Vertices.GetM() > 0)
      {
	xmin_ = this->Vertices(0);
	xmax_ = this->Vertices(this->Vertices.GetM()-1);
      }
    
    // PML layers
    if (this->pml_add_type != PML_NO)
      {
	if (nb_layers_pml == 0)
	  {
	    int N = this->Vertices.GetM();
	    Real_wp h_left(0), h_right(0);
	    if (this->Vertices.GetM() >= 2)
	      {
		h_left = Vertices(1) - Vertices(0);
		h_right = Vertices(N-1) - Vertices(N-2);
	      }
	    else
	      {
		abort();
	      }
	    
	    Real_wp h = 0.5*(h_left + h_right);
	    nb_layers_pml = max(1, toInteger(round(thicknessPML/h)));
	  }
	
	this->AddPMLElements(nb_layers_pml);
      }
  }
  
  
  //! creates an extracted mesh with elements between n0 and n1
  template<class T>
  void Mesh<Dimension1, T>::CreateSubmesh(Mesh<Dimension1, T>& sub_mesh, int n0, int n1) const
  {
    if (n1 <= n0)
      {
        sub_mesh.Clear();
        return;
      }
    
    int nb_elt = n1 - n0;
    sub_mesh.Vertices.Reallocate(nb_elt+1);
    sub_mesh.elements.Reallocate(nb_elt);
    sub_mesh.nodl = 1;
    for (int i = n0; i < n1; i++)
      sub_mesh.nodl += this->GetOrderElement(i);
        
    sub_mesh.order = order;
    sub_mesh.Nodle.Reallocate(nb_elt, order+1);
    sub_mesh.GlobDofNumber_Subdomain.Reallocate(sub_mesh.nodl);
    sub_mesh.Vertices(0) = Vertices(n0);
    sub_mesh.nodl = 0;
    for (int i = n0; i < n1; i++)
      {
        sub_mesh.Vertices(i-n0+1) = Vertices(i+1);
        sub_mesh.elements(i-n0).Init(i-n0, i-n0+1, elements(i).GetReference());
        if (elements(i).IsPML())
          sub_mesh.elements(i-n0).SetPML();
        
        int r = this->GetOrderElement(i);
        for (int j = 0; j <= r; j++)
          {
            sub_mesh.Nodle(i-n0, j) = sub_mesh.nodl;
            sub_mesh.GlobDofNumber_Subdomain(sub_mesh.nodl) = Nodle(i, j);
            if (j < r)
              sub_mesh.nodl++;
          }
      }
    
    sub_mesh.nodl++;
    
    sub_mesh.periodic_ext = false;
  }
  
  
  //! clears the mesh
  template<class T>
  void Mesh<Dimension1, T>::Clear()
  {
    Nodle.Clear(); Vertices.Clear();
    elements.Clear();
    nodl = 0;
    GlobDofNumber_Subdomain.Clear();
  }
  
  
  //! Adds elements on the left and/or right side (PML elements)
  template<class T>
  void Mesh<Dimension1, T>::AddPMLElements(int nb_div)
  {
    int nb_new_elt(0);
    if ((pml_add_type == PML_NEGATIVE_SIDE) || (pml_add_type == PML_POSITIVE_SIDE))
      nb_new_elt = nb_div;
    
    if (pml_add_type == PML_BOTH_SIDES)
      nb_new_elt = 2*nb_div;
    
    if (nb_new_elt <= 0)
      return;
    
    Vector<T> OldVertex(Vertices);
    Vector<Edge<Dimension1> > OldElt(elements);
    
    Vertices.Reallocate(OldVertex.GetM() + nb_new_elt);
    elements.Reallocate(OldElt.GetM() + nb_new_elt);
    T dx = thicknessPML/nb_div;
    T xmin = OldVertex(0) - thicknessPML;
    int ref = OldElt(0).GetReference();
    if (ref_pml > 0)
      ref = ref_pml;
    
    int nb = 0;
    // elements in the left PML
    if ((pml_add_type == PML_NEGATIVE_SIDE) || (pml_add_type == PML_BOTH_SIDES))
      for (int i = 0; i < nb_div; i++)
	{
	  Vertices(i) = xmin + i*dx;
	  elements(i).Init(i, i+1, ref);
	  elements(i).SetPML();
	  nb++;
	}
    
    // physical elements
    for (int i = 0; i < OldElt.GetM(); i++)
      {
	Vertices(nb) = OldVertex(i);
	elements(nb).Init(nb, nb+1, OldElt(i).GetReference());
	ref = OldElt(i).GetReference();
	nb++;
      }

    if (ref_pml > 0)
      ref = ref_pml;
    
    // elements in the right PML
    T xmax = OldVertex(OldVertex.GetM()-1);
    if ((pml_add_type == PML_POSITIVE_SIDE) || (pml_add_type == PML_BOTH_SIDES))
      {
	for (int i = 0; i < nb_div; i++)
	  {
	    Vertices(nb) = xmax + i*dx;
	    elements(nb).Init(nb, nb+1, ref);
	    elements(nb).SetPML();
	    nb++;
	  }
	
	xmax += thicknessPML;
      }
        
    Vertices(nb) = xmax;
  }
  
  
  //! prints details of class Mesh1D
  template<class T>
  ostream& operator <<(ostream& out, const Mesh<Dimension1, T>& mesh)
  {
    out<<"The mesh contains "<<mesh.GetNbVertices()<<" vertices and "
       <<mesh.GetNbElt()<<" elements "<<endl;
    out<<"Coordinate of vertices "<<endl;
    for (int i = 0; i < mesh.GetNbVertices(); i++)
      out<<mesh.Vertex(i)<<endl;
    
    out<<"Global numbering "<<endl;
    for (int i = 0; i < mesh.GetNbElt(); i++)
      {
	out << "Element " << i << endl;
	for (int  j = 0; j <= mesh.GetOrder(); j++)
	  out<<mesh.GetNumberDof(i, j)<<"  ";
	
	out << endl;
      }
    return out;
  }
  
}

#define MONTJOIE_FILE_MESH_1D_CXX
#endif
