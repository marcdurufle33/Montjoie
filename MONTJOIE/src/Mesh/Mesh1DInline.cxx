#ifndef MONTJOIE_FILE_MESH1D_INLINE_CXX

namespace Montjoie
{
  
  //! returns coordinates of the points
  template<class T>
  inline const Vector<T>& Mesh<Dimension1, T>::Vertex() const
  {
    return Vertices;
  }
  
  
  //! returns coordinate of the point i
  template<class T>
  inline const T& Mesh<Dimension1, T>::Vertex(int i) const
  {
    return Vertices(i);
  }
  
  
  //! returns coordinate of the point i
  template<class T>
  inline T& Mesh<Dimension1, T>::Vertex(int i)
  {
    return Vertices(i);
  }
  
  
  //! returns element i
  template<class T>
  inline const Edge<Dimension1>& Mesh<Dimension1, T>::Element(int i) const
  {
    return elements(i);
  }
  
  
  //! returns element i
  template<class T>
  inline Edge<Dimension1>& Mesh<Dimension1, T>::Element(int i)
  {
    return elements(i);
  }
  
  
  //! local to global number map
  template<class T>
  inline int Mesh<Dimension1, T>::GetNumberDof(int num_elt,int num_dof) const
  {
    return Nodle(num_elt,num_dof);
  }
  
  
  //! returns the two extremities of edge num_elt
  template<class T> template<class Vector1>
  inline void Mesh<Dimension1, T>::GetVerticesElement(int num_elt, Vector1& points) const
  {
    points(0) = Vertices(num_elt);
    points(1) = Vertices(num_elt+1);
  }
  
  
  //! returns the left extremity of physical domain
  template<class T>
  inline const T& Mesh<Dimension1, T>::GetXmin() const
  {
    return Vertices(0);
  }


  //! returns the right extremity of physical domain
  template<class T>
  inline const T& Mesh<Dimension1, T>::GetXmax() const
  {
    return Vertices(Vertices.GetM()-1);
  }
  
  
  //! returns the order associated with the mesh
  template<class T>
  inline int Mesh<Dimension1, T>::GetOrder() const
  {
    return order;
  }

  
  //! returns the order associated with the mesh
  template<class T>
  inline int Mesh<Dimension1, T>::GetOrderElement(int i) const
  {
    return order;
  }
  
  
  //! returns the number of degrees of freedom of the mesh
  template<class T>
  inline int Mesh<Dimension1, T>::GetNbDof() const
  {
    return nodl;
  }
  
  
  //! returns the number of elements in the mesh
  template<class T>
  inline int Mesh<Dimension1, T>::GetNbElt() const
  {
    return elements.GetM();
  }
  
  
  //! returns the number of vertices of the mesh
  template<class T>
  inline int Mesh<Dimension1, T>::GetNbVertices() const
  {
    return Vertices.GetM();
  }
  
  
  //! sets the order associated with the mesh
  template<class T>
  inline void Mesh<Dimension1, T>::SetOrder(int r)
  {
    order = r;
  }

  
  //! periodic condition is set
  template<class T>
  inline void Mesh<Dimension1, T>::SetPeriodicExtremity(bool per)
  {
    periodic_ext = per; 
  }
  
  
  //! returns the thickness of PML
  template<class T>
  inline const T& Mesh<Dimension1, T>::GetThicknessPML() const
  {
    return thicknessPML;
  }
  
  
  //! sets the thickness of PML
  template<class T>
  inline void Mesh<Dimension1, T>::SetThicknessPML(const T& delta)
  {
    thicknessPML = delta;
  }
  
  
  //! returns how PML have to be added to the mesh
  template<class T>
  inline int Mesh<Dimension1, T>::GetTypeAdditionPML(int) const
  {
    return pml_add_type;
  }
  
  
  //! returns the number of cells in the direction of PML
  template<class T>
  inline int Mesh<Dimension1, T>::GetNbLayersPML() const
  {
    return nb_layers_pml;
  }
  
  
  //! sets PML to be added to the mesh
  template<class T>
  inline void Mesh<Dimension1, T>::SetAdditionPML(int type_add, int nb_layers, int ref)
  {
    pml_add_type = type_add;
    nb_layers_pml = nb_layers;
    ref_pml = ref;
  }
  
}

#define MONTJOIE_FILE_MESH1D_INLINE_CXX
#endif
