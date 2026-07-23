#ifndef MONTJOIE_FILE_TETRAHEDRON_HCURL_LOBATTO_CXX

namespace Montjoie
{
  
  TetrahedronHcurlLobatto::TetrahedronHcurlLobatto()
    : elt_hcurl(tet_h1, *this)
  {
    this->use_piola_transform = false;
  }


  void TetrahedronHcurlLobatto::ConstructNumberMap(NumberMap& nmap, int dg) const
  {
    // dg formulation
    nmap.SetFormulationDG(dg);

    if (dg == ElementReference_Base::DISCONTINUOUS)
      return TetrahedronReference<2>::ConstructNumberMap(nmap, dg);
    else if (dg == ElementReference_Base::HDG)
      {
	nmap.SetNbDofVertex(this->order, 0);
	nmap.SetNbDofEdge(this->order, 0);
	nmap.SetNbDofTriangle(this->order, (this->order+1)*(this->order+2));
	nmap.SetNbDofTetrahedron(this->order, 0);
	
	// rotation of dofs on faces
	element_tri_surf->FindHcurlRotationTri(nmap, 0);
	
	return;
      }
    else
      {
	cout << "Not implemented for continuous elements" << endl;
	abort();
      }
  }
    

  // construction of finite element
  void TetrahedronHcurlLobatto
  ::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
			   int rsurf_tri, int rsurf_quad,
			   int type_surf_tri, int type_surf_quad, int gauss_z)
  {
    tet_h1.ConstructFiniteElement(r, rgeom, rquad, type_quad,
				  rsurf_tri, rsurf_quad, type_surf_tri, type_surf_quad, gauss_z);

    this->Fb_geom = tet_h1.GetGeometricTetrahedron();

    this->element_tri_surf = new TriangleHcurlSecondFamily();
    this->element_tri_surf->ConstructFiniteElement(r);
    
    elt_hcurl.ConstructFiniteElement();

    this->element_scal_vol = &tet_h1;
  }
  

  void TetrahedronHcurlLobatto::ComputeValuesPhiRef(const R3& x, VectR3& phi) const
  {
    elt_hcurl.ComputeValuesPhiRef(x, phi);
  }

  
  void TetrahedronHcurlLobatto::ComputeCurlPhiRef(const R3&, VectR3&) const
  {
    cout << "not useful without Piola transform" << endl;
    abort();
  }

  void TetrahedronHcurlLobatto::GetValuePhiOnQuadraturePoint(int k, VectR3& phi) const
  {
    elt_hcurl.GetValuePhiOnQuadraturePoint(k, phi);
  }
  
  void TetrahedronHcurlLobatto::ApplyCh(const VectReal_wp& U, VectReal_wp& V) const
  {
    elt_hcurl.ApplyChGen(U, V);
  }


  void TetrahedronHcurlLobatto::ApplyCh(const VectComplex_wp& U, VectComplex_wp& V) const
  {
    elt_hcurl.ApplyChGen(U, V);
  }
  
  
  ostream& operator <<(ostream& out, const TetrahedronHcurlLobatto& e)
  {
    return out;
  }

} // namespace Montjoie

#define MONTJOIE_FILE_TETRAHEDRON_HCURL_LOBATTO_CXX
#endif
