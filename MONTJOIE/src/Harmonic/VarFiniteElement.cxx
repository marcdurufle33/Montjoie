#ifndef MONTJOIE_FILE_VAR_FINITE_ELEMENT_CXX

namespace Montjoie
{  
  
  /********************
   * VarFiniteElement *
   ********************/
  
#ifdef MONTJOIE_WITH_NODAL_H1

#ifdef MONTJOIE_WITH_ONE_DIM
  //! Inits static attibutes containing the list of available finite elements
  void VarFiniteElementEnum<Dimension1, 1>::InitStaticData()
  {
    list_element.Reallocate(1);
    list_element(0)["EDGE_GAUSS"] = _EdgeGauss;
    list_element(0)["EDGE_LOBATTO"] = _EdgeLobatto;
    list_element(0)["EDGE_RADAU"] = _EdgeRadau;
    list_element(0)["EDGE_HIERARCHIC"] = _EdgeHierarchic;
  }


  //! returns a new pointer to a finite element of a given id
  ElementReference<Dimension1, 1>* VarFiniteElementEnum<Dimension1, 1>
  ::GetNewReferenceElement(int id)
  {    
    switch (id)
      {
      case _EdgeGauss :
	return new EdgeGaussReference();
      case _EdgeLobatto :
	return new EdgeLobattoReference();
      case _EdgeRadau :
	return new EdgeRadauReference();
      case _EdgeHierarchic :
	return new EdgeHierarchicReference();
      default :
	{
	  cout << "Unknown finite element" << endl;
	  abort();
	}
      }
    
    return NULL;
  }
#endif

#ifdef MONTJOIE_WITH_TWO_DIM
  //! Inits static attibutes containing the list of available finite elements
  void VarFiniteElementEnum<Dimension2, 1>::InitStaticData()
  {
    list_element.Reallocate(2);
    list_element(0)["TRIANGLE_CLASSICAL"] = _TriangleClassical;
    list_element(1)["TRIANGLE_CLASSICAL"] = _QuadrangleGauss;

    list_element(0)["TRIANGLE_LOBATTO"] = _TriangleLobatto;
    list_element(1)["TRIANGLE_LOBATTO"] = _QuadrangleLobatto;
    list_element(0)["QUADRANGLE_LOBATTO"] = _TriangleLobatto;
    list_element(1)["QUADRANGLE_LOBATTO"] = _QuadrangleLobatto;

    list_element(0)["TRIANGLE_HIERARCHIC"] = _TriangleHierarchic;
    list_element(1)["TRIANGLE_HIERARCHIC"] = _QuadrangleHierarchic;

    list_element(0)["TRIANGLE_DG_LOBATTO"] = _TriangleLobatto;
    list_element(1)["TRIANGLE_DG_LOBATTO"] = _QuadrangleLobatto;

    list_element(0)["TRIANGLE_DG_LOBATTO_GAUSS"] = _TriangleClassical;
    list_element(1)["TRIANGLE_DG_LOBATTO_GAUSS"] = _QuadrangleGauss;

    list_element(0)["TRIANGLE_DG_CLASSICAL"] = _TriangleClassical;
    list_element(1)["TRIANGLE_DG_CLASSICAL"] = _QuadrangleDgGauss;
    list_element(0)["QUADRANGLE_DG_GAUSS"] = _TriangleClassical;
    list_element(1)["QUADRANGLE_DG_GAUSS"] = _QuadrangleDgGauss;
    
    list_element(0)["TRIANGLE_DG_ORTHO"] = _TriangleDgOrtho;
    list_element(1)["TRIANGLE_DG_ORTHO"] = _QuadrangleDgGauss;

    list_element(0)["TRIANGLE_RADAU"] = _TriangleLobatto;
    list_element(1)["TRIANGLE_RADAU"] = _QuadrangleRadau;
  }


  //! returns a new pointer to a finite element of a given id  
  ElementReference<Dimension2, 1>* VarFiniteElementEnum<Dimension2, 1>
  ::GetNewReferenceElement(int id)
  {    
    switch (id)
      {
      case _TriangleClassical :
	return new TriangleClassical();
      case _QuadrangleGauss :
	return new QuadrangleGauss();
      case _TriangleLobatto :
	return new TriangleLobatto();
      case _QuadrangleLobatto :
	return new QuadrangleLobatto();
#ifdef MONTJOIE_WITH_HP_H1
      case _TriangleHierarchic :
	return new TriangleHierarchic();
      case _QuadrangleHierarchic :
	return new QuadrangleHierarchic();
#endif
#ifdef MONTJOIE_WITH_NODAL_DG
      case _QuadrangleDgGauss :
	return new QuadrangleDgGauss();
#endif
#ifdef MONTJOIE_WITH_ORTHO_DG
      case _TriangleDgOrtho :
	return new TriangleDgOrtho();
#endif
      case _QuadrangleRadau :
	{
	  QuadrangleRadau* elt;
	  elt = new QuadrangleRadau();
	  elt->SetRadauInterpolationPoints();
	  return elt;
	}
      default :
	{
	  cout << "Unknown finite element" << endl;
	  abort();
	}
      }
    
    return NULL;
  }

  //! returns a new pointer to a surface finite element of order r and of a given id  
  ElementReference<Dimension1, 1>* VarFiniteElementEnum<Dimension2, 1>
  ::GetNewSurfaceElement(int r, int rgeom, int rquad, int id)
  {    
    switch (id)
      {
      case _TriangleClassical :
      case _TriangleDgOrtho :
	{
	  EdgeGaussReference* edge = new EdgeGaussReference();
	  edge->ConstructFiniteElement(r, rgeom, rquad, EdgeGauss::GAUSS);
	  return edge;
	}
      case _TriangleLobatto :
	{
	  EdgeLobattoReference* edge = new EdgeLobattoReference();
	  edge->ConstructFiniteElement(r, rgeom, rquad);	  
	  return edge;
	}
      case _TriangleHierarchic :
	{
	  EdgeHierarchicReference* edge = new EdgeHierarchicReference();
	  edge->ConstructFiniteElement(r, rgeom, rquad);	  
	  return edge;
	}
      default :
	{
	  cout << "Unknown finite element" << endl;
	  abort();
	}
      }
    
    return NULL;
  }
#endif

#ifdef MONTJOIE_WITH_THREE_DIM
  //! Inits static attibutes containing the list of available finite elements
  void VarFiniteElementEnum<Dimension3, 1>::InitStaticData()
  {
    list_element.Reallocate(4);
    list_element(0)["TETRAHEDRON_CLASSICAL"] = _TetrahedronClassical;
    list_element(1)["TETRAHEDRON_CLASSICAL"] = _PyramidClassical;
    list_element(2)["TETRAHEDRON_CLASSICAL"] = _WedgeClassical;
    list_element(3)["TETRAHEDRON_CLASSICAL"] = _HexahedronGauss;

    list_element(0)["TETRAHEDRON_LOBATTO"] = _TetrahedronClassical;
    list_element(1)["TETRAHEDRON_LOBATTO"] = _PyramidLobatto;
    list_element(2)["TETRAHEDRON_LOBATTO"] = _WedgeLobatto;
    list_element(3)["TETRAHEDRON_LOBATTO"] = _HexahedronLobatto;
    list_element(0)["HEXAHEDRON_LOBATTO"] = _TetrahedronClassical;
    list_element(1)["HEXAHEDRON_LOBATTO"] = _PyramidLobatto;
    list_element(2)["HEXAHEDRON_LOBATTO"] = _WedgeLobatto;
    list_element(3)["HEXAHEDRON_LOBATTO"] = _HexahedronLobatto;

    list_element(0)["TETRAHEDRON_HIERARCHIC"] = _TetrahedronHierarchic;
    list_element(1)["TETRAHEDRON_HIERARCHIC"] = _PyramidHierarchic;
    list_element(2)["TETRAHEDRON_HIERARCHIC"] = _WedgeHierarchic;
    list_element(3)["TETRAHEDRON_HIERARCHIC"] = _HexahedronHierarchic;

    list_element(0)["TETRAHEDRON_DG_LOBATTO"] = _TetrahedronClassical;
    list_element(1)["TETRAHEDRON_DG_LOBATTO"] = _PyramidLobatto;
    list_element(2)["TETRAHEDRON_DG_LOBATTO"] = _WedgeLobatto;
    list_element(3)["TETRAHEDRON_DG_LOBATTO"] = _HexahedronLobatto;

    list_element(0)["TETRAHEDRON_DG_LOBATTO_GAUSS"] = _TetrahedronClassical;
    list_element(1)["TETRAHEDRON_DG_LOBATTO_GAUSS"] = _PyramidClassical;
    list_element(2)["TETRAHEDRON_DG_LOBATTO_GAUSS"] = _WedgeClassical;
    list_element(3)["TETRAHEDRON_DG_LOBATTO_GAUSS"] = _HexahedronGauss;
    
    list_element(0)["TETRAHEDRON_DG_CLASSICAL"] = _TetrahedronClassical;
    list_element(1)["TETRAHEDRON_DG_CLASSICAL"] = _PyramidClassical;
    list_element(2)["TETRAHEDRON_DG_CLASSICAL"] = _WedgeDgClassical;
    list_element(3)["TETRAHEDRON_DG_CLASSICAL"] = _HexahedronDgGauss;
    
    list_element(0)["HEXAHEDRON_DG_GAUSS"] = _TetrahedronClassical;
    list_element(1)["HEXAHEDRON_DG_GAUSS"] = _PyramidClassical;
    list_element(2)["HEXAHEDRON_DG_GAUSS"] = _WedgeDgClassical;
    list_element(3)["HEXAHEDRON_DG_GAUSS"] = _HexahedronDgGauss;

    list_element(0)["TETRAHEDRON_DG_LEGENDRE"] = _TetrahedronDgOrtho;
    list_element(1)["TETRAHEDRON_DG_LEGENDRE"] = _PyramidDgLegendre;
    list_element(2)["TETRAHEDRON_DG_LEGENDRE"] = _WedgeDgLegendre;
    list_element(3)["TETRAHEDRON_DG_LEGENDRE"] = _HexahedronDgLegendre;

    list_element(0)["TETRAHEDRON_DG_OPTIMAL"] = _TetrahedronClassical;
    list_element(1)["TETRAHEDRON_DG_OPTIMAL"] = _PyramidDgOrtho;
    list_element(2)["TETRAHEDRON_DG_OPTIMAL"] = _WedgeDgOrtho;
    list_element(3)["TETRAHEDRON_DG_OPTIMAL"] = _HexahedronDgGauss;

    list_element(0)["TETRAHEDRON_DG_ORTHO"] = _TetrahedronDgOrtho;
    list_element(1)["TETRAHEDRON_DG_ORTHO"] = _PyramidDgOrtho;
    list_element(2)["TETRAHEDRON_DG_ORTHO"] = _WedgeDgOrtho;
    list_element(3)["TETRAHEDRON_DG_ORTHO"] = _HexahedronDgGauss;
  }


  //! returns a new pointer to a finite element of a given id  
  ElementReference<Dimension3, 1>* VarFiniteElementEnum<Dimension3, 1>
  ::GetNewReferenceElement(int id)
  {    
    switch (id)
      {
      case _TetrahedronClassical :
	return new TetrahedronClassical();
      case _PyramidClassical :
	return new PyramidClassical();
      case _WedgeClassical :
	return new WedgeClassical();
      case _HexahedronGauss :
	return new HexahedronGauss();
      case _PyramidLobatto :
	return new PyramidLobatto();
      case _WedgeLobatto :
	return new WedgeLobatto();
      case _HexahedronLobatto :
	return new HexahedronLobatto();
#ifdef MONTJOIE_WITH_HP_H1
      case _TetrahedronHierarchic :
	return new TetrahedronHierarchic();
      case _PyramidHierarchic :
	return new PyramidHierarchic();
      case _WedgeHierarchic :
	return new WedgeHierarchic();
      case _HexahedronHierarchic :
	return new HexahedronHierarchic();
#endif
#ifdef MONTJOIE_WITH_NODAL_DG
      case _WedgeDgClassical :
	return new WedgeDgClassical();
      case _HexahedronDgGauss :
	return new HexahedronDgGauss();
#endif
#ifdef MONTJOIE_WITH_ORTHO_DG
      case _TetrahedronDgOrtho :
	return new TetrahedronDgOrtho();
      case _PyramidDgOrtho :
	return new PyramidDgOrtho();
      case _WedgeDgOrtho :
	return new WedgeDgOrtho();
#endif
#ifdef MONTJOIE_WITH_LEGENDRE_DG
      case _PyramidDgLegendre :
	return new PyramidDgLegendre();
      case _WedgeDgLegendre :
	return new WedgeDgLegendre();
      case _HexahedronDgLegendre :
	return new HexahedronDgLegendre();
#endif
     default :
	{
	  cout << "Unknown finite element" << endl;
	  abort();
	}
      }
    
    return NULL;
  }

  //! returns a new pointer to a finite element of order r and of a given id  
  ElementReference<Dimension2, 1>* VarFiniteElementEnum<Dimension3, 1>
  ::GetNewSurfaceElement(int r, int rgeom, int rquad, int id)
  {    
    switch (id)
      {
      case _TetrahedronClassical :
	{
          TriangleClassical* elt = new TriangleClassical();
          elt->ConstructFiniteElement(r, rgeom, rquad);
          return elt;
        }
#ifdef MONTJOIE_WITH_HP_H1
      case _TetrahedronHierarchic :
	{
          TriangleHierarchic* elt = new TriangleHierarchic();
          elt->SetBasisType(TriangleHierarchic::INVARIANT_BASIS);
          elt->ConstructFiniteElement(r, rgeom, rquad);
          return elt;
        }
#endif
#ifdef MONTJOIE_WITH_ORTHO_DG
      case _TetrahedronDgOrtho :
	{
          TriangleDgOrtho* elt = new TriangleDgOrtho();
          elt->ConstructFiniteElement(r, rgeom, rquad,
                                      TriangleQuadrature::QUADRATURE_GAUSS);
          return elt;
        }
#endif
#ifdef MONTJOIE_WITH_NODAL_DG
      case _PyramidClassical :
      case _PyramidHierarchic :
      case _PyramidDgOrtho :
      case _PyramidDgLegendre :
	{
          QuadrangleDgGauss* elt = new QuadrangleDgGauss();
          elt->ConstructFiniteElement(r, rgeom, rquad);
          return elt;
        }
#endif
      case _PyramidLobatto :
        {
          QuadrangleLobatto* elt = new QuadrangleLobatto();
          elt->ConstructFiniteElement(r, rgeom, rquad);
          return elt;
        }          
      default :
	{
	  cout << "Unknown finite element" << endl;
	  abort();
	}
      }

    return NULL;
  }
#endif

#endif

#ifdef MONTJOIE_WITH_NODAL_HCURL

#ifdef MONTJOIE_WITH_TWO_DIM
  //! Inits static attibutes containing the list of available finite elements
  void VarFiniteElementEnum<Dimension2, 2>::InitStaticData()
  {
    list_element.Reallocate(2);
    list_element(0)["TRIANGLE_FIRST_FAMILY"] = _TriangleHcurlFirstFamily;
    list_element(1)["TRIANGLE_FIRST_FAMILY"] = _QuadrangleHcurlFirstFamily;
    list_element(0)["QUADRANGLE_FIRST_FAMILY"] = _TriangleHcurlFirstFamily;
    list_element(1)["QUADRANGLE_FIRST_FAMILY"] = _QuadrangleHcurlFirstFamily;

    list_element(0)["QUADRANGLE_GAUSS_FIRST_FAMILY"] = _TriangleHcurlFirstFamily;
    list_element(1)["QUADRANGLE_GAUSS_FIRST_FAMILY"] = _QuadrangleHcurlGaussFirstFamily;
    
    list_element(0)["TRIANGLE_OPTIMAL_FIRST_FAMILY"] = _TriangleHcurlOptimalFirstFamily;
    list_element(1)["TRIANGLE_OPTIMAL_FIRST_FAMILY"] = _QuadrangleHcurlOptimalFirstFamily;

    list_element(0)["TRIANGLE_SECOND_FAMILY"] = _TriangleHcurlSecondFamily;
    list_element(1)["TRIANGLE_SECOND_FAMILY"] = _QuadrangleHcurlLobatto;

    list_element(0)["QUADRANGLE_HCURL_LOBATTO"] = _TriangleHcurlSecondFamily;
    list_element(1)["QUADRANGLE_HCURL_LOBATTO"] = _QuadrangleHcurlLobatto;

    list_element(0)["TRIANGLE_HP_SECOND_FAMILY"] = _TriangleHcurlHierarchic;
    
    list_element(0)["TRIANGLE_HP_FIRST_FAMILY"] = _TriangleHcurlOptimalHpFirstFamily;
    list_element(1)["TRIANGLE_HP_FIRST_FAMILY"] = _QuadrangleHcurlHpFirstFamily;
    
    list_element(0)["TRIANGLE_OPTIMAL_HP_FIRST_FAMILY"] = _TriangleHcurlOptimalHpFirstFamily;
    list_element(1)["TRIANGLE_OPTIMAL_HP_FIRST_FAMILY"] = _QuadrangleHcurlOptimalHpFirstFamily;
  }


  //! returns a new pointer to a finite element of a given id
  ElementReference<Dimension2, 2>* VarFiniteElementEnum<Dimension2, 2>
  ::GetNewReferenceElement(int id)
  {
    switch (id)
      {
      case _TriangleHcurlFirstFamily:
	return new TriangleHcurlFirstFamily();
      case _QuadrangleHcurlFirstFamily:
	return new QuadrangleHcurlFirstFamily();
      case _TriangleHcurlOptimalFirstFamily:
	return new TriangleHcurlOptimalFirstFamily();
      case _QuadrangleHcurlGaussFirstFamily:
	return new QuadrangleHcurlGaussFirstFamily();
      case _QuadrangleHcurlOptimalFirstFamily:
	return new QuadrangleHcurlOptimalFirstFamily();
      case _TriangleHcurlSecondFamily:
	return new TriangleHcurlSecondFamily();
      case _QuadrangleHcurlLobatto:
	return new QuadrangleHcurlLobatto();
#ifdef MONTJOIE_WITH_HP_HCURL
      case _TriangleHcurlHierarchic:
	return new TriangleHcurlHierarchic();
      case _QuadrangleHcurlHpFirstFamily:
	return new QuadrangleHcurlHpFirstFamily();
      case _TriangleHcurlOptimalHpFirstFamily:
	return new TriangleHcurlOptimalHpFirstFamily();
      case _QuadrangleHcurlOptimalHpFirstFamily:
	return new QuadrangleHcurlOptimalHpFirstFamily();	
#endif
      default:
	{
	  cout << "Unknown finite element" << endl;
	  abort();
	}
      }
    
    return NULL;
  }

  
  //! returns a new pointer to a finite element of order r and of a given id
  ElementReference<Dimension1, 1>* VarFiniteElementEnum<Dimension2, 2>
  ::GetNewSurfaceElement(int r, int rgeom, int rquad, int id)
  {
    switch (id)
      {
      case _TriangleHcurlFirstFamily:
	return new EdgeGaussReference();
      default:
	{
	  cout << "Unknown finite element" << endl;
	  abort();
	}
      }
    
    return NULL;
  }
#endif

#ifdef MONTJOIE_WITH_THREE_DIM
  //! Inits static attibutes containing the list of available finite elements
  void VarFiniteElementEnum<Dimension3, 2>::InitStaticData()
  {
    list_element.Reallocate(4);
    list_element(0)["TETRAHEDRON_FIRST_FAMILY"] = _TetrahedronHcurlFirstFamily;
    list_element(1)["TETRAHEDRON_FIRST_FAMILY"] = _PyramidHcurlFirstFamily;
    list_element(2)["TETRAHEDRON_FIRST_FAMILY"] = _WedgeHcurlFirstFamily;
    list_element(3)["TETRAHEDRON_FIRST_FAMILY"] = _HexahedronHcurlFirstFamily;

    list_element(0)["HEXAHEDRON_FIRST_FAMILY"] = _TetrahedronHcurlFirstFamily;
    list_element(1)["HEXAHEDRON_FIRST_FAMILY"] = _PyramidHcurlFirstFamily;
    list_element(2)["HEXAHEDRON_FIRST_FAMILY"] = _WedgeHcurlFirstFamily;
    list_element(3)["HEXAHEDRON_FIRST_FAMILY"] = _HexahedronHcurlFirstFamily;

    list_element(0)["TETRAHEDRON_OPTIMAL_FIRST_FAMILY"] = _TetrahedronHcurlOptimalFirstFamily;
    list_element(1)["TETRAHEDRON_OPTIMAL_FIRST_FAMILY"] = _PyramidHcurlOptimalFirstFamily;
    list_element(2)["TETRAHEDRON_OPTIMAL_FIRST_FAMILY"] = _WedgeHcurlOptimalFirstFamily;
    list_element(3)["TETRAHEDRON_OPTIMAL_FIRST_FAMILY"] = _HexahedronHcurlOptimalFirstFamily;

    list_element(0)["HEXAHEDRON_HCURL_LOBATTO"] = _TetrahedronHcurlLobatto;
    list_element(1)["HEXAHEDRON_HCURL_LOBATTO"] = _PyramidHcurlLobatto;
    list_element(2)["HEXAHEDRON_HCURL_LOBATTO"] = _WedgeHcurlLobatto;
    list_element(3)["HEXAHEDRON_HCURL_LOBATTO"] = _HexahedronHcurlLobatto;

    list_element(0)["TETRAHEDRON_HP_FIRST_FAMILY"] = _TetrahedronHcurlOptimalHpFirstFamily;
    list_element(1)["TETRAHEDRON_HP_FIRST_FAMILY"] = _PyramidHcurlHpFirstFamily;
    list_element(2)["TETRAHEDRON_HP_FIRST_FAMILY"] = _WedgeHcurlHpFirstFamily;
    list_element(3)["TETRAHEDRON_HP_FIRST_FAMILY"] = _HexahedronHcurlHpFirstFamily;

    list_element(0)["TETRAHEDRON_OPTIMAL_HP_FIRST_FAMILY"] = _TetrahedronHcurlOptimalHpFirstFamily;
    list_element(1)["TETRAHEDRON_OPTIMAL_HP_FIRST_FAMILY"] = _PyramidHcurlOptimalHpFirstFamily;
    list_element(2)["TETRAHEDRON_OPTIMAL_HP_FIRST_FAMILY"] = _WedgeHcurlOptimalHpFirstFamily;
    list_element(3)["TETRAHEDRON_OPTIMAL_HP_FIRST_FAMILY"] = _HexahedronHcurlOptimalHpFirstFamily;
  }


  //! returns a new pointer to a finite element of a given id
  ElementReference<Dimension3, 2>* VarFiniteElementEnum<Dimension3, 2>
  ::GetNewReferenceElement(int id)
  {
    switch (id)
      {
      case _TetrahedronHcurlFirstFamily:
	return new TetrahedronHcurlFirstFamily();
      case _PyramidHcurlFirstFamily:
	return new PyramidHcurlFirstFamily();
      case _WedgeHcurlFirstFamily:
	return new WedgeHcurlFirstFamily();
      case _HexahedronHcurlFirstFamily:
	return new HexahedronHcurlFirstFamily();
      case _TetrahedronHcurlOptimalFirstFamily:
	return new TetrahedronHcurlOptimalFirstFamily();
      case _PyramidHcurlOptimalFirstFamily:
	return new PyramidHcurlOptimalFirstFamily();
      case _WedgeHcurlOptimalFirstFamily:
	return new WedgeHcurlOptimalFirstFamily();
      case _HexahedronHcurlOptimalFirstFamily:
	return new HexahedronHcurlOptimalFirstFamily();
      case _TetrahedronHcurlLobatto:
	return new TetrahedronHcurlLobatto();
      case _HexahedronHcurlLobatto:
	return new HexahedronHcurlLobatto();
#ifdef MONTJOIE_WITH_HP_HCURL
      case _TetrahedronHcurlOptimalHpFirstFamily:
	return new TetrahedronHcurlOptimalHpFirstFamily();
      case _PyramidHcurlOptimalHpFirstFamily:
	return new PyramidHcurlOptimalHpFirstFamily();
      case _WedgeHcurlOptimalHpFirstFamily:
	return new WedgeHcurlOptimalHpFirstFamily();
      case _HexahedronHcurlOptimalHpFirstFamily:
	return new HexahedronHcurlOptimalHpFirstFamily();
      case _PyramidHcurlHpFirstFamily:
	return new PyramidHcurlHpFirstFamily();
      case _WedgeHcurlHpFirstFamily:
	return new WedgeHcurlHpFirstFamily();
      case _HexahedronHcurlHpFirstFamily:
	return new HexahedronHcurlHpFirstFamily();
#endif
      default:
	{
	  cout << "Unknown finite element, id = " << id << endl;
	  abort();
	}
      }
    
    return NULL;
  }


  //! returns a new pointer to a finite element of order r and of a given id
  ElementReference<Dimension2, 2>* VarFiniteElementEnum<Dimension3, 2>
  ::GetNewSurfaceElement(int r, int rgeom, int rquad, int id)
  {
    switch (id)
      {
      case _TetrahedronHcurlFirstFamily:
	{
	  TriangleHcurlFirstFamily* elt = new TriangleHcurlFirstFamily();
	  elt->ConstructFiniteElement(r, rgeom, rquad);
	  return elt;
	}
      case _TetrahedronHcurlLobatto:
	{
	  TriangleHcurlSecondFamily* elt = new TriangleHcurlSecondFamily();
	  elt->ConstructFiniteElement(r, rgeom, rquad);
	  return elt;
	}	
      case _PyramidHcurlFirstFamily :
	{
          QuadrangleHcurlFirstFamily* elt = new QuadrangleHcurlFirstFamily();
          elt->ConstructFiniteElement(r, rgeom, rquad);
          return elt;
        }
      case _PyramidHcurlLobatto :
	{
	  QuadrangleHcurlLobatto* elt = new QuadrangleHcurlLobatto();
          elt->ConstructFiniteElement(r, rgeom, rquad);
          return elt;
	}
      default:
	{
	  cout << "Unknown finite element, id = " << id << endl;
	  abort();
	}
      }
    
    return NULL;
  }
#endif

#endif


#ifdef MONTJOIE_WITH_NODAL_HDIV

#ifdef MONTJOIE_WITH_TWO_DIM
  //! Inits static attibutes containing the list of available finite elements
  void VarFiniteElementEnum<Dimension2, 3>::InitStaticData()
  {
    list_element.Reallocate(2);
    list_element(0)["TRIANGLE_FIRST_FAMILY"] = _TriangleHdivFirstFamily;
    list_element(1)["TRIANGLE_FIRST_FAMILY"] = _QuadrangleHdivFirstFamily;
    list_element(0)["QUADRANGLE_FIRST_FAMILY"] = _TriangleHdivFirstFamily;
    list_element(1)["QUADRANGLE_FIRST_FAMILY"] = _QuadrangleHdivFirstFamily;

    list_element(0)["TRIANGLE_OPTIMAL_FIRST_FAMILY"] = _TriangleHdivOptimalFirstFamily;
    list_element(1)["TRIANGLE_OPTIMAL_FIRST_FAMILY"] = _QuadrangleHdivOptimalFirstFamily;
    
    list_element(0)["TRIANGLE_HP_FIRST_FAMILY"] = _TriangleHdivOptimalHpFirstFamily;
    list_element(1)["TRIANGLE_HP_FIRST_FAMILY"] = _QuadrangleHdivHpFirstFamily;
    
    list_element(0)["TRIANGLE_OPTIMAL_HP_FIRST_FAMILY"] = _TriangleHdivOptimalHpFirstFamily;
    list_element(1)["TRIANGLE_OPTIMAL_HP_FIRST_FAMILY"] = _QuadrangleHdivOptimalHpFirstFamily;
  }


  //! returns a new pointer to a finite element of a given id
  ElementReference<Dimension2, 3>* VarFiniteElementEnum<Dimension2, 3>
  ::GetNewReferenceElement(int id)
  {
    switch (id)
      {
      case _TriangleHdivFirstFamily:
	return new TriangleHdivFirstFamily();
      case _QuadrangleHdivFirstFamily:
	return new QuadrangleHdivFirstFamily();
      case _TriangleHdivOptimalFirstFamily:
	return new TriangleHdivOptimalFirstFamily();
      case _QuadrangleHdivOptimalFirstFamily:
	return new QuadrangleHdivOptimalFirstFamily();
#ifdef MONTJOIE_WITH_HP_HDIV
      case _QuadrangleHdivHpFirstFamily:
	return new QuadrangleHdivHpFirstFamily();
      case _TriangleHdivOptimalHpFirstFamily:
	return new TriangleHdivOptimalHpFirstFamily();
      case _QuadrangleHdivOptimalHpFirstFamily:
	return new QuadrangleHdivOptimalHpFirstFamily();	
#endif
      default:
	{
	  cout << "Unknown finite element" << endl;
	  abort();
	}
      }
    
    return NULL;
  }


  //! returns a new pointer to a finite element of order r and of a given id
  ElementReference<Dimension1, 1>* VarFiniteElementEnum<Dimension2, 3>
  ::GetNewSurfaceElement(int r, int rgeom, int rquad, int id)
  {
    switch (id)
      {
      case _TriangleHdivFirstFamily:
	return new EdgeGaussReference();
      default:
	{
	  cout << "Unknown finite element" << endl;
	  abort();
	}
      }
    
    return NULL;
  }
#endif

#ifdef MONTJOIE_WITH_THREE_DIM
  //! Inits static attibutes containing the list of available finite elements
  void VarFiniteElementEnum<Dimension3, 3>::InitStaticData()
  {
    list_element.Reallocate(4);
    list_element(0)["TETRAHEDRON_FIRST_FAMILY"] = _TetrahedronHdivFirstFamily;
    list_element(1)["TETRAHEDRON_FIRST_FAMILY"] = _PyramidHdivFirstFamily;
    list_element(2)["TETRAHEDRON_FIRST_FAMILY"] = _WedgeHdivFirstFamily;
    list_element(3)["TETRAHEDRON_FIRST_FAMILY"] = _HexahedronHdivFirstFamily;

    list_element(0)["TETRAHEDRON_OPTIMAL_FIRST_FAMILY"] = _TetrahedronHdivOptimalFirstFamily;
    list_element(1)["TETRAHEDRON_OPTIMAL_FIRST_FAMILY"] = _PyramidHdivOptimalFirstFamily;
    list_element(2)["TETRAHEDRON_OPTIMAL_FIRST_FAMILY"] = _WedgeHdivOptimalFirstFamily;
    list_element(3)["TETRAHEDRON_OPTIMAL_FIRST_FAMILY"] = _HexahedronHdivOptimalFirstFamily;

    list_element(0)["TETRAHEDRON_HP_FIRST_FAMILY"] = _TetrahedronHdivOptimalHpFirstFamily;
    list_element(1)["TETRAHEDRON_HP_FIRST_FAMILY"] = _PyramidHdivHpFirstFamily;
    list_element(2)["TETRAHEDRON_HP_FIRST_FAMILY"] = _WedgeHdivHpFirstFamily;
    list_element(3)["TETRAHEDRON_HP_FIRST_FAMILY"] = _HexahedronHdivHpFirstFamily;

    list_element(0)["TETRAHEDRON_OPTIMAL_HP_FIRST_FAMILY"] = _TetrahedronHdivOptimalHpFirstFamily;
    list_element(1)["TETRAHEDRON_OPTIMAL_HP_FIRST_FAMILY"] = _PyramidHdivOptimalHpFirstFamily;
    list_element(2)["TETRAHEDRON_OPTIMAL_HP_FIRST_FAMILY"] = _WedgeHdivOptimalHpFirstFamily;
    list_element(3)["TETRAHEDRON_OPTIMAL_HP_FIRST_FAMILY"] = _HexahedronHdivOptimalHpFirstFamily;
  }


  //! returns a new pointer to a finite element of a given id
  ElementReference<Dimension3, 3>* VarFiniteElementEnum<Dimension3, 3>
  ::GetNewReferenceElement(int id)
  {
    switch (id)
      {
      case _TetrahedronHdivFirstFamily:
	return new TetrahedronHdivFirstFamily();
      case _PyramidHdivFirstFamily:
	return new PyramidHdivFirstFamily();
      case _WedgeHdivFirstFamily:
	return new WedgeHdivFirstFamily();
      case _HexahedronHdivFirstFamily:
	return new HexahedronHdivFirstFamily();
      case _TetrahedronHdivOptimalFirstFamily:
	return new TetrahedronHdivOptimalFirstFamily();
      case _PyramidHdivOptimalFirstFamily:
	return new PyramidHdivOptimalFirstFamily();
      case _WedgeHdivOptimalFirstFamily:
	return new WedgeHdivOptimalFirstFamily();
      case _HexahedronHdivOptimalFirstFamily:
	return new HexahedronHdivOptimalFirstFamily();
#ifdef MONTJOIE_WITH_HP_HDIV
      case _TetrahedronHdivOptimalHpFirstFamily:
	return new TetrahedronHdivOptimalHpFirstFamily();
      case _PyramidHdivOptimalHpFirstFamily:
	return new PyramidHdivOptimalHpFirstFamily();
      case _WedgeHdivOptimalHpFirstFamily:
	return new WedgeHdivOptimalHpFirstFamily();
      case _HexahedronHdivOptimalHpFirstFamily:
	return new HexahedronHdivOptimalHpFirstFamily();
      case _PyramidHdivHpFirstFamily:
	return new PyramidHdivHpFirstFamily();
      case _WedgeHdivHpFirstFamily:
	return new WedgeHdivHpFirstFamily();
      case _HexahedronHdivHpFirstFamily:
	return new HexahedronHdivHpFirstFamily();
#endif
      default:
	{
	  cout << "Unknown finite element" << endl;
	  abort();
	}
      }
    
    return NULL;
  }


  //! returns a new pointer to a finite element of order r and of a given id
  ElementReference<Dimension2, 1>* VarFiniteElementEnum<Dimension3, 3>
  ::GetNewSurfaceElement(int r, int rgeom, int rquad, int id)
  {
    switch (id)
      {
      case _TetrahedronHdivFirstFamily:
	return new TriangleClassical();
      default:
	{
	  cout << "Unknown finite element" << endl;
	  abort();
	}
      }
    
    return NULL;
  }
#endif

#endif
  
  
  /*************************
   * VarFiniteElement_Base *
   *************************/

  //! default constructor
  VarFiniteElement_Base::VarFiniteElement_Base()
  {
    use_piola_transform = true;
    nb_elt_ = 0; nb_surf_ = 0; nb_surf_neigh_ = 0;
  }


  void VarFiniteElement_Base::SetInputData(const string& description_field, const VectString& parameters)
  {
    if (description_field == "UsePiolaTransform")
      {
	if (parameters(0) == "YES")
	  this->use_piola_transform = true;
	else
	  this->use_piola_transform = false;
      }
  }
  
  
#ifdef MONTJOIE_WITH_ONE_DIM
  //! Constructs finite element and affects them to each element
  void VarFiniteElement_Base
  ::AddFiniteElement1D(const string& name_elt, int order,
		       const Vector<bool>& change_elt, Mesh<Dimension1>& mesh)
  {
    VarFiniteElement<Dimension1>& var = static_cast<VarFiniteElement<Dimension1>& >(*this);
    
    // we assume that name_elt is new (and not partially contained in stored_reference_element)
    // resizing the array
    int nb_new_elt = 1;
    int nb_old_elt = var.stored_reference_element.GetM();
    var.stored_reference_element.Resize(nb_old_elt + nb_new_elt);
    var.id_stored_reference_element.Resize(nb_old_elt + nb_new_elt);
    
    // constructing the new elements
    int id_elt = VarFiniteElementEnum<Dimension1, 1>::GetIdentityNumber(name_elt, 0);
    nb_new_elt = nb_old_elt;
    var.id_stored_reference_element(nb_new_elt) = id_elt;
    var.stored_reference_element(nb_new_elt)
      = VarFiniteElementEnum<Dimension1, 1>::GetNewReferenceElement(id_elt);
    
    var.stored_reference_element(nb_new_elt)->ConstructFiniteElement(order);
    
    // filling reference_element
    this->nb_elt_ = mesh.GetNbElt();
    var.reference_element.Reallocate(mesh.GetNbElt());
    for (int i = 0; i < mesh.GetNbElt(); i++)
      if (change_elt(i))
	var.reference_element(i) = var.stored_reference_element(nb_new_elt);
  }
#endif


#ifdef MONTJOIE_WITH_TWO_DIM
  //! retrieves the type of integration on boundaries and checks compatibility
  void VarFiniteElement_Base
  ::GetTypeIntegrationBoundary(const Mesh<Dimension2>& mesh, const MeshNumbering<Dimension2>& mesh_num,
                               const Vector<MeshNumbering<Dimension2>* >& other_mesh_num,
			       int& type_integration_edge,
			       int& type_integration_tri, int& type_integration_quad)
  {
    VarFiniteElement<Dimension2>& var = static_cast<VarFiniteElement<Dimension2>& >(*this);
    type_integration_edge = -1;
    type_integration_tri = -1;
    type_integration_quad = -1;
    int N = other_mesh_num.GetM() + 1;
    for (int i = 0; i < mesh.GetNbElt(); i++)
      for (int n = 0; n < N; n++)
        {
          int t = var.GetReferenceElement(i, n).GetTypeIntegrationEdge();
          if (type_integration_edge == -1)
            type_integration_edge = t;
          else if (type_integration_edge != t)
            {
              cout << "Use the same quadrature rule for edges" << endl;
              abort();
            }
        }
  }
#endif


#ifdef MONTJOIE_WITH_THREE_DIM
  //! retrieves the type of integration on boundaries and checks compatibility
  void VarFiniteElement_Base
  ::GetTypeIntegrationBoundary(const Mesh<Dimension3>& mesh, const MeshNumbering<Dimension3>& mesh_num,
                               const Vector<MeshNumbering<Dimension3>* >& other_mesh_num,
			       int& type_integration_edge,
			       int& type_integration_tri, int& type_integration_quad)
  {
    VarFiniteElement<Dimension3>& var = static_cast<VarFiniteElement<Dimension3>& >(*this);
    type_integration_edge = -1;
    type_integration_tri = -1;
    type_integration_quad = -1;
    int N = other_mesh_num.GetM() + 1;
    for (int i = 0; i < mesh.GetNbElt(); i++)
      for (int n = 0; n < N; n++)
        {
          int tri = var.GetReferenceElement(i, n).GetTypeIntegrationTriangle();
          int quad = var.GetReferenceElement(i, n).GetTypeIntegrationQuadrangle();
          int t = mesh.GetTypeElement(i);
          if (t < 3)
            {
              // element with a triangular face
              if (type_integration_tri == -1)
                type_integration_tri = tri;
              else if (type_integration_tri != tri)
                {
                  cout << "Use the same quadrature rule for triangles" << endl;
                  abort();
                }
            }
          
          if (t > 0)
            {
              // element with a quadrilateral face
              if (type_integration_quad == -1)
                type_integration_quad = quad;
              else if (type_integration_quad != quad)
                {
                  cout << "Use the same quadrature rule for quadrangles" << endl;
                  abort();
                }
            }
        }
  }
#endif
  
  
  /********************
   * VarFiniteElement *
   ********************/
  
  
  //! destructor
  template<class Dimension>
  VarFiniteElement<Dimension>::~VarFiniteElement()
  {
    ClearFiniteElement();
  }
  
  
  //! releases memory
  template<class Dimension>
  void VarFiniteElement<Dimension>::ClearFiniteElement()
  {
    for (int i = 0; i < this->stored_reference_element.GetM(); i++)
      {
	if (stored_reference_element(i) != NULL)
	  {
	    delete stored_reference_element(i);
	    stored_reference_element(i) = NULL;
	  }
      }

    type_stored_reference_element.Clear();
    id_stored_reference_element.Clear();
    stored_reference_element.Clear();
    reference_element.Clear();
    this->nb_elt_ = 0; this->nb_surf_ = 0; this->nb_surf_neigh_ = 0;
  }
  
  
  //! adds memory used by the current object
  template<class Dimension>
  void VarFiniteElement<Dimension>::GetMemoryUsed(map<string, size_t>& var) const
  {
    size_t n = sizeof(Vector<int>)*stored_reference_element.GetM();
    n += reference_element.GetMemorySize();
    n += reference_neighbor_element.GetMemorySize();
    n += id_stored_reference_element.GetMemorySize();
    for (int i = 0; i < stored_reference_element.GetM(); i++)
      n += stored_reference_element(i)->GetMemorySize();
    
    var["FiniteElement"] = n;
  }


  //! returns the identity number of an element (according with its name and type)
  /*!
    name : name of the finite element
    t : type of element (0, 1, 2 or 3 for tet/pyramid/prism/hex)
    type : nature of element (1: H1, 2 : Hcurl, 3 : Hdiv)
   */
  template<class Dimension>
  int VarFiniteElement<Dimension>::GetIdentityNumber(const string& name, int t, int type)
  {
    switch (type)
      {
      case 1:
        return VarFiniteElementEnum<Dimension, 1>::GetIdentityNumber(name, t);
      case 2:
        return VarFiniteElementEnum<Dimension, 2>::GetIdentityNumber(name, t);
      case 3:
        return VarFiniteElementEnum<Dimension, 3>::GetIdentityNumber(name, t);
      default:
        cout << "Type must be equal to 1, 2 or 3" << endl;
        abort();
      }
  }
  
    
  //! returns the pointer of the surface finite element if already stored in the class
  template<class Dimension>
  ElementReference_Dim<typename Dimension::DimensionBoundary>* VarFiniteElement<Dimension>
  ::GetSurfaceFiniteElementPtr(int r, int id_elt, int type_elt) const
  {
    ElementReferenceTrace* null_ptr;
    null_ptr = NULL; 
    
    for (int i = 0; i < id_stored_surface_element.GetM(); i++)
      if (id_stored_surface_element(i) == id_elt)
        {
          if (stored_surface_element(i)->GetOrder() == r)
            if (type_stored_surface_element(i) == type_elt)
              return stored_surface_element(i);
        }

    return null_ptr;
  }


  //! returns the pointer of the reference finite element if already stored in the class
  template<class Dimension>
  ElementReference_Dim<Dimension>* VarFiniteElement<Dimension>
  ::GetReferenceFiniteElementPtr(int r, int id_elt, int type_elt) const
  {
    ElementReference_Dim<Dimension>* null_ptr;
    null_ptr = NULL; 
    
    for (int i = 0; i < id_stored_reference_element.GetM(); i++)
      if (id_stored_reference_element(i) == id_elt)
        {
          if (stored_reference_element(i)->GetOrder() == r)
            if (type_stored_reference_element(i) == type_elt)
              return stored_reference_element(i);
        }

    return null_ptr;
  }


  //! no surface element in 1-D
  template<>
  ElementReference_Dim<Dimension1>*
  VarFiniteElement<Dimension1>::GetNewSurfaceElement(int r, int rgeom, int rquad, int id_elt, int type_elt)
  {
    return NULL;
  }
  

  //! create a new surface element
  template<class Dimension>
  ElementReference_Dim<typename Dimension::DimensionBoundary>*
  VarFiniteElement<Dimension>
  ::GetNewSurfaceElement(int r, int rgeom, int rquad, int id_elt, int type_elt)
  {
    switch (type_elt)
      {
      case 1:
        return VarFiniteElementEnum<Dimension, 1>::GetNewSurfaceElement(r, rgeom, rquad, id_elt);
      case 2:
        return VarFiniteElementEnum<Dimension, 2>::GetNewSurfaceElement(r, rgeom, rquad, id_elt);
      case 3:
        return VarFiniteElementEnum<Dimension, 3>::GetNewSurfaceElement(r, rgeom, rquad, id_elt);
      default:
        cout << "Type must be equal to 1, 2 or 3" << endl;
        abort();
      }
  }
  

  //! create a new reference element
  template<>
  ElementReference_Dim<Dimension1>* VarFiniteElement<Dimension1>
  ::GetNewReferenceElement(int id_elt, int type_elt)
  {
    return VarFiniteElementEnum<Dimension1, 1>::GetNewReferenceElement(id_elt);
  }
  

  //! create a new reference element
  template<class Dimension>
  ElementReference_Dim<Dimension>* VarFiniteElement<Dimension>
  ::GetNewReferenceElement(int id_elt, int type_elt)
  {
    switch (type_elt)
      {
      case 1:
        return VarFiniteElementEnum<Dimension, 1>::GetNewReferenceElement(id_elt);
      case 2: 
        return VarFiniteElementEnum<Dimension, 2>::GetNewReferenceElement(id_elt);
      case 3: 
        return VarFiniteElementEnum<Dimension, 3>::GetNewReferenceElement(id_elt);
      default:
        cout << "Type must be equal to 1, 2 or 3" << endl;
        abort();
      } 
  }


  //! constructs the finite element
  template<>
  void VarFiniteElement<Dimension2>
  ::ConstructFiniteElement(int r, int rgeom, int rquad, int type_elt,
                           ElementReference_Dim<Dimension2>* elt)
  {
    switch (type_elt)
      {
      case 1 :
        static_cast<FaceReference<1>& >(*elt).
          ConstructFiniteElement(r, rgeom, rquad);
        break;
      case 2 :
        static_cast<FaceReference<2>& >(*elt).
          ConstructFiniteElement(r, rgeom, rquad);
        break;
      case 3 :
        static_cast<FaceReference<3>& >(*elt).
          ConstructFiniteElement(r, rgeom, rquad);
        break;
      }
  }


  //! constructs the finite element
  template<>
  void VarFiniteElement<Dimension3>
  ::ConstructFiniteElement(int r, int rgeom, int rquad, int type_elt,
                           ElementReference_Dim<Dimension3>* elt)
  {
    switch (type_elt)
      {
      case 1 :
        static_cast<VolumeReference<1>& >(*elt).
          ConstructFiniteElement(r, rgeom, rquad);
        break;
      case 2 :
        static_cast<VolumeReference<2>& >(*elt).
          ConstructFiniteElement(r, rgeom, rquad);
        break;
      case 3 :
        static_cast<VolumeReference<3>& >(*elt).
          ConstructFiniteElement(r, rgeom, rquad);
        break;
      }
  }

  
  //! Constructs surface finite element and affects them to each face
  template<class Dimension>
  void VarFiniteElement<Dimension>
  ::AddSurfaceFiniteElement(const string& name_elt, TinyVector<IVect, 4>& order,
			    Vector<bool>& change_elt, Mesh<Dimension>& mesh,
			    MeshNumbering<Dimension>& mesh_num, int r_over_quad, 
			    int discontinuous_formulation, int type)
  {
    Vector<int> other_type; Vector<TinyVector<IVect, 4> > other_order;
    Vector<string> other_name; Vector<int> r_over_q(1); r_over_q(0) = r_over_quad;
    Vector<MeshNumbering<Dimension>* > other_mesh_num;
    AddSurfaceFiniteElement(name_elt, order, type, other_name, other_order, other_type,
                            change_elt, mesh, mesh_num, other_mesh_num, r_over_q,
                            discontinuous_formulation);
  }
  

  //! no surface element in 1-D
  template<>
  void VarFiniteElement<Dimension1>
  ::AddSurfaceFiniteElement(const string& name_elt, TinyVector<IVect, 4>& order, int type,
                            const Vector<string>& other_name,
                            Vector<TinyVector<IVect, 4> >& other_order, Vector<int>& other_type,
			    Vector<bool>& change_elt, Mesh<Dimension1>& mesh,
			    MeshNumbering<Dimension1>& mesh_num,
                            Vector<MeshNumbering<Dimension1>* >& other_mesh_num,
                            Vector<int>& r_over_quad, int discontinuous_formulation)
  {
    abort();
  }
  

  //! Constructs surface finite element and affects them to each face
  template<class Dimension>
  void VarFiniteElement<Dimension>
  ::AddSurfaceFiniteElement(const string& name_elt, TinyVector<IVect, 4>& order, int type,
                            const Vector<string>& other_name,
                            Vector<TinyVector<IVect, 4> >& other_order, Vector<int>& other_type,
			    Vector<bool>& change_elt, Mesh<Dimension>& mesh,
			    MeshNumbering<Dimension>& mesh_num,
                            Vector<MeshNumbering<Dimension>* >& other_mesh_num,
                            Vector<int>& r_over_quad, int discontinuous_formulation)
  {
    for (int t = 0; t < 2; t++)
      RemoveDuplicate(order(t));

    for (int n = 0; n < other_order.GetM(); n++)
      for (int t = 0; t < 2; t++)
        RemoveDuplicate(other_order(n)(t));
    
    // we compute the maximal order for each unknown
    ElementReferenceTrace* null_ptr;
    null_ptr = NULL; 
    int N = other_order.GetM() + 1;
    Vector<Array2D<ElementReferenceTrace* > > ptr_elt(N);
    for (int n = 0; n < N; n++)
      {
        int rmax = 0;
        if (n == 0)
          for (int t = 0; t < 2; t++)
            for (int k = 0; k < order(t).GetM(); k++)
              rmax = max(rmax, order(t)(k));
        else
          for (int t = 0; t < 2; t++)
            for (int k = 0; k < other_order(n-1)(t).GetM(); k++)
              rmax = max(rmax, other_order(n-1)(t)(k));

        // initializing the array
        ptr_elt(n).Reallocate(2, rmax+1);
        ptr_elt(n).Fill(null_ptr);
      }
    
    // constructing the new elements
    int rgeom = mesh.GetGeometryOrder();
    for (int n = 0; n < N; n++)
      for (int t = 0; t < 2; t++)
        {
          int Nr = 0;
          if (n == 0)
            Nr = order(t).GetM();
          else
            Nr = other_order(n-1)(t).GetM();
          
          if (Nr > 0)
            {
              int id_elt = 0, type_elt = 0;
              if (n == 0)
                {
                  type_elt = type;
                  id_elt = this->GetIdentityNumber(name_elt, t, type);
                }
              else
                {
                  type_elt = other_type(n-1);
                  id_elt = this->GetIdentityNumber(other_name(n-1), t, other_type(n-1));
                }
              
              for (int k = 0; k < Nr; k++)
                {
                  int r = 0;
                  if (n == 0)
                    r = order(t)(k);
                  else
                    r = other_order(n-1)(t)(k);
                  
                  int r1 = max(r, 1);
                  int rquad = r + r_over_quad(n);

                  ptr_elt(n)(t, r) = this->GetSurfaceFiniteElementPtr(r, id_elt, type_elt);
                  if (ptr_elt(n)(t, r) == null_ptr)
                    {
                      id_stored_surface_element.PushBack(id_elt);
                      ElementReferenceTrace* elt =
                        this->GetNewSurfaceElement(r, min(r1, rgeom), rquad, id_elt, type_elt);
                      stored_surface_element.PushBack(elt);
                      type_stored_surface_element.PushBack(type_elt);
                      ptr_elt(n)(t, r) = elt;
                    }
                }
	    }
	}
    
    // filling reference_element
    this->nb_surf_ = mesh.GetNbBoundary();
    surface_element.Reallocate(N*mesh.GetNbBoundary());
    for (int i = 0; i < mesh.GetNbBoundary(); i++)
      if (change_elt(i))
        for (int n = 0; n < N; n++)
          {
            int t = mesh.Boundary(i).GetHybridType();
            int r = mesh_num.GetOrderQuadrature(i);
            if (n > 0)
              r = other_mesh_num(n-1)->GetOrderQuadrature(i);
            
            surface_element(i + n*this->nb_surf_) = ptr_elt(n)(t, r);
            if (ptr_elt(n)(t, r) == NULL)
              {
                cout << "This case should not appear" << endl;
                abort();
              }
          }
  }


  //! Constructs finite element and affects them to each element
  template<class Dimension>
  void VarFiniteElement<Dimension>
  ::AddFiniteElement(const string& name_elt, TinyVector<IVect, 4>& order,
		     Vector<bool>& change_elt, Mesh<Dimension>& mesh,
		     MeshNumbering<Dimension>& mesh_num, int r_over_q, 
		     int discontinuous_formulation, int type)
  {
    Vector<int> other_type; Vector<TinyVector<IVect, 4> > other_order;
    Vector<string> other_name; Vector<int> r_over_quad(1); r_over_quad(0) = r_over_q;
    Vector<MeshNumbering<Dimension>* > other_mesh_num;
    Vector<int> other_discont_form;
    this->AddFiniteElement(name_elt, order, type, other_name, other_order, other_type,
                           change_elt, mesh, mesh_num,
                           other_mesh_num, r_over_quad, discontinuous_formulation,
                           other_discont_form);
  }


  template<>
  void VarFiniteElement<Dimension1>
  ::AddFiniteElement(const string& name_elt, TinyVector<IVect, 4>& order, int type,
                     const Vector<string>& other_name,
                     Vector<TinyVector<IVect, 4> >& other_order, const Vector<int>& other_type,
		     Vector<bool>& change_elt, Mesh<Dimension1>& mesh,
		     MeshNumbering<Dimension1>& mesh_num,
                     Vector<MeshNumbering<Dimension1>* >& other_mesh_num,
                     Vector<int>& r_over_quad, int discontinuous_formulation,
                     IVect& other_discont_form)
  {
    AddFiniteElement1D(name_elt, order(0)(0), change_elt, mesh);
  }

  
  //! Constructs finite element and affects them to each element
  template<class Dimension>
  void VarFiniteElement<Dimension>
  ::AddFiniteElement(const string& name_elt, TinyVector<IVect, 4>& order, int type,
                     const Vector<string>& other_name,
                     Vector<TinyVector<IVect, 4> >& other_order, const Vector<int>& other_type,
		     Vector<bool>& change_elt, Mesh<Dimension>& mesh,
		     MeshNumbering<Dimension>& mesh_num,
                     Vector<MeshNumbering<Dimension>* >& other_mesh_num,
                     Vector<int>& r_over_quad, int discontinuous_formulation,
                     IVect& other_discont_form)
  {
    for (int t = 0; t < 4; t++)
      RemoveDuplicate(order(t));

    for (int n = 0; n < other_order.GetM(); n++)
      for (int t = 0; t < 4; t++)
        RemoveDuplicate(other_order(n)(t));

    // we compute the maximal order for each unknown
    ElementReference_Dim<Dimension>* null_ptr;
    null_ptr = NULL; 
    int N = other_order.GetM() + 1;
    Vector<Array2D<ElementReference_Dim<Dimension>* > > ptr_elt(N);
    for (int n = 0; n < N; n++)
      {
        int rmax = 0;
        if (n == 0)
          for (int t = 0; t < 4; t++)
            for (int k = 0; k < order(t).GetM(); k++)
              rmax = max(rmax, order(t)(k));
        else
          for (int t = 0; t < 4; t++)
            for (int k = 0; k < other_order(n-1)(t).GetM(); k++)
              rmax = max(rmax, other_order(n-1)(t)(k));
        
        // initializing the array
        ptr_elt(n).Reallocate(4, rmax+1);
        ptr_elt(n).Fill(null_ptr);
      }

    // constructing the new elements
    int rgeom = mesh.GetGeometryOrder();
    for (int n = 0; n < N; n++)
      for (int t = 0; t < 4; t++)
        {
          int Nr = 0;
          if (n == 0)
            Nr = order(t).GetM();
          else
            Nr = other_order(n-1)(t).GetM();

          int dg_form = discontinuous_formulation;
          if (n > 0)
            dg_form = other_discont_form(n-1);
          
          if (Nr > 0)
            {
              int id_elt = 0, type_elt = 0;
              if (n == 0)
                {
                  type_elt = type;
                  id_elt = this->GetIdentityNumber(name_elt, t, type);
                }
              else
                {
                  type_elt = other_type(n-1);
                  id_elt = this->GetIdentityNumber(other_name(n-1), t, other_type(n-1));
                }

              for (int k = 0; k < Nr; k++)
                {
                  int r = 0;
                  if (n == 0)
                    r = order(t)(k);
                  else
                    r = other_order(n-1)(t)(k);

                  int r1 = max(r, 1);
                  int rquad = r + r_over_quad(n);

                  ElementReference_Dim<Dimension>* elt;
                  elt = this->GetReferenceFiniteElementPtr(r, id_elt, type_elt);
                  if (elt == null_ptr)
                    {
                      id_stored_reference_element.PushBack(id_elt);
                      elt = this->GetNewReferenceElement(id_elt, type_elt);
                      
                      stored_reference_element.PushBack(elt);
                      type_stored_reference_element.PushBack(type_elt);
                      ptr_elt(n)(t, r) = elt;

                      elt->SetPiolaTransform(use_piola_transform);
                      this->ConstructFiniteElement(r, min(r1, rgeom), rquad, type_elt, elt);
                      elt->SetMesh(mesh);
                      elt->SetVariableOrder(mesh, mesh_num);
	      
                      if (dg_form == ElementReference_Base::CONTINUOUS)
                        if (elt-> DiscontinuousElement())
                          {
                            cout << "A discontinuous element cannot be used with a continuous formulation" << endl;
                            abort();
                          }
                    }
                  else
                    ptr_elt(n)(t, r) = elt;

                  if (n == 0)
                    elt->ConstructNumberMap(mesh_num.number_map, dg_form);
                  else
                    elt->ConstructNumberMap(other_mesh_num(n-1)->number_map, dg_form);
                }
            }
        }
    
    // filling reference_element
    this->nb_elt_ = mesh.GetNbElt();
    reference_element.Reallocate(mesh.GetNbElt()*N);
    for (int i = 0; i < mesh.GetNbElt(); i++)
      if (change_elt(i))
        for (int n = 0; n < N; n++)
          {
            int t = mesh.GetTypeElement(i);
            int r = mesh_num.GetOrderElement(i);
            if (n > 0)
              r = other_mesh_num(n-1)->GetOrderElement(i);

            reference_element(i + n*this->nb_elt_) = ptr_elt(n)(t, r);
            if (ptr_elt(n)(t, r) == NULL)
              {
                cout << "This case should not appear" << endl;
                abort();
              }
          }
    
    // filling reference_neighbor_element
    this->nb_surf_neigh_ = mesh_num.GetNbNeighborElt();
    reference_neighbor_element.Reallocate(mesh_num.GetNbNeighborElt()*N);
    for (int i = 0; i < mesh_num.GetNbNeighborElt(); i++)
      for (int n = 0; n < N; n++)
        {
          int t(0), r(0);
          if (n == 0)
            {
              t = mesh_num.GetTypeNeighborElement(i);
              r = mesh_num.GetOrderNeighborElement(i);
            }
          else
            {
              t = other_mesh_num(n-1)->GetTypeNeighborElement(i);
              r = other_mesh_num(n-1)->GetOrderNeighborElement(i);
            }
          
          reference_neighbor_element(i + n*this->nb_surf_neigh_) = ptr_elt(n)(t, r);
          if (ptr_elt(n)(t, r) == NULL)
            {
              cout << "This case should not appear" << endl;
              abort();
            }        
        }
  }

  template<>
  void VarFiniteElement<Dimension1>
  ::UpdateInterpolationElement(Mesh<Dimension1>& mesh, MeshNumbering<Dimension1>& mesh_num)
  {
  }
  
  
  //! Updates finite elements (when order of quadrature are changed on faces)
  template<class Dimension>
  void VarFiniteElement<Dimension>
  ::UpdateInterpolationElement(Mesh<Dimension>& mesh, MeshNumbering<Dimension>& mesh_num)
  {
    for (int i = 0; i < stored_reference_element.GetM(); i++)
      stored_reference_element(i)->SetVariableOrder(mesh, mesh_num);  
  }
  

  template<>
  void VarFiniteElement<Dimension1>
  ::CopyFiniteElement(const VarFiniteElement<Dimension1>& var,
		      const Mesh<Dimension1>& mesh,
                      const MeshNumbering<Dimension1>& mesh_num,
                      const Vector<MeshNumbering<Dimension1>* >& other_mesh_num,
                      int type, const Vector<int>& other_type)
  {
    abort();
  }
  

  //! copies finite element from another object
  template<class Dimension>
  void VarFiniteElement<Dimension>
  ::CopyFiniteElement(const VarFiniteElement<Dimension>& var,
		      const Mesh<Dimension>& mesh,
                      const MeshNumbering<Dimension>& mesh_num,
                      const Vector<MeshNumbering<Dimension>* >& other_mesh_num,
                      int type, const Vector<int>& other_type)
  {
    int N = other_mesh_num.GetM() + 1;
    // volumic finite elements
    {
      int rmax = 1;
      for (int k = 0; k < var.stored_reference_element.GetM(); k++)
        rmax = max(rmax, var.stored_reference_element(k)->GetOrder());
      
      ElementReference_Dim<Dimension>* null_ptr;
      null_ptr = NULL;
      TinyMatrix<Vector<ElementReference_Dim<Dimension>* >, General, 4, 4> elt_target;
      for (int t = 0; t < 4; t++)
        for (int p = 0; p < 4; p++)
          {
            elt_target(t, p).Reallocate(rmax+1);
            elt_target(t, p).Fill(null_ptr);
          }
      
      for (int p = 0; p < var.stored_reference_element.GetM(); p++)
        {
          int t = var.stored_reference_element(p)->GetHybridType();
          int r = var.stored_reference_element(p)->GetOrder();
          int q = var.type_stored_reference_element(p);
          elt_target(t, q)(r) = var.stored_reference_element(p);
        }
      
      // shallow copy : only pointers are copied
      this->nb_elt_ = mesh.GetNbElt();
      reference_element.Reallocate(mesh.GetNbElt()*N);
      for (int n = 0; n < N; n++)
        for (int i = 0; i < mesh.GetNbElt(); i++)
          {
            int t = mesh.GetTypeElement(i);
            int r = 0, q = 0;
            if (n == 0)
              {
                r = mesh_num.GetOrderElement(i);
                q = type;
              }
            else
              {
                r = other_mesh_num(n-1)->GetOrderElement(i);
                q = other_type(n-1);
              }
            
            reference_element(i + n*this->nb_elt_) = elt_target(t, q)(r);
            if (elt_target(t, q)(r) == NULL)
              {
                cout << "This case should not appear" << endl;
                abort();
              }
          }
    }
    
    // surface finite elements
    if (var.stored_surface_element.GetM() > 0)
      {
	int rmax = 1;
	for (int k = 0; k < var.stored_surface_element.GetM(); k++)
	  rmax = max(rmax, var.stored_surface_element(k)->GetOrder());
	
	ElementReferenceTrace* null_ptr;
	null_ptr = NULL;
	TinyMatrix<Vector<ElementReferenceTrace* >, General, 2, 4> elt_target;
	for (int t = 0; t < 2; t++)
          for (int p = 0; p < 4; p++)
            {
              elt_target(t, p).Reallocate(rmax+1);
              elt_target(t, p).Fill(null_ptr);
            }
	
	for (int p = 0; p < var.stored_surface_element.GetM(); p++)
	  {
	    int t = var.stored_surface_element(p)->GetHybridType();
	    int r = var.stored_surface_element(p)->GetOrder();
            int q = var.type_stored_surface_element(p);
	    elt_target(t, q)(r) = var.stored_surface_element(p);
	  }
	
	// shallow copy : only pointers are copied
        this->nb_surf_ = mesh.GetNbBoundary();
	surface_element.Reallocate(mesh.GetNbBoundary()*N);
        for (int n = 0; n < N; n++)
          for (int i = 0; i < mesh.GetNbBoundary(); i++)
            {
              int t = mesh.Boundary(i).GetHybridType();
              int r = 0, q = 0;
              if (n == 0)
                {
                  r = mesh_num.GetOrderQuadrature(i);
                  q = type;
                }
              else
                {
                  r = other_mesh_num(n-1)->GetOrderQuadrature(i);
                  q = other_type(n-1);
                }
              
              surface_element(i + n*this->nb_surf_) = elt_target(t, q)(r);
              if (elt_target(t, q)(r) == NULL)
                {
                  cout << "This case should not appear" << endl;
                  abort();
                }
            }
      }    
  }  
  

  //! Sets finite element to use for a set of elements
  template<class Dimension>
  void VarFiniteElement<Dimension>
  ::SetFiniteElement(const string& name_elt, const Vector<bool>& change_elt, int type, int n)
  {
    int rmax = 1;
    for (int k = 0; k < stored_reference_element.GetM(); k++)
      rmax = max(rmax, stored_reference_element(k)->GetOrder());
    
    ElementReference_Dim<Dimension>* null_ptr;
    null_ptr = NULL;
    TinyVector<Vector<ElementReference_Dim<Dimension>* >, 4> elt_target;
    for (int t = 0; t < 4; t++)
      {
	elt_target(t).Reallocate(rmax+1);
	elt_target(t).Fill(null_ptr);
	int id_elt = this->GetIdentityNumber(name_elt, t, type);
	for (int p = 0; p < stored_reference_element.GetM(); p++)
	  if (id_stored_reference_element(p) == id_elt)
	    elt_target(t)(stored_reference_element(p)->GetOrder()) = stored_reference_element(p);
      }
    
    for (int i = 0; i < reference_element.GetM(); i++)
      if (change_elt(i))
	{
	  int t = reference_element(i)->GetHybridType();
	  int r = reference_element(i)->GetOrder();
	  reference_element(i + n*this->nb_elt_) = elt_target(t)(r);
	  if (elt_target(t)(r) == NULL)
	    {
	      cout << "This case should not appear" << endl;
	      abort();
	    }
	}
  }
  
  
  //! returns the array of finite elements stored in the class
  template<class Dimension>
  void VarFiniteElement<Dimension>
  ::GetReferenceElement(Vector<const ElementReference_Dim<Dimension>* >& elt) const
  {
    elt.Reallocate(stored_reference_element.GetM());
    for (int i = 0; i < elt.GetM(); i++)
      elt(i) = stored_reference_element(i);
  }


  //! returns the array of surface finite elements stored in the class
  template<class Dimension>
  void VarFiniteElement<Dimension>
  ::GetSurfaceFiniteElement(Vector<const ElementReferenceTrace* >& elt) const
  {
    elt.Reallocate(stored_surface_element.GetM());
    for (int i = 0; i < elt.GetM(); i++)
      elt(i) = stored_surface_element(i);
  }


  template<>
  void VarFiniteElement<Dimension1>::SetAxisymGeometry(bool axisym)
  {
  }

  
  //! informs that the geometry is axisymmetric
  template<class Dimension>
  void VarFiniteElement<Dimension>::SetAxisymGeometry(bool axisym)
  {
    for (int i = 0; i < stored_reference_element.GetM(); i++)
      stored_reference_element(i)->GetGeometricElement().SetAxisymmetricGeometry();
  }
  
}

#define MONTJOIE_FILE_VAR_FINITE_ELEMENT_CXX
#endif

