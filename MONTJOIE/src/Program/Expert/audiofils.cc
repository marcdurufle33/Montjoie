#include "Montjoie.hxx"

using namespace Montjoie;

class ComputationFieldH : public InputDataProblem_Base
{
public:
  //! component number to be displayed on grids (-1 -> all components)
  //! 0 : Hx, 1 : Hy, 2: Hz, -2 : Hz without decibels
  int choice_field_output;

  //! precision for output
  int precision_output_file;

  //! height of measurement
  Real_wp h;
  
  //! predefined display grids
  Vector<GridInterpolationFull<Dimension2> > var_grid;

  //! list of polylines
  Vector<VectR2> liste_wire;

  // list of walls
  bool smart_loop;
  VectR2 ptA_walls, ptB_walls;
  VectReal_wp transmission_wall;
  VectR2 pos_borne; VectReal_wp height_borne;
  bool allow_reflection; bool reflection_floor;

  Vector<IVect> liste_wall_vertex;
  IVect wall_first_vertex, wall_second_vertex;

  //! current for each wire
  VectComplex_wp current_wire;

  Real_wp level_field_zeroDB, level_dB;
  Real_wp epsilon;

  Real_wp rho_metal, thickness_metal, mu_metal;
  Real_wp pos_metal;
  Matrix<Real_wp> coef_metal;
  bool ratio_pade; Complex_wp coef_pade;

  Real_wp xmin_fourier, xmax_fourier, ymin_fourier, ymax_fourier;
  int nb_points_fourier; Real_wp freq;
  
  // algorithm to find the minimum field in the room
  int type_algo_min; R2 target; R2 minBox, maxBox; int nb_points_x , nb_points_y;
  enum {AUTO, MINIMUM_BOX, TARGET_POINT};

  int type_algo_fourier;
  enum {INTERPOLATION, DIRECT, TEST_ERROR};
  TinyVector<int, 4> nb_pts_interpolation; Real_wp Longueur_interp;
  bool interp_spline_L;

  bool dirichlet;
  bool scalar_app, one_layer;
public:
  ComputationFieldH()
  {
    allow_reflection = false;
    reflection_floor = true;
    smart_loop = true;
    choice_field_output = 2;
    precision_output_file = OutputTypeEnum::SINGLE_PRECISION;
    type_algo_min = AUTO;

    h = 1.2;
    
    // 0 dB correspond a 400 mA
    level_field_zeroDB = 0.4;

    // level that should be reached (minimal value)
    level_dB = -3.0;
    
    // pour eviter un log(0)
    epsilon = 1e-15;

    type_algo_fourier = INTERPOLATION;
    nb_pts_interpolation.Init(11, 7, 11, 5); Longueur_interp = 30.0;
    interp_spline_L = true;
    
    rho_metal = -1.725e-8; thickness_metal = 0.01; mu_metal = 1.0;
    pos_metal = -0.1; ratio_pade = false;

    xmin_fourier = -60; xmax_fourier = 80;
    ymin_fourier = -60; ymax_fourier = 80;
    nb_points_fourier = 1000;

    dirichlet = false;
    scalar_app = false;
    one_layer = false;
  }

  // on lit le fichier de donnees
  void SetInputData(const string& keyword, const Vector<string>& parameters)
  {
    if (keyword == "SismoLine")
      {
	// output on a line
	int i = var_grid.GetM();
	var_grid.Resize(var_grid.GetM()+1);
	var_grid(i).SetInputData(keyword, parameters);
      }
    else if (keyword == "SismoPlane")
      {
	// output on a line
	int i = var_grid.GetM();
	var_grid.Resize(var_grid.GetM()+1);
	var_grid(i).SetInputData(keyword, parameters);
      }
    else if (keyword == "ComponentField")
      {
	choice_field_output = to_num<int>(parameters(0));
      }
    else if (keyword == "HeightMeasurement")
      {
	h = to_num<Real_wp>(parameters(0));
      }
    else if (keyword == "Metal")
      {
        ratio_pade = false;
        pos_metal = to_num<Real_wp>(parameters(0));
        if (parameters(1) == "FILE")
          {
            rho_metal = 1.0;
            coef_metal.ReadText(parameters(2));
          }
        else if (parameters(1) == "PADE")
          {
            ratio_pade = true;
            rho_metal = 1.0;
            coef_pade = to_num<Complex_wp>(parameters(2));
          }
        else
          {
            rho_metal = to_num<Real_wp>(parameters(1));
            mu_metal = to_num<Real_wp>(parameters(2));
            thickness_metal = to_num<Real_wp>(parameters(3));
          }
      }
    else if (keyword == "Frequency")
      {
        freq = to_num<Real_wp>(parameters(0));
      }
    else if (keyword == "NbPointsFourier")
      {
        nb_points_fourier = to_num<int>(parameters(0));        
      }
    else if (keyword == "FourierBox")
      {
        xmin_fourier = to_num<Real_wp>(parameters(0));
        xmax_fourier = to_num<Real_wp>(parameters(1));
        ymin_fourier = to_num<Real_wp>(parameters(2));
        ymax_fourier = to_num<Real_wp>(parameters(3));        
      }
    else if (keyword == "FourierAlgo")
      {
        if (parameters(0) == "Direct")
          type_algo_fourier = DIRECT;
        else if (parameters(0) == "TestError")
          type_algo_fourier = TEST_ERROR;
        else if (parameters(0) == "Interpolation")
          type_algo_fourier = INTERPOLATION;

        if (parameters.GetM() > 1)
          {
            if (parameters(1) == "Spline")
              interp_spline_L = true;
            else
              interp_spline_L = false;
          }
        
        if (parameters.GetM() > 2)
          Longueur_interp = to_num<Real_wp>(parameters(2));
      }
    else if (keyword == "NumberPointsInterpolation")
      {
        int N = min(nb_pts_interpolation.GetM(), int(parameters.GetM()));
        for (int i = 0; i < N; i++)
          nb_pts_interpolation(i) = to_num<int>(parameters(i));
      }
    else if (keyword == "AddWall")
      {
        R2 ptA, ptB;
        ptA(0) = to_num<Real_wp>(parameters(0));
        ptA(1) = to_num<Real_wp>(parameters(1));
        ptB(0) = to_num<Real_wp>(parameters(2));
        ptB(1) = to_num<Real_wp>(parameters(3));
        ptA_walls.PushBack(ptA);
        ptB_walls.PushBack(ptB);
        
        if (parameters(4) == "Faible")
          {
            // 3-5 dB
            transmission_wall.PushBack(0.63);
          }
        else if (parameters(4) == "Moyen")
          {
            // 6-15 dB
            transmission_wall.PushBack(0.3);
          }
        else if (parameters(4) == "Fort")
          {
            // 6-15 dB
            transmission_wall.PushBack(0.05);
          }
        else
          {
            Real_wp coef = to_num<Real_wp>(parameters(4));
            transmission_wall.PushBack(coef);
          }
      }
    else if (keyword == "AddBorne")
      {
        R2 pos;
        pos(0) = to_num<Real_wp>(parameters(0));
        pos(1) = to_num<Real_wp>(parameters(1));
        Real_wp h = 3.0;
        if (parameters.GetM() > 2)
          h = to_num<Real_wp>(parameters(2));
        
        pos_borne.PushBack(pos);
        height_borne.PushBack(h);
      }
    else if (keyword == "SmartLoop")
      {
        if (parameters(0) == "YES")
          smart_loop = true;
        else
          {
            smart_loop = false;
            if (parameters.GetM() > 1)
              {
                if (parameters(1) == "Reflection")
                  allow_reflection = true;
                else if (parameters(1) == "Floor")
                  {
                    allow_reflection = false;
                    reflection_floor = true;
                  }
                else
                  allow_reflection = false;
              }
          }
      }
    else if (keyword == "AddWire")
      {
	current_wire.PushBack(to_num<Complex_wp>(parameters(0)));
	if (parameters(1) == "PolyLine")
	  {
	    VectR2 poly_line(parameters.GetM()-2);
	    for (int k = 0; k < poly_line.GetM(); k++)
	      {
		Complex_wp z = to_num<Complex_wp>(parameters(2+k));
		poly_line(k).Init(real(z), imag(z));
	      }

	    liste_wire.PushBack(poly_line);
	  }
	else if (parameters(1) == "MultiLoop")
	  {
	    int n0 = to_num<int>(parameters(2));
	    int Nseg = to_num<int>(parameters(3));
	    int nf = to_num<int>(parameters(4));
	    Real_wp L0 = to_num<Real_wp>(parameters(5));
	    Real_wp Lseg = to_num<Real_wp>(parameters(6));
	    if (parameters(7) == "X")
	      {
		Real_wp x0 = to_num<Real_wp>(parameters(8));
		Real_wp xf = to_num<Real_wp>(parameters(9));
		Real_wp y0 = to_num<Real_wp>(parameters(10));
		Real_wp yf = to_num<Real_wp>(parameters(11));

		VectR2 poly_line;
		// premier segment
		poly_line.PushBack(R2(x0, y0)); poly_line.PushBack(R2(x0, yf));

		// segments pour aller jusqu'a xf (aller)
		Real_wp x = x0;
		int N = n0 + Nseg + nf;
		for (int i = 0; i < N; i++)
		  {
		    Real_wp L = Lseg;
		    if ((i < n0) || (i >= n0+Nseg))
		      L = L0;
		    
		    x += L;
		    if (i%2 == 0)
		      {
			poly_line.PushBack(R2(x, yf));
			poly_line.PushBack(R2(x, y0));
		      }
		    else
		      {
			poly_line.PushBack(R2(x, y0));
			poly_line.PushBack(R2(x, yf));
		      }		    
		  }
		
		if (abs(x-xf) > 1e-6)
		  {
		    cout << "Abscisse finale " << x+L0 << " differe de " << xf << endl;
		    abort();
		  }
		
		// segments de retour
		for (int i = N-1; i >= 0; i--)
		  {
		    Real_wp L = Lseg;
		    if ((i < n0) || (i >= n0+Nseg))
		      L = L0;
		    
		    x -= L;
		    if (i%2 == 0)
		      {
			poly_line.PushBack(R2(x, y0));
			if (i > 0)
			  poly_line.PushBack(R2(x, yf));
		      }
		    else
		      {
			poly_line.PushBack(R2(x, yf));
			poly_line.PushBack(R2(x, y0));
		      }		    
		  }
		
		DISP(poly_line);
                
		// on rajoute la boucle ainsi cree
		liste_wire.PushBack(poly_line);
	      }
	    else
	      {
		cout << "Pas implemente" << endl;
		abort();
	      }
	  }
	else
	  {
	    cout << "Type de boucle inconnu " << parameters(1) << endl;
	    abort();
	  }
      }
    else if (keyword == "AlgorithmMinimum")
      {
	if (parameters(0) == "AUTO")
	  type_algo_min = AUTO;
	else if (parameters(0) == "BOX")
	  {
	    type_algo_min = MINIMUM_BOX;
	    minBox.Init(to_num<Real_wp>(parameters(1)), to_num<Real_wp>(parameters(2)));
	    maxBox.Init(to_num<Real_wp>(parameters(3)), to_num<Real_wp>(parameters(4)));
	    nb_points_x = to_num<int>(parameters(5));
	    nb_points_y = to_num<int>(parameters(6));
	  }
	else if (parameters(0) == "POINT")
	  {
	    type_algo_min = TARGET_POINT;
	    target.Init(to_num<Real_wp>(parameters(1)), to_num<Real_wp>(parameters(2)));
	  }
	else
	  {
	    cout << "Algorithme de recherche de minimum inconnu " << parameters(0) << endl;
	    abort();
	  }
      }
    else if (keyword == "LevelDecibel")
      level_dB = to_num<Real_wp>(parameters(0));
    else if (keyword == "Model")
      {
        dirichlet = false;
        scalar_app = false;
        one_layer = false;
        
        if (parameters(0) == "Scalar")
          scalar_app = true;
        
        if (parameters.GetM() > 1)
          {
            if (parameters(1) == "Dirichlet")
              dirichlet = true;
            else if (parameters(1) == "OneLayer")
              one_layer = true;
          }        
      }
    else
      {
	cout << "Mot-cle inconnu : " << keyword << "endl";
	abort();
      }
  }

  int GetPositionPoint(const R2& pt1, const R2& pt2, const R2& pointA, const R2& pointB,
                       R2& intersec, int& pt_on_vertex)
  {
    // calcul de l'intersection
    Real_wp dx = pointB(0)-pointA(0);
    Real_wp dy = pointB(1)-pointA(1);
    Real_wp dxj = pt2(0)-pt1(0);
    Real_wp dyj = pt2(1)-pt1(1);
    Real_wp delta = dx*dyj-dy*dxj;
    Real_wp x, y, zero;
    SetComplexZero(zero);
    pt_on_vertex = 0;
    if (abs(delta) > 1e-10)
      {
	x = (dxj*dx*(pointA(1)-pt1(1))+pt1(0)*dyj*dx-pointA(0)*dy*dxj)/delta;
	y = -(dyj*dy*(pointA(0)-pt1(0))+pt1(1)*dxj*dy-pointA(1)*dx*dyj)/delta;
        intersec.Init(x, y);        
	if (((x-pt1(0))*(x-pt2(0))) <= 1e-14)
          if (((y-pt1(1))*(y-pt2(1))) <= 1e-14)
            if (((x-pointA(0))*(x-pointB(0))) <= 1e-14)
              if (((y-pointA(1))*(y-pointB(1))) <= 1e-14)
                {
                  if (intersec == pointA)
                    pt_on_vertex = 1;

                  if (intersec == pointB)
                    pt_on_vertex = 2;
                  
                  return 1;
                }
      }
    
    return 0;
  }
  
  void DetectVertexWall()
  {
    DISP(R2::threshold);
    int N = ptA_walls.GetM();
    VectR2 list_vertices(2*N);
    wall_first_vertex.Reallocate(N);
    wall_second_vertex.Reallocate(N);
    wall_first_vertex.Fill(-1);
    wall_second_vertex.Fill(-1);    
    for (int i = 0; i < N; i++)
      {
        list_vertices(2*i) = ptA_walls(i);
        list_vertices(2*i+1) = ptB_walls(i);
      }

    IVect permut(2*N); permut.Fill();
    Sort(list_vertices, permut);

    R2 ptA = list_vertices(0);
    int k = 0, n = 0;
    liste_wall_vertex.Clear();
    while (k < 2*N)
      {
        int k2 = k;
        while ((k2 < 2*N) && (list_vertices(k2) == ptA))
          k2++;

        IVect num_wall(k2-k);
        for (int q = k; q < k2; q++)
          {
            int p = permut(q);
            num_wall(q-k) = p/2;
            int m = p%2;
            if (m == 0)
              wall_first_vertex(num_wall(q-k)) = n;
            else
              wall_second_vertex(num_wall(q-k)) = n;
          }
        
        liste_wall_vertex.PushBack(num_wall);
        DISP(n); DISP(ptA); DISP(liste_wall_vertex(n));
        n++; k = k2;
        if (k2 < 2*N)
          ptA = list_vertices(k2);
      }

    DISP(wall_first_vertex);
    DISP(wall_second_vertex);
  }


  Real_wp GetCoefTransmission(const R2& ptA, const R2& ptB, const IVect& num_excluded)
  {
    R2 intersec;
    int nb_vertices = liste_wall_vertex.GetM();
    Vector<bool> WallToBeTreated(nb_vertices);

    Real_wp coef(1);
    WallToBeTreated.Fill(true);
    for (int i = 0; i < num_excluded.GetM(); i++)
      WallToBeTreated(num_excluded(i)) = false;
    
    for (int j = 0; j < ptA_walls.GetM(); j++)
      if (WallToBeTreated(j))
        {
          int pt_on_vertex = 0;
          int s = GetPositionPoint(ptA, ptB, ptA_walls(j), ptB_walls(j), intersec, pt_on_vertex);
          int n1 = wall_first_vertex(j);
          int n2 = wall_second_vertex(j);
          if (pt_on_vertex == 1)
            {
              for (int q = 0; q < liste_wall_vertex(n1).GetM(); q++)
                WallToBeTreated(liste_wall_vertex(n1)(q)) = false;
            }                  
          else if (pt_on_vertex == 2)
            {
              for (int q = 0; q < liste_wall_vertex(n2).GetM(); q++)
                WallToBeTreated(liste_wall_vertex(n2)(q)) = false;
            }
          
          if (s == 1)
            coef *= transmission_wall(j);
        }

    return coef;
  }

  void ComputeNormaleWall(VectR2& normale)
  {
    int nb_walls = ptA_walls.GetM();
    normale.Reallocate(nb_walls);
    for (int j = 0; j < nb_walls; j++)
      {
        R2 u = ptA_walls(j) - ptB_walls(j);
        R2 v(u(1), -u(0));
        Real_wp coef = 1.0/Norm2(v); Mlt(coef, v);
        normale(j) = v;
      }
  }

  void ComputeFieldAura(const VectR2& Points2D, VectReal_wp& val)
  {
    Real_wp coef_reflexion = 0.25;
    VectR2 normale;
    if (allow_reflection)
      {
        ComputeNormaleWall(normale);
        DISP(normale);
      }

    Real_wp freq = 2.4e9;
    Real_wp omega = 2.0*pi_wp*freq / 3e8; DISP(omega);
    IVect num_excluded;
    for (int k = 0; k < Points2D.GetM(); k++)
      {
        R2 pt = Points2D(k); R2 res;
        Real_wp Intens(0); 
        for (int i = 0; i < pos_borne.GetM(); i++)
          {
            if (!allow_reflection)
              {
                IVect num_excluded;
                Real_wp coef = GetCoefTransmission(pt, pos_borne(i), num_excluded);
                Real_wp r2 = pt.DistanceSquare(pos_borne(i));
                if (r2 < 0.1)
                  r2 = 0.1;

                Real_wp r = sqrt(r2);                
                Complex_wp fieldU(0, 0);
                Complex_wp phase = exp(Iwp*omega*r) / r;
                fieldU += phase*coef;
                
                // reflexion plancher
                if (reflection_floor)
                  {
                    Real_wp z = height_borne(i) + h;
                    r = sqrt(r*r + z*z);
                    phase = exp(Iwp*omega*r) / r;
                    fieldU -= phase*coef;
                  }
                
                Intens += absSquare(fieldU);
              }
            else
              {
                Complex_wp fieldU(0, 0);
                
                // partie directe 
                IVect num_excluded;
                Real_wp coef = GetCoefTransmission(pt, pos_borne(i), num_excluded);
                Real_wp r = pt.Distance(pos_borne(i));
                if (r < 0.3)
                  r = 0.3;
                
                Complex_wp phase = exp(Iwp*omega*r) / r;
                fieldU += phase*coef;

                // reflexion plancher
                if (reflection_floor)
                  {
                    Real_wp z = height_borne(i) + h;
                    r = sqrt(r*r + z*z);
                    phase = exp(Iwp*omega*r) / r;
                    fieldU += phase*coef;
                  }
                
                // partie reflexion
                for (int j = 0; j < ptA_walls.GetM(); j++)
                  {
                    R2 v = pt - ptA_walls(j);
                    Real_wp alpha = 2.0*DotProd(normale(j), v);
                    R2 reflexB = pt - alpha*normale(j);
                    //if (k == 0)
                    //{
                        //DISP(j); DISP(pos_borne(i)); DISP(pt); DISP(ptA_walls(j)); DISP(ptB_walls(j)); DISP(reflexB);
                        int pt_on_vertex; R2 intersecB;
                        int s = GetPositionPoint(pos_borne(i), reflexB, ptA_walls(j), ptB_walls(j), intersecB, pt_on_vertex);
                        //DISP(intersecB); DISP(s);
                        if (s == 1)
                          {
                            //DISP(pt_on_vertex);
                            if (pt_on_vertex == 1)
                              {
                                int n1 = wall_first_vertex(j);
                                num_excluded = liste_wall_vertex(n1);
                              }
                            else if (pt_on_vertex == 2)
                              {
                                int n2 = wall_second_vertex(j);
                                num_excluded = liste_wall_vertex(n2);
                              }
                            else
                              {
                                num_excluded.Clear();
                                num_excluded.Reallocate(1);
                                num_excluded(0) = j;
                              }
                            
                            Real_wp coef1 = GetCoefTransmission(intersecB, pos_borne(i), num_excluded);
                            Real_wp coef2 = GetCoefTransmission(intersecB, pt, num_excluded);
                            Real_wp r = intersecB.Distance(pos_borne(i)) + intersecB.Distance(pt);
                            if (r < 0.3)
                              r = 0.3;

                            //DISP(coef1); DISP(coef2); DISP(r);
                            Complex_wp phase = exp(Iwp*omega*r) / r;
                            fieldU += phase*coef1*coef2*coef_reflexion;
                          }
                        //}
                  }
                
                // on rajoute le resultat
                Intens += absSquare(fieldU);
              }
          }
        
        val(k) = sqrt(Intens);
      }
  }
  
  // calcul du champ H en espace libre a une hauteur h0
  void ComputeFieldH(const VectR2& Points2D, const Vector<VectR2>& pt_wire,
                     const VectComplex_wp& I0_wire, const Real_wp& h0,
                     Vector<TinyVector<Complex_wp, 3> >& fieldH)
  {
    fieldH.Reallocate(Points2D.GetM());
    fieldH.Fill(TinyVector<Complex_wp, 3>(Real_wp(0), Real_wp(0), Real_wp(0)));
    Real_wp coef = 1.0/(4.0*pi_wp);
    R2 ur, vr, vec_u;
    for (int n = 0; n < pt_wire.GetM(); n++)
      for (int i = 0; i < pt_wire(n).GetM()-1; i++)
	{
	  vec_u(0) = pt_wire(n)(i+1)(0) - pt_wire(n)(i)(0);
	  vec_u(1) = pt_wire(n)(i+1)(1) - pt_wire(n)(i)(1);
	  
	  Real_wp L = Norm2(vec_u);
	  if (L == Real_wp(0))
	    {
	      cout << "Fil de longueur nulle " << endl;
	      cout << "polyline = " << pt_wire(n) << endl;
	      abort();
	    }
	  
	  Real_wp invNorm = 1.0/L;
	  ur.Init(invNorm*vec_u(0), invNorm*vec_u(1));
	  vr.Init(-ur(1), ur(0));
	  //DISP(i); DISP(ur); DISP(vr); DISP(L); DISP(current_wire(n));
	  
	  for (int k = 0; k < Points2D.GetM(); k++)
	    {
	      R2 diff = Points2D(k)-pt_wire(n)(i);
	      Real_wp beta = DotProd(diff, ur);
	      Real_wp gamma = DotProd(diff, vr);
	      Real_wp offset = h0*h0 + gamma*gamma;
	      Real_wp valInt = (L - beta) / sqrt(offset + square(L-beta)) + beta / sqrt(offset + beta*beta);
	      Complex_wp vloc = I0_wire(n) * valInt / offset;
	      //DISP(k); DISP(diff); DISP(beta); DISP(gamma); DISP(offset); DISP(valInt);
	      fieldH(k)(0) -= h0*vloc*vr(0);
	      fieldH(k)(1) -= h0*vloc*vr(1);
	      fieldH(k)(2) += gamma*vloc;
	    }
	}
    
    for (int k = 0; k < Points2D.GetM(); k++)
      fieldH(k) *= coef;
  }

  // on met u dans v en rajoutant la condition de periodicite
  void FillPeriodic(int nb_points_x, int nb_points_y, const VectComplex_wp& u, VectComplex_wp& v)
  {
    for (int i = 0; i < nb_points_x; i++)
      for (int j = 0; j < nb_points_y; j++)
        v(i*(nb_points_y+1) + j) = u(i*nb_points_y + j);

    for (int i = 0; i < nb_points_x; i++)
      v(i*(nb_points_y+1) + nb_points_y) = u(i*nb_points_y);

    for (int j = 0; j < nb_points_y; j++)
      v(nb_points_x*(nb_points_y+1) + j) = u(j);

    v((nb_points_x+1)*(nb_points_y+1) - 1) = u(0);
  }

  // calcul de H quand J est une gaussienne
  void ComputeGaussianFieldH(const VectR2& Points2D, const Real_wp& h, const Real_wp& omega,
                             Vector<TinyVector<Complex_wp, 3> >& fieldH)
  {
    int r = 12;
    Globatto<Real_wp> lob;
    lob.ConstructQuadrature(r);
    
    fieldH.Reallocate(Points2D.GetM());
    Real_wp alpha = log(1e6), beta = sqrt(alpha/pi_wp);
    R3 ptRp, vec_v, ptR; TinyVector<Complex_wp, 3> Jc, vloc;
    Matrix3_3sym_Complex_wp hessian_phi, dyadic_G;
    R3_Complex_wp grad_phi; Complex_wp phi;
    Complex_wp czero(0, 0); Real_wp omega2 = omega*omega; R3 center(0, 0, 1);
    for (int j1 = 0; j1 < lob.GetNbPointsQuad(); j1++)
      for (int j2 = 0; j2 < lob.GetNbPointsQuad(); j2++)
        for (int j3 = 0; j3 < lob.GetNbPointsQuad(); j3++)
          {
            ptRp.Init(2.0*lob.Points(j1)-1.0, 2.0*lob.Points(j2)-1.0, 2.0*lob.Points(j3)-1.0);
            Real_wp coef = 8.0*lob.Weights(j1)*lob.Weights(j2)*lob.Weights(j3)
              *beta*exp(-alpha*DotProd(ptRp, ptRp));
            
            ptRp += center;
            Jc.Init(Complex_wp(coef, 0), czero, czero);
	  
            for (int k = 0; k < Points2D.GetM(); k++)
              {
                R3 pt3D(Points2D(k)(0), Points2D(k)(1), h);
                TransparencySolver_Base::ComputeGreenKernel(pt3D, ptRp, omega, phi,
                                                            grad_phi, hessian_phi);              
                
                // G = \phi(x,y) I + 1/k^2 \nabla_y \nabla_y \phi(x,y)
                //Copy(hessian_phi, dyadic_G); Mlt(Complex_wp(1.0/omega2), dyadic_G);
                
                //dyadic_G(0, 0) += phi;
                //dyadic_G(1, 1) += phi;
                //dyadic_G(2, 2) += phi;
                
                //Mlt(dyadic_G, Jc, vloc);
                //fieldE(k) += vloc;
              
                TimesProd(grad_phi, Jc, vloc);
                fieldH(k) += vloc;
	    }
	}
  }

  // calcul du noyau de Green sur une liste de points
  void ComputeGreenFieldH(const VectR2& Points2D, const Real_wp& h, const Complex_wp& omega,
                          Vector<TinyVector<Complex_wp, 3> >& fieldH)
  {
    R3 X(0, 0, 1);
    R3 polar(1, 0, 0);
    fieldH.Reallocate(Points2D.GetM());
    for (int k = 0; k < Points2D.GetM(); k++)
      {
        R3 Y(Points2D(k)(0), Points2D(k)(1), h);

        // in 3-D \phi(x,y) = exp(ik|x-y|) / (4 pi |x-y| )
        R3 xMinusy = X - Y;

        Real_wp T = Norm2(xMinusy); Complex_wp arg = Iwp*omega*T;
        Real_wp inv_T = 1.0/T,inv_T2 = inv_T*inv_T;

        Complex_wp phi = exp(arg) * inv_T / (4*pi_wp); Complex_wp phi_div_T = phi * inv_T;

        // evaluation of its gradient
        R3_Complex_wp grad_phi; grad_phi = xMinusy;    
        Complex_wp alpha = (-Iwp*omega + inv_T)*phi_div_T; Mlt( alpha , grad_phi);    

        TimesProd(grad_phi, polar, fieldH(k));
      }
  }

  // calcul de log(y/Hz) sur les points d'interpolation
  Real_wp FindPolynomialCoef(const VectComplex_wp& Hz, const VectReal_wp& PointsX, const VectReal_wp& PointsY, const Real_wp& xc, Matrix<Complex_wp>& Hi, TinyVector<bool, 2>& sign_Hi)
  {
    int i0 = -1, imax = PointsX.GetM();
    int j0 = -1, jmax = PointsY.GetM();
    Real_wp ymax = Longueur_interp, xmax = xc + Longueur_interp;
    for (int i = 0; i < PointsX.GetM(); i++)
      {
        if (abs(PointsX(i)) <= 1e-12)
          i0 = i;
        
        if (PointsX(i) > xmax)
          {
            imax = i;        
            break;
          }
      }

    for (int j = 0; j < PointsY.GetM(); j++)
      {
        if (abs(PointsY(j)) <= 1e-12)
          j0 = j;
        
        if (PointsY(j) > ymax)
          {
            jmax = j;        
            break;
          }
      }
    
    //DISP(i0); DISP(j0); DISP(imax); DISP(jmax);
    //DISP(PointsX(imax)); DISP(PointsY(jmax));

    // pour evaluer la derivee de y/Hz
    int order = 4;
    Globatto<Real_wp> lob;
    VectReal_wp pts(order+1);

    for (int i = 0; i <= order; i++)
      pts(i) = PointsY(j0+i);
    
    lob.AffectPoints(pts);
    VectReal_wp dphi(order+1);
    for (int i = 0; i <= order; i++)
      dphi(i) = lob.EvaluatePhiGrad(i, pts(0));

    int Ny = PointsY.GetM();

    // boucle sur tous les points x pour evaluer
    // log(Re(y/Hz)) sur des points py reguliers
    int N = nb_pts_interpolation(0);
    Real_wp py_max = log(1.0 + ymax);
    VectReal_wp Yi(N), pYi(N);
    for (int i = 0; i < N; i++)
      Yi(i) = exp(py_max * Real_wp(i) / (N-1)) - 1.0;
    
    Matrix<Complex_wp> Hy(imax-i0+1, N);
    bool pos_real = true, pos_imag = true;
    for (int i = i0; i <= imax; i++)
      {
        SplineInterpolation<Complex_wp> spline;
        VectReal_wp xdiv(jmax-j0+1);
        VectComplex_wp ydiv(jmax-j0+1);
        
        // valeur pour y = 0 en calculant la derivee
        Complex_wp val = 0;
        for (int j = 0; j <= order; j++)
          val += dphi(j)*Hz(i*Ny+j+j0);
        
        xdiv(0) = 0.0; ydiv(0) = log(abs(real(1.0/val))) + Iwp*log(abs(imag(1.0/val)));
        if (real(1.0/val) < 0)
          pos_real = false;

        if (imag(1.0/val) < 0)
          pos_imag = false;

        // autres valeurs
        for (int j = j0+1; j <= jmax; j++)
          {
            //xdiv(j-j0) = log(1.0 + PointsY(j));
            xdiv(j-j0) = PointsY(j);
            ydiv(j-j0) = log(abs(real(PointsY(j) / Hz(i*Ny + j)))) + Iwp*log(abs(imag(PointsY(j)/Hz(i*Ny+j))));
          }

        // interpolation par spline cubique
        spline.Init(xdiv, ydiv);
        for (int j = 0; j < N; j++)
          Hy(i-i0, j) = spline.Evaluate(Yi(j));
      }
    
    int N1 = nb_pts_interpolation(1), N2 = nb_pts_interpolation(2);
    VectReal_wp Xi(N1+N2-1);
    Real_wp p1_max = log(xc+1.0), p2_max = log(1.0 + (xmax-xc));
    for (int i = 0; i < N1; i++)
      Xi(i) = xc - exp(p1_max * Real_wp(N1-i-1) / (N1-1)) + 1.0;
    
    for (int i = 0; i < N2; i++)
      Xi(N1-1+i) = xc + exp(p2_max * Real_wp(i) / (N2-1)) - 1.0;
    
    //DISP(Xi); DISP(Yi); 
    Hi.Reallocate(N1+N2-1, N);
    for (int j = 0; j < N; j++)
      {
        SplineInterpolation<Complex_wp> spline;
        VectReal_wp xdiv(imax-i0+1);
        VectComplex_wp ydiv(imax-i0+1);

        // autres valeurs
        for (int i = i0; i <= imax; i++)
          {
            //xdiv(i-i0) = log(1.0 + PointsX(i));
            xdiv(i-i0) = PointsX(i);
            ydiv(i-i0) = Hy(i-i0, j);
          }
        
        // interpolation par spline cubique
        spline.Init(xdiv, ydiv);
        for (int i = 0; i < Xi.GetM(); i++)
          {
            Hi(i, j) = spline.Evaluate(Xi(i));
          }        
      }
    
    Hi.Write("Hi.dat");
    sign_Hi(0) = pos_real; sign_Hi(1) = pos_imag;
    //DISP(Hi);
    
    // on calcule l'erreur maximale avec cette interpolation
    order = 4; 
    pts.Reallocate(order+1); pts.Fill(); Mlt(Real_wp(1)/order, pts);
    lob.AffectPoints(pts); Real_wp dy = py_max / (N-1);
    Real_wp dx1 = p1_max / (N1-1), dx2 = p2_max / (N2-1);
    Real_wp err_max = 0;

    for (int i = i0; i <= imax; i++)
      for (int j = j0+1; j <= jmax; j++)
        {
          Real_wp x = PointsX(i), y = PointsY(j);
          Real_wp py = log(1.0 + y);
          if (py > py_max)
            continue;
          
          Real_wp px=0;
          if (x < xc)
            px = -log(1.0 + abs(x-xc));
          else
            px = log(1.0 + abs(x-xc));
          
          if (px > p2_max)
            continue;
          
          // calcul des fcts de base en y
          int iy = toInteger(round(py / py_max*N));
          int iymin = iy -order/2, iymax = iy+order/2;
          if (iymin < 0)
            {
              iymax += -iymin;
              iymin = 0;
            }

          if (iymax >= N)
            {
              iymin -= (iymax-N+1);
              iymax = N-1;
            }
          
          Real_wp yloc = Real_wp(py - iymin*dy) / (order*dy);
          VectReal_wp phi_y;
          lob.ComputeValuesPhiRef(yloc, phi_y);

          // calcul des fcts de base en x
          VectReal_wp phi_x;
          int ixmin = 0, ixmax = order, ix;
          if (x < xc)
            {
              ix = N1 - toInteger(round(-px/p1_max*N1));
              ixmin = ix - order/2; ixmax = ix + order/2;
              if (ixmin < 0)
                {
                  ixmax += -ixmin;
                  ixmin = 0;
                }
              
              if (ixmax >= N1)
                {
                  ixmin -= (ixmax-N1+1);
                  ixmax = N1-1;
                }
              
              Real_wp xloc = (px + p1_max - ixmin*dx1) / (order*dx1);
              lob.ComputeValuesPhiRef(xloc, phi_x);
            }
          else
            {
              ix = toInteger(round(px/p2_max*N2));
              ixmin = ix - order/2; ixmax = ix + order/2;
              if (ixmin < 0)
                {
                  ixmax += -ixmin;
                  ixmin = 0;
                }
              
              if (ixmax >= N2)
                {
                  ixmin -= (ixmax-N2+1);
                  ixmax = N2-1;
                }
              
              Real_wp xloc = (px - ixmin*dx2) / (order*dx2);
              lob.ComputeValuesPhiRef(xloc, phi_x);
              ixmin += N1-1; ixmax += N1-1;
            }
          
          // on effectue l'interpolation
          Complex_wp val(0, 0);
          for (ix = ixmin; ix <= ixmax; ix++)
            for (iy = iymin; iy <= iymax; iy++)
              val += phi_x(ix-ixmin)*phi_y(iy-iymin)*Hi(ix, iy);
          
          Real_wp vr = exp(real(val));
          if (!pos_real)
            vr = -vr;

          Real_wp vi = exp(imag(val));
          if (!pos_imag)
            vi = -vi;
          
          //DISP(x); DISP(y); DISP(px); DISP(py); DISP(phi_x); DISP(phi_y);
          //DISP(ixmin); DISP(ixmax); DISP(iymin); DISP(iymax);
          Complex_wp vloc = y/Complex_wp(vr, vi);
          Real_wp err = abs(vloc - Hz(i*Ny+j)) / abs(vloc); //DISP(err); DISP(vloc); DISP(Hz(i*Ny+j));
          if (abs(vloc) > 1e-6)
            err_max = max(err, err_max);
        }
    
    //DISP(pos_real); DISP(pos_imag); DISP(err_max);
    
    return err_max;
    
    /*
    // interpolation polynomiale
    order = 6; int nb_pol = (order+1)*(order+2)/2;
    Matrix<Real_wp> VDM(N*N, nb_pol);
    VectComplex_wp coefY(N*N);
    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++)
        {
          int npoint = i*N + j, num = 0;
          Real_wp x = px_max*i/(N-1), y = py_max*j/(N-1);
          for (int k = 0; k <= order; k++)
            for (int m = 0; m <= order-k; m++)
              {
                VDM(npoint, num) = pow(x, k)*pow(y, m);
                num++;
              }
          
          coefY(npoint) = Hi(i, j);
        }
    
    VectReal_wp tau;
    GetQR(VDM, tau);
    
    SolveQR(VDM, tau, coefY);
    DISP(coefY);

    Matrix<Complex_wp> EvalG(N, N);
    Real_wp err = 0;
    for (int i = 0; i < N; i++)
      for (int j = 0; j < N; j++)
        {
          int npoint = i*N + j, num = 0;
          Real_wp x = px_max*i/(N-1), y = py_max*j/(N-1);
          Complex_wp val(0, 0);
          for (int k = 0; k <= order; k++)
            for (int m = 0; m <= order-k; m++)
              {
                val += coefY(num)*pow(x, k)*pow(y, m);
                num++;
              }

          EvalG(i, j) = val;
          err = max(err, abs(val - Hi(i, j)) / abs(val));
        }
    
    DISP(err);
    EvalG.Write("EvalG.dat"); */
  }

  void ComputeFieldFourier(const Vector<VectR2>& pt_wire, const VectComplex_wp& I0_wire,
                           const Real_wp& xmin, const Real_wp& xmax, const Real_wp& ymin,
                           const Real_wp& ymax, int nb_points_x, int nb_points_y,
                           VectReal_wp& PointsX, VectReal_wp& PointsY,
                           VectComplex_wp& Hx_grid, VectComplex_wp& Hy_grid, VectComplex_wp& Hz_grid)
  {
    Real_wp delta_fourier = 10;
    
    GridInterpolationFull<Dimension2> grid;
    grid.SetPlaneOutput(xmin, xmax, ymin, ymax, nb_points_x+1, nb_points_y+1);
    
    Linspace(xmin, xmax, nb_points_x+1, PointsX);
    Linspace(ymin, ymax, nb_points_y+1, PointsY);
    
    VectR2 Points2D(nb_points_x * nb_points_y);
    for (int i = 0; i < nb_points_x; i++)
      for (int j = 0; j < nb_points_y; j++)
        Points2D(i*nb_points_y + j).Init(PointsX(i), PointsY(j));

    // troncature pour forcer la periodicite
    VectReal_wp tronc(nb_points_x*nb_points_y);
    for (int i = 0; i < Points2D.GetM(); i++)
      {
        Real_wp troncX = 1.0;
        if (Points2D(i)(0) < xmin + delta_fourier)
          {
            Real_wp x = (Points2D(i)(0)-xmin)/delta_fourier;
            troncX = x*x*(-2.0*x+3.0);
          }

        if (Points2D(i)(0) > xmax - delta_fourier)
          {
            Real_wp x = (xmax-Points2D(i)(0))/delta_fourier;
            troncX = x*x*(-2.0*x+3.0);
          }

        Real_wp troncY = 1.0;
        if (Points2D(i)(1) < ymin + delta_fourier)
          {
            Real_wp x = (Points2D(i)(1)-ymin)/delta_fourier;
            troncY = x*x*(-2.0*x+3.0);
          }

        if (Points2D(i)(1) > ymax - delta_fourier)
          {
            Real_wp x = (ymax-Points2D(i)(1))/delta_fourier;
            troncY = x*x*(-2.0*x+3.0);
          }

        tronc(i) = troncX*troncY;
      }
    
    Real_wp h0(pos_metal);
    Vector<TinyVector<Complex_wp, 3> > fieldH;

    // calcul du champ incident a la position du metal
    ComputeFieldH(Points2D, pt_wire, I0_wire, h0, fieldH);
    
    // initialisation du calcul de fft
    FftInterface<Complex_wp> fft2d;    
    fft2d.Init(nb_points_x, nb_points_y);

    //typedef __float128 Real_q;
    typedef Real_wp Real_q;
    typedef complex<Real_q> Complex_q;
    Real_wp pi_q = acos(Real_q(-1));
    Complex_wp Iq(0, 1);
    
    Real_q c0(299792458);
    //freq = c0/2;
    Complex_q eps0(1, 0.0), mu0(1);
    Complex_q kwave = Real_wp(2)*pi_q*freq*sqrt(eps0*mu0) / c0; //DISP(kwave);
    

    // pulsations spatiales ou on evalue la fft de H
    Vector<Real_q> nu_x, nu_y;
    nu_x.Reallocate(nb_points_x);
    nu_y.Reallocate(nb_points_y);
    for (int i = 0; i < nb_points_x/2; i++)
      {
        nu_x(i) = 2.0*pi_q*Real_q(i)/(xmax-xmin);
        nu_x(nb_points_x-1-i) = -2.0*pi_q*Real_q(i+1)/(xmax-xmin);
      }
    
    for (int i = 0; i < nb_points_y/2; i++)
      {
        nu_y(i) = 2.0*pi_q*Real_q(i)/(ymax-ymin);
        nu_y(nb_points_y-1-i) = -2.0*pi_q*Real_q(i+1)/(ymax-ymin);
      }

    //DISP(nu_x); DISP(nu_y);
    Vector<Complex_q> nu_z(nb_points_x * nb_points_y);
    Vector<Complex_q> nup_z(nb_points_x * nb_points_y);
    Vector<Complex_q> nu3_z(nb_points_x * nb_points_y);
    //Complex_wp eps_m(4.0, 0.2), mu_m(2.0);
    //Complex_wp eps_ext(2.0, 0.3), mu_ext(1.5);

    Real_q omega = 2.0*pi_q*freq/c0, sigma = PhysicalConstant::impedance0/rho_metal;
    Complex_q eps_m = Real_q(1) + Iwp*sigma / omega, mu_m(mu_metal); //DISP(eps_m);
    Complex_q eps_ext(1), mu_ext(1);
    
    Complex_q kwave_m = Real_wp(2)*pi_q*freq*sqrt(eps_m*mu_m) / c0; //DISP(kwave_m);
    Complex_q kwave_ext = Real_wp(2)*pi_q*freq*sqrt(eps_ext*mu_ext) / c0;
    for (int i = 0; i < nb_points_x; i++)
      for (int j = 0; j < nb_points_y; j++)
        {
          Real_q kx2_ky2 = nu_x(i)*nu_x(i) + nu_y(j)*nu_y(j);
          nu_z(i*nb_points_y + j) = sqrt(kwave*kwave - kx2_ky2);
          if (imag(nu_z(i*nb_points_y + j)) < 0)
            nu_z(i*nb_points_y + j) = -nu_z(i*nb_points_y + j);
          
          nup_z(i*nb_points_y + j) = sqrt(kwave_m*kwave_m - kx2_ky2);
          if (imag(nup_z(i*nb_points_y + j)) > 0)
            nup_z(i*nb_points_y + j) = -nup_z(i*nb_points_y + j);

          nu3_z(i*nb_points_y + j) = sqrt(kwave_ext*kwave_ext - kx2_ky2);
          if (imag(nu3_z(i*nb_points_y + j)) > 0)
            nu3_z(i*nb_points_y + j) = -nu3_z(i*nb_points_y + j);
        }

    VectComplex_wp  Hx, Hy, Hz;    
    Hz.Reallocate(nb_points_x*nb_points_y);
    
    Hx = Hz; Hy = Hz;
    Hx_grid.Reallocate((nb_points_x+1)*(nb_points_y+1));
    Hy_grid.Reallocate((nb_points_x+1)*(nb_points_y+1));
    Hz_grid.Reallocate((nb_points_x+1)*(nb_points_y+1));
    for (int i = 0; i < Hz.GetM(); i++)
      {
        Hz(i) = fieldH(i)(2); Hx(i) = fieldH(i)(0); Hy(i) = fieldH(i)(1);
        Hz(i) *= tronc(i);
      }

    FillPeriodic(nb_points_x, nb_points_y, Hz, Hz_grid);
    //WriteMatlab(Hz_grid, grid, "Hz0.dat", precision_output_file);
    
    // on calcule les composantes dans le domaine spectral
    fft2d.ApplyForward(Hx);
    fft2d.ApplyForward(Hy);
    fft2d.ApplyForward(Hz);

    Hx.Zero(); Hy.Zero();

    // calcul du champ diffracte sur un autre h
    Real_wp h1(h);
    Real_wp z0(h0), z1(h0 - thickness_metal);
    Real_wp dz = h1 - h0;
    Real_wp Err1(0), Err2(0), Err3(0), Err4(0), Err5(0), Err6(0), Err7(0), Err8(0);
    int nb_points = 0;
    VectComplex_wp ratio(nb_points_x); ratio.Zero();
    VectReal_wp kx_over_k(nb_points_x); kx_over_k.Zero();
    bool ratio_file = false;
    if (coef_metal.GetM() > 0)
      ratio_file = true;

    for (int i = 0; i < nb_points_x; i++)
      for (int j = 0; j < nb_points_y; j++)
        {
          Complex_q HxI, HyI, HzI;
          HxI = Hx(i*nb_points_y + j);
          HyI = Hy(i*nb_points_y + j);
          HzI = Hz(i*nb_points_y + j);

          Real_q kx = nu_x(i), ky = nu_y(j);
          Complex_q kz = nu_z(i*nb_points_y+j);
          if (abs(imag(kz)*dz) > 20.0)
            {
              Hx(i*nb_points_y + j) = 0;
              Hy(i*nb_points_y + j) = 0;
              Hz(i*nb_points_y + j) = 0;
              continue;
            }
          else
            nb_points++;
          
          //if (j == 0) { DISP(i); DISP(j); DISP(kx); DISP(kx/kwave); DISP(kz); }
          
          // on force kx HxI + ky HyI - kz HzI = 0
          if ( (abs(kx)  == Real_q(0) ) && (abs(ky) == Real_q(0)))
            {
              //DISP(HzI);
              HzI = Complex_q(0, 0);
            }
          else
            {
              if (abs(kx) >= abs(ky))
                HxI = (kz*HzI - ky*HyI) / kx;
              else
                HyI = (kz*HzI - kx*HxI) / ky;
            }
          
          Complex_q ExI = (ky*HzI+kz*HyI)/eps0;
          Complex_q EyI = (-kz*HxI-kx*HzI)/eps0;

          //DISP(HxI); DISP(HyI); DISP(HzI);
          //DISP(ExI); DISP(EyI);
          
          Complex_q phase = exp(Iq*kz*dz);
          Complex_q HxD, HyD, HzD;
          
          Complex_q kz1 = nup_z(i*nb_points_y + j);
          Complex_q kz2 = -kz1;
          Complex_q kz3 = nu3_z(i*nb_points_y + j);
          Complex_q phase1 = exp(Iq*kz1*(z1-z0));
          Complex_q phase2 = exp(Iq*kz2*(z1-z0));

          //DISP(phase);
          
          // cas scalaire (Helmholtz)
          if (ratio_file)
            {
              Real_wp kxy = sqrt(kx*kx + ky*ky);
              if (kxy > 5.0)
                HzD = 0.0;
              else
                {
                  int n = toInteger(round(kxy/0.025));
                  Complex_wp alpha = coef_metal(n, 1) + Iwp*coef_metal(n, 2);
                  HzD = alpha*HzI;
                }
            }
          else if (ratio_pade)
            {
              Real_wp kxy = sqrt(kx*kx + ky*ky);
              Complex_wp alpha = coef_pade / (kxy - coef_pade);
              HzD = alpha*HzI;
            }
          else if (scalar_app)
            {
              Complex_q u_inc = HzI;
              Complex_q alpha = mu_m * nu_z(i*nb_points_y+j) / (nup_z(i*nb_points_y+j) * eps0);
              
              if (dirichlet)
                HzD = -u_inc;
              else if (one_layer)
                {
                  Complex_q A = -(Real_wp(1) + alpha)/(Real_wp(1)-alpha)*u_inc;
                  HzD = -A; HxD = 0; HyD = 0;
                  if ((j == 0) && (i == 0))
                    {
                      DISP(alpha);
                      DISP(-(Real_wp(1) + alpha)/(Real_wp(1)-alpha));
                      DISP(HzD);
                    }
                }
              else
                {
                  Complex_q alpha_p = mu_ext * nup_z(i*nb_points_y+j) / (nu3_z(i*nb_points_y+j) * mu_m);
                  Complex_q beta1 = (Real_wp(1) + alpha_p) / (Real_wp(2)*alpha_p * phase1);
                  Complex_q beta2 = (-Real_wp(1) + alpha_p) / (Real_wp(2)*alpha_p * phase2);
                  Complex_q c1 = beta1 + beta2, c2 = beta1 - beta2;
                  Complex_q A = - (alpha*c1 + c2) / (c2 - alpha*c1) * u_inc;
                  HzD = -A; HxD = 0; HyD = 0;                
                }
            }
          else
            {
              // cas vectoriel (Maxwell)
              if (dirichlet)
                {
                  TinyMatrix<Complex_q, General, 2, 2> A;
                  TinyVector<Complex_q, 2> x, b;
                  Complex_q cx = -kx/kz, cy = -ky/kz;
                  
                  A(0, 0) = ky*cx/eps0; A(0, 1) = (-kz+ky*cy)/eps0;
                  A(1, 0) = (kz-kx*cx)/eps0; A(1, 1) = -kx*cy/eps0;

                  b(0) = -ExI; b(1) = -EyI;

                  GetInverse(A);
                  Mlt(A, b, x);

                  HxD = x(0);
                  HyD = x(1);
                  HzD = (-kx*HxD-ky*HyD)/kz;
                  
                  Complex_q ExD = (ky*HzD - kz*HyD) / eps0;
                  Complex_q EyD = (kz*HxD - kx*HzD) / eps0;

                  Real_wp err1 = abs(ExD + ExI);
                  Real_wp err2 = abs(EyD + EyI);

                  Err1 = max(Err1, err1); Err2 = max(Err2, err2);
                }
              else if ((abs(phase1) > 1e8) || (abs(phase2) > 1e8) || one_layer)
                {
                  TinyMatrix<Complex_q, General, 4, 4> A, B;
                  Complex_q kzp = nup_z(i*nb_points_y + j);
                  Complex_q cxp = -kx/kzp, cyp = -ky/kzp;          
                  Complex_q cx = -kx/kz, cy = -ky/kz;
                  
                  A(0, 0) = 1.0; A(0, 2) = -1.0;
                  A(1, 1) = 1.0; A(1, 3) = -1.0;
                  A(2, 0) = ky*cxp/eps_m; A(2, 1) = (ky*cyp-kzp)/eps_m;
                  A(2, 2) = -ky*cx/eps0; A(2, 3) = (-ky*cy+kz)/eps0;
                  A(3, 0) = (kzp-cxp*kx)/eps_m; A(3, 1) = -kx*cyp/eps_m;
                  A(3, 2) = (cx*kx-kz)/eps0; A(3, 3) = cy*kx/eps0;
                  
                  TinyVector<Complex_q, 4> b, x;
                  
                  b(0) = HxI; b(1) = HyI;
                  b(2) = ExI;
                  b(3) = EyI;

                  //DISP(A); DISP(b);
                  //B = A;
                  //GetInverse(A);
                  //Mlt(A, b, x);

                  // resolution avec schur
                  TinyMatrix<Complex_q, General, 2, 2> schur, schur0;
                  TinyVector<Complex_q, 2> bs, xs;
                  schur(0, 0) = A(2, 0) + A(2, 2); schur(0, 1) = A(2, 1) + A(2, 3);
                  schur(1, 0) = A(3, 0) + A(3, 2); schur(1, 1) = A(3, 1) + A(3, 3);

                  bs(0) = b(2) - A(2, 0)*b(0) - A(2, 1)*b(1);
                  bs(1) = b(3) - A(3, 0)*b(0) - A(3, 1)*b(1);
                  schur0 = schur;
                  Complex_q invDet = Real_q(1) / (schur(0, 0)*schur(1, 1) - schur(0, 1)*schur(1, 0));
                  schur(0, 0) = schur0(1, 1) *invDet; schur(0, 1) = -schur0(0, 1)*invDet;
                  schur(1, 0) = -schur0(1, 0) *invDet; schur(1, 1) = schur0(0, 0)*invDet;
                  Mlt(schur, bs, xs);

                  x(2) = xs(0); x(3) = xs(1);
                  x(0) = x(2) + b(0); x(1) = x(3) + b(1);
                  
                  HxD = x(2);
                  HyD = x(3);
                  HzD = (-kx*HxD-ky*HyD)/kz;
                  if ((j == 1) && (i == 1))
                    {
                      DISP(kzp);
                      DISP(HzD); DISP(HzD/HzI);
                    }
                  //DISP(HzD);
                  Complex_q ExD = (ky*HzD - kz*HyD) / eps0;
                  Complex_q EyD = (kz*HxD - kx*HzD) / eps0;

                  Complex_q Hx1 = x(0), Hy1 = x(1);
                  Complex_q Hz1 = (-kx*Hx1 - ky*Hy1) / kz1;

                  Complex_q Ex1 = (ky*Hz1 - kz1*Hy1) / eps_m;
                  Complex_q Ey1 = (kz1*Hx1 - kx*Hz1) / eps_m;

                  Real_wp err1 = abs(HxD + HxI - Hx1);
                  Real_wp err2 = abs(HyD + HyI - Hy1);
                  Real_wp err3 = abs(ExD + ExI - Ex1);
                  Real_wp err4 = abs(EyD + EyI - Ey1);

                  /*if (abs(err3) > 1e-5)
                    {
                      DISP(A); DISP(schur0); DISP(kx); DISP(ky); DISP(kz); DISP(kz1); DISP(bs);
                      DISP(x); DISP(b); DISP(schur); DISP(err1); DISP(err2); DISP(err3); DISP(err4);
                      exit(0);
                      }*/
                  
                  Err1 = max(Err1, err1); Err2 = max(Err2, err2);
                  Err3 = max(Err3, err3); Err4 = max(Err4, err4);
                }
              else
                {

                  // cas vectoriel (slab)
                  // inconnues : Hx,  Hy, Hx1, Hy1, Hx2, Hy2, Hx3, Hy3
                  Matrix<Complex_q> A(8, 8);
                  Vector<Complex_q> b(8);
                  
                  A.Zero(); b.Zero();
                  
                  Complex_q cx = -kx/kz, cy = -ky/kz;
                  Complex_q cx1 = -kx/kz1, cy1 = -ky/kz1;
                  Complex_q cx2 = -kx/kz2, cy2 = -ky/kz2;
                  Complex_q cx3 = -kx/kz3, cy3 = -ky/kz3;
                  
                  // premiere ligne : Hx + Hx^inc = Hx1 + Hx2
                  A(0, 0) = 1.0; A(0, 2) = -1.0; A(0, 4) = -1.0; b(0) = -HxI;
                  
                  // seconde ligne : Hy + Hy^inc = Hy1 + Hy2
                  A(1, 1) = 1.0; A(1, 3) = -1.0; A(1, 5) = -1.0; b(1) = -HyI;
                  
                  // troisieme ligne (ky Hz - kz Hy)/eps + ExI = Ex1 + Ex2
                  A(2, 0) = ky*cx/eps0; A(2, 1) = (ky*cy-kz) / eps0;
                  A(2, 2) = -ky*cx1/eps_m; A(2, 3) = -(ky*cy1 - kz1) / eps_m;
                  A(2, 4) = -ky*cx2/eps_m; A(2, 5) = -(ky*cy2 - kz2) / eps_m;
                  b(2) = -ExI;
                  
                  // quatrieme ligne (kz Hx - kx Hz)/eps + EyI = Ey1 + Ey2
                  A(3, 0) = (kz-kx*cx)/eps0; A(3, 1) = -kx*cy / eps0;
                  A(3, 2) = -(kz1-kx*cx1)/eps_m; A(3, 3) = kx*cy1 / eps_m;
                  A(3, 4) = -(kz2-kx*cx2)/eps_m; A(3, 5) = kx*cy2 / eps_m;
                  b(3) = -EyI;
                  
                  // cinquieme ligne : (Hx1 + Hx2) phase = Hx3
                  A(4, 2) = phase1; A(4, 4) = phase2; A(4, 6) = -1.0;
                  
                  // sixieme ligne : (Hy1 + Hy2) phase = Hy3
                  A(5, 3) = phase1; A(5, 5) = phase2; A(5, 7) = -1.0;
                  
                  // septieme ligne : (Ex1 + Ex2) phase = Ex3
                  A(6, 2) = ky*cx1*phase1/eps_m; A(6, 3) = (ky*cy1-kz1)*phase1 / eps_m;
                  A(6, 4) = ky*cx2*phase2/eps_m; A(6, 5) = (ky*cy2 - kz2)*phase2 / eps_m;
                  A(6, 6) = -ky*cx3/eps_ext; A(6, 7) = -(ky*cy3 - kz3) / eps_ext;
                  
                  // huitieme ligne : (Ey1 + Ey2) phase = Ey3
                  A(7, 2) = (kz1-kx*cx1)*phase1/eps_m; A(7, 3) = -kx*cy1*phase1 / eps_m;
                  A(7, 4) = (kz2-kx*cx2)*phase2/eps_m; A(7, 5) = -kx*cy2*phase2 / eps_m;
                  A(7, 6) = -(kz3-kx*cx3)/eps_ext; A(7, 7) = kx*cy3 / eps_ext;
                  
                  Vector<int> pivot(8);
                  GetLU(A, pivot);
                  
                  Vector<Complex_q> x(b);
                  SolveLU(A, pivot, x);
                  
                  HxD = x(0); HyD = x(1);
                  HzD = (-kx*HxD-ky*HyD)/kz;
                  
                  Complex_q Hx1 = x(2), Hy1 = x(3), Hx2 = x(4), Hy2 = x(5), Hx3 = x(6), Hy3 = x(7);
                  
                  Complex_q Hz1 = (-kx*Hx1 - ky*Hy1) / kz1;
                  Complex_q Hz2 = (-kx*Hx2 - ky*Hy2) / kz2;
                  Complex_q Hz3 = (-kx*Hx3 - ky*Hy3) / kz3;
                  
                  Complex_q ExD = (ky*HzD - kz*HyD) / eps0;
                  Complex_q EyD = (kz*HxD - kx*HzD) / eps0;
                  
                  Complex_q Ex1 = (ky*Hz1 - kz1*Hy1) / eps_m;
                  Complex_q Ey1 = (kz1*Hx1 - kx*Hz1) / eps_m;
                  
                  Complex_q Ex2 = (ky*Hz2 - kz2*Hy2) / eps_m;
                  Complex_q Ey2 = (kz2*Hx2 - kx*Hz2) / eps_m;
                  
                  Complex_q Ex3 = (ky*Hz3 - kz3*Hy3) / eps_ext;
                  Complex_q Ey3 = (kz3*Hx3 - kx*Hz3) / eps_ext;
                  
                  Real_wp err1 = abs(HxD + HxI - (Hx1 + Hx2));
                  Real_wp err2 = abs(HyD + HyI - (Hy1 + Hy2));
                  Real_wp err3 = abs(ExD + ExI - (Ex1 + Ex2));
                  Real_wp err4 = abs(EyD + EyI - (Ey1 + Ey2));
                  
                  Real_wp err5 = abs(Hx1*phase1 + Hx2*phase2 - Hx3);
                  Real_wp err6 = abs(Hy1*phase1 + Hy2*phase2 - Hy3);
                  Real_wp err7 = abs(Ex1*phase1 + Ex2*phase2 - Ex3);
                  Real_wp err8 = abs(Ey1*phase1 + Ey2*phase2 - Ey3);
                  
                  Err1 = max(Err1, err1); Err2 = max(Err2, err2);
                  Err3 = max(Err3, err3); Err4 = max(Err4, err4);
                  Err5 = max(Err5, err5); Err6 = max(Err6, err6);
                  Err7 = max(Err7, err7); Err8 = max(Err8, err8);
                }
            }

          if (j == 0)
            {
              if (abs(HzI) > 1e-30)
                ratio(i) = HzD / HzI;
              else
                ratio(i) = 1e-30;
              
              kx_over_k(i) = kx;
            }
          
          Hx(i*nb_points_y + j) = HxD * phase;
          Hy(i*nb_points_y + j) = HyD * phase;
          Hz(i*nb_points_y + j) = HzD * phase;

          //DISP(Hz(i*nb_points_y+j));
        }

    //DISP(ratio);
    //ratio.Write("ratio.dat");
    //kx_over_k.WriteText("kx.dat");
    
    //DISP(nb_points);
    //DISP(Err1); DISP(Err2); DISP(Err3); DISP(Err4);
    //DISP(Err5); DISP(Err6); DISP(Err7); DISP(Err8);
    
    FillPeriodic(nb_points_x, nb_points_y, Hz, Hz_grid);
    //WriteMatlab(Hz_grid, grid, "HzChap.dat", precision_output_file);
    
    // on revient dans le domaine reel
    fft2d.ApplyInverse(Hx);
    fft2d.ApplyInverse(Hy);
    fft2d.ApplyInverse(Hz);

    //ComputeGreenFieldH(Points2D, h1, kwave, fieldH);
    ComputeFieldH(Points2D, pt_wire, I0_wire, h1, fieldH);    
    
    // champ diffracte
    FillPeriodic(nb_points_x, nb_points_y, Hz, Hz_grid);
    //WriteMatlab(Hz_grid, grid, "Hz1.dat", precision_output_file);
    
    // on rajoute le champ incident
    for (int i = 0; i < Hz.GetM(); i++)
      {
        Hx(i) += fieldH(i)(0);
        Hy(i) += fieldH(i)(1);
        Hz(i) += fieldH(i)(2);
        //R3 pt(Points2D(i)(0), Points2D(i)(1), h1-1.0);
        //Real_wp r = Norm2(pt);
        //Hz(i) += exp(Iwp*kwave*r) / (4.0*pi_wp*r);
      }

    //DISP(Hz);
    
    // champ total
    FillPeriodic(nb_points_x, nb_points_y, Hx, Hx_grid);
    FillPeriodic(nb_points_x, nb_points_y, Hy, Hy_grid);
    FillPeriodic(nb_points_x, nb_points_y, Hz, Hz_grid);
    //WriteMatlab(Hz_grid, grid, "HzF.dat", precision_output_file);
  }

  void ComputeSlopeInterp(const Matrix<Complex_wp>& Hi, const Real_wp& xc, 
                          Vector<SplineInterpolation<Complex_wp> >& spline_slope)
  {
    Real_wp ymax = Longueur_interp, xmax = xc + Longueur_interp;
    Real_wp py_max = log(1.0 + ymax);
    Real_wp p1_max = log(xc+1.0), p2_max = log(1.0 + (xmax-xc));

    int N = nb_pts_interpolation(0);
    int N1 = nb_pts_interpolation(1), N2 = nb_pts_interpolation(2);
    
    Real_wp dy = py_max / (N-1);
    Real_wp dx1 = p1_max / (N1-1), dx2 = p2_max / (N2-1);

    int order = 4; VectReal_wp pts;
    Globatto<Real_wp> lob;
    pts.Reallocate(order+1); pts.Fill(); Mlt(Real_wp(1)/order, pts);
    lob.AffectPoints(pts); lob.ComputeGradPhi();
    
    VectComplex_wp SlopeCoefY; VectReal_wp xdiv(N1);
    SlopeCoefY.Reallocate(N1);
    for (int i = 0; i < N1; i++)
      {
        Complex_wp vloc = 0;
        for (int j = 0; j <= order; j++)
          vloc += lob.GradPhi(j, order)*Hi(i, N-order-1+j);
        
        Complex_wp SlopeY = vloc / (order*dy);
        xdiv(i) = -(N1-1-i)*dx1; SlopeCoefY(i) = SlopeY;
      }

    //Hi.Write("Hi.dat");
    //DISP(xdiv); DISP(SlopeCoefY);

    VectReal_wp teta_div; VectComplex_wp SlopeCoef, radius_div, Hz_div;    
    for (int i = 0; i <= N-2; i++)
      {
        Complex_wp vloc = 0;
        for (int j = 0; j <= order; j++)
          vloc += lob.GradPhi(j, order)*Hi(N1-2+N2-order+j, i);
        
        Complex_wp SlopeX = vloc / (order*dx2);
        int iymin = i-order/2, iymax = i+order/2;
        if (iymin < 0)
          {
            iymin = 0;
            iymax = order;
          }
        
        if (iymax >= N)
          {
            iymin = N-1-order;
            iymax = N-1;
          }
        
        vloc = Complex_wp(0, 0);
        for (int j = 0; j <= order; j++)
          vloc += lob.GradPhi(j, i-iymin)*Hi(N1+N2-2, iymin+j);
        
        Complex_wp SlopeY = vloc / (order*dy);
        Real_wp y = i*dy, x = p2_max;
        Real_wp teta = atan(y / x);

        Real_wp ux = cos(teta), uy = sin(teta);
        Complex_wp Slope = ux*SlopeX + uy*SlopeY;
        teta_div.PushBack(teta); SlopeCoef.PushBack(Slope);
        radius_div.PushBack(sqrt(x*x + y*y));
        Hz_div.PushBack(Hi(N1+N2-2, i));
        //DISP(i); DISP(SlopeY); DISP(SlopeX); DISP(teta); DISP(Slope);
      }

    for (int i = N2-1; i >= 0; i--)
      {
        Complex_wp vloc = 0;
        for (int j = 0; j <= order; j++)
          vloc += lob.GradPhi(j, order)*Hi(N1-1+i, N-order-1+j);
        
        Complex_wp SlopeY = vloc / (order*dy);
        int ixmin = i-order/2, ixmax = i+order/2;
        if (ixmin < 0)
          {
            ixmin = 0;
            ixmax = order;
          }
        
        if (ixmax >= N2)
          {
            ixmin = N2-1-order;
            ixmax = N2-1;
          }
        
        vloc = Complex_wp(0, 0);
        for (int j = 0; j <= order; j++)
          vloc += lob.GradPhi(j, i-ixmin)*Hi(N1-1+ixmin+j, N-1);
        
        Complex_wp SlopeX = vloc / (order*dx2);
        Real_wp y = py_max, x = i*dx2;
        Real_wp teta = 0.5*pi_wp;
        if (i > 0)
          teta = atan(y / x);

        Real_wp ux = cos(teta), uy = sin(teta);
        Complex_wp Slope = ux*SlopeX + uy*SlopeY;
        teta_div.PushBack(teta); SlopeCoef.PushBack(Slope);
        radius_div.PushBack(sqrt(x*x + y*y));
        Hz_div.PushBack(Hi(N1-1+i, N-1));
        //DISP(i); DISP(SlopeY); DISP(SlopeX); DISP(teta); DISP(Slope);
      }

    //DISP(teta_div); DISP(SlopeCoef); DISP(radius_div); DISP(Hz_div);
    spline_slope.Reallocate(4);
    spline_slope(0).Init(xdiv, SlopeCoefY);
    spline_slope(1).Init(teta_div, SlopeCoef);
    spline_slope(2).Init(teta_div, radius_div);
    spline_slope(3).Init(teta_div, Hz_div);
  }

  void ComputeInterpolationField(const Vector<VectR2>& pt_wire, const VectComplex_wp& I0_wire,
                                 const Vector<Vector<int> >& numL, const Vector<Matrix<Complex_wp> >& matHi, const TinyVector<bool, 2>& sign_Hi,
                                 const VectR2& Points2D, VectComplex_wp& Hz)
  {
    /* VectR2 Points2D(4);
    Points2D(0).Init(-2.8, 7.8);
    Points2D(1).Init(2.8, 7.8);
    Points2D(2).Init(-2.8, -7.8);
    Points2D(3).Init(2.8, -7.8); */
    
    Hz.Reallocate(Points2D.GetM());
    Hz.Zero();

    int order = 4; VectReal_wp pts; Globatto<Real_wp> lob;
    pts.Reallocate(order+1); pts.Fill(); Mlt(Real_wp(1)/order, pts);
    lob.AffectPoints(pts); //DISP(sign_Hi); 
    
    // boucle sur les fils
    R2 ur, vr, vec_u;
    int N = nb_pts_interpolation(0);
    int N1 = nb_pts_interpolation(1), N2 = nb_pts_interpolation(2);
    for (int n = 0; n < pt_wire.GetM(); n++)
      for (int i = 0; i < pt_wire(n).GetM()-1; i++)
        {
          R2 ptM = 0.5*(pt_wire(n)(i) + pt_wire(n)(i+1));
          vec_u(0) = pt_wire(n)(i+1)(0) - pt_wire(n)(i)(0);
          vec_u(1) = pt_wire(n)(i+1)(1) - pt_wire(n)(i)(1);
          
          Real_wp L = Norm2(vec_u); Real_wp xc = 0.5*L;
          Real_wp invNorm = 1.0/L;
          int num_l = numL(n)(i);
	  ur.Init(invNorm*vec_u(0), invNorm*vec_u(1));
	  vr.Init(ur(1), -ur(0));

          Real_wp ymax = Longueur_interp, xmax = xc + Longueur_interp;
          Real_wp py_max = log(1.0 + ymax);
          Real_wp p1_max = log(xc+1.0), p2_max = log(1.0 + (xmax-xc));
          Real_wp dy = py_max / (N-1);
          Real_wp dx1 = p1_max / (N1-1), dx2 = p2_max / (N2-1);
          
          // calcul des pentes
          Vector<SplineInterpolation<Complex_wp> > spline_slope;
          ComputeSlopeInterp(matHi(num_l), xc, spline_slope);
          
          // boucle sur les points
          for (int k = 0; k < Points2D.GetM(); k++)
            {
              // on calcule x, y dans le repere du fil considere
              R2 diff = Points2D(k) - ptM;
	      Real_wp x = abs(DotProd(diff, ur));
	      Real_wp y = DotProd(diff, vr);
              
              Real_wp py = log(1.0 + abs(y));
              Real_wp px = 0;
              if (x < xc)
                px = -log(1.0 + abs(x-xc));
              else
                px = log(1.0 + abs(x-xc));
              
              // on recupere les bornes iymin, iymax
              int iy = toInteger(round(py / py_max*N));
              int iymin = iy -order/2, iymax = iy+order/2;
              if (iymin < 0)
                {
                  iymax = order;
                  iymin = 0;
                }
              
              if (iymax >= N)
                {
                  iymin = N-1-order;
                  iymax = N-1;
                }
              
              Real_wp yloc = Real_wp(py - iymin*dy) / (order*dy);
              yloc = min(1.0, yloc); yloc = max(0.0, yloc);
              
              // calcul des fonctions de base en y
              VectReal_wp phi_y, phi_x;
              lob.ComputeValuesPhiRef(yloc, phi_y);
              
              // calcul des fonctions de base en x
              int ixmin = 0, ixmax = order, ix; Real_wp xloc;
              if (x < xc)
                {
                  ix = N1 - toInteger(round(-px/p1_max*N1));
                  ixmin = ix - order/2; ixmax = ix + order/2;
                  if (ixmin < 0)
                    {
                      ixmax = order;
                      ixmin = 0;
                    }
                  
                  if (ixmax >= N1)
                    {
                      ixmin = N1-1-order;
                      ixmax = N1-1;
                    }
                  
                  xloc = (px + p1_max - ixmin*dx1) / (order*dx1);
                  xloc = min(1.0, xloc); xloc = max(0.0, xloc);
                  
                  lob.ComputeValuesPhiRef(xloc, phi_x);
                }
              else
                {
                  ix = toInteger(round(px/p2_max*N2));
                  ixmin = ix - order/2; ixmax = ix + order/2;
                  if (ixmin < 0)
                    {
                      ixmax = order;
                      ixmin = 0;
                    }
                  
                  if (ixmax >= N2)
                    {
                      ixmin = N2-1-order;
                      ixmax = N2-1;
                    }
                  
                  xloc = (px - ixmin*dx2) / (order*dx2);
                  xloc = min(1.0, xloc); xloc = max(0.0, xloc);
                  
                  lob.ComputeValuesPhiRef(xloc, phi_x);
                  ixmin += N1-1; ixmax += N1-1;
                }
              
              Complex_wp val(0, 0);
              if ((py > py_max) || (px > p2_max))
                {
                  // cas ou on est en dehors du domaine d'interpolation
                  // il faut trouver la pente
                  //DISP(py); DISP(px); DISP(py_max); DISP(p2_max);
                  //DISP(xloc); DISP(yloc);
                  
                  if (px <= 0)
                    {
                      Complex_wp coef = spline_slope(0).Evaluate(px);
                      for (ix = ixmin; ix <= ixmax; ix++)
                        val += phi_x(ix-ixmin)*matHi(num_l)(ix, N-1);
                      
                      val += coef*(py - py_max);
                    }
                  else
                    {
                      Real_wp teta = atan(py / px); Real_wp r = sqrt(px*px + py*py);
                      Complex_wp coef = spline_slope(1).Evaluate(teta);
                      Real_wp rmax = real(spline_slope(2).Evaluate(teta));
                      val = spline_slope(3).Evaluate(teta);
                      val += coef*(r - rmax);
                    }
                }
              else
                {
                  // cas ou on est a l'interieur du domaine d'interpolation
                  // on effectue l'interpolation                  
                  for (ix = ixmin; ix <= ixmax; ix++)
                    for (iy = iymin; iy <= iymax; iy++)
                      val += phi_x(ix-ixmin)*phi_y(iy-iymin)*matHi(num_l)(ix, iy);
                }

              Real_wp vr = exp(real(val));
              if (!sign_Hi(0))
                vr = -vr;
              
              Real_wp vi = exp(imag(val));
              if (!sign_Hi(1))
                vi = -vi;
              
              //DISP(x); DISP(y); DISP(px); DISP(py); DISP(phi_x); DISP(phi_y);
              //DISP(ixmin); DISP(ixmax); DISP(iymin); DISP(iymax);
              Complex_wp vloc = y/Complex_wp(vr, vi);              
              Hz(k) += vloc*I0_wire(n);
            }
        }

    //DISP(Hz);
  }
  
  void TestFourier()
  {
    if (type_algo_fourier == DIRECT)
      {
        int nb_points_x = nb_points_fourier, nb_points_y = nb_points_fourier;
        Real_wp xmin(xmin_fourier), xmax(xmax_fourier), ymin(ymin_fourier), ymax(ymax_fourier);
        
        VectReal_wp PointsX, PointsY; VectComplex_wp Hx_grid, Hy_grid, Hz_grid;
        ComputeFieldFourier(liste_wire, current_wire, xmin, xmax, ymin, ymax, nb_points_x,
                            nb_points_y, PointsX, PointsY, Hx_grid, Hy_grid, Hz_grid);

        VectReal_wp val(Hz_grid.GetM());
        for (int k = 0; k < Hz_grid.GetM(); k++)
          val(k) = 20.0*log10(abs(Hz_grid(k)/level_field_zeroDB)+epsilon);
        
        GridInterpolationFull<Dimension2> grid;
        grid.SetPlaneOutput(xmin, xmax, ymin, ymax, nb_points_x+1, nb_points_y+1);
        
        WriteMatlab(val, grid, "FieldH_Grid0_Z.dat", precision_output_file);
        
        for (int k = 0; k < Hz_grid.GetM(); k++)
          {
            Real_wp Ht = sqrt(absSquare(Hx_grid(k)) + absSquare(Hy_grid(k)));
            val(k) = 20.0*log10(abs(Ht/level_field_zeroDB)+epsilon);
          }
        
        WriteMatlab(val, grid, "FieldH_Grid0_XY.dat", precision_output_file);

      }
    else if (type_algo_fourier == INTERPOLATION)
      {
        // on calcule les differentes longueurs
        VectReal_wp tabL;
        for (int i = 0; i < liste_wire.GetM(); i++)
          for (int j = 0; j < liste_wire(i).GetM()-1; j++)
            tabL.PushBack(0.5*liste_wire(i)(j).Distance(liste_wire(i)(j+1)));
        
        Sort(tabL);
        
        VectReal_wp listeL(1);
        listeL(0) = tabL(0); 
        for (int i = 1; i < tabL.GetM(); i++)
          {
            if (abs(tabL(i) - tabL(i-1)) > 1e-6)
              listeL.PushBack(tabL(i));
          }
        
        Vector<Vector<int> > numL(liste_wire.GetM());
        for (int i = 0; i < liste_wire.GetM(); i++)
          {
            numL(i).Reallocate(liste_wire(i).GetM()-1);
            for (int j = 0; j < liste_wire(i).GetM()-1; j++)
              {
                Real_wp L = 0.5*liste_wire(i)(j).Distance(liste_wire(i)(j+1));
                for (int k = 0; k < listeL.GetM(); k++)
                  if (abs(listeL(k) - L) <= 2e-6)
                    {
                      numL(i)(j) = k;
                      break;
                    }
              }
          }
        
        DISP(listeL); DISP(numL);
        VectReal_wp Lcomp; int Ns = nb_pts_interpolation(3);
        if (!interp_spline_L)
          Lcomp = listeL;
        else
          {            
            Lcomp.Reallocate(Ns);
            Real_wp lmin = -2.0, lmax = log(listeL(listeL.GetM()-1));
            Real_wp dl = (lmax - lmin)/(Ns-1); dl = 0.65;
            for (int num_l = 0; num_l < Ns; num_l++)
              Lcomp(num_l) = exp(lmin + num_l*dl);
          }

        DISP(Lcomp);
        // on calcule les coefficients sur les points d'interpolation
        Vector<Matrix<Complex_wp> > matHi(Lcomp.GetM());
        TinyVector<bool, 2> sign_Hi;
        for (int num_l = 0; num_l < Lcomp.GetM(); num_l++)
          {
            Real_wp L = Lcomp(num_l);
            Vector<VectR2> elem_wire(1); elem_wire(0).Reallocate(2);
            VectComplex_wp I0_wire(1); I0_wire(0) = 1.0;
            elem_wire(0)(0).Init(-L, 0);
            elem_wire(0)(1).Init(L, 0);
            
            int nb_points_x = nb_points_fourier, nb_points_y = nb_points_fourier;
            Real_wp xmin(xmin_fourier), xmax(xmax_fourier), ymin(ymin_fourier), ymax(ymax_fourier);
            xmin -= L; xmax += L;
            
            VectReal_wp PointsX, PointsY; VectComplex_wp Hx_grid, Hy_grid, Hz_grid;
            ComputeFieldFourier(elem_wire, I0_wire, xmin, xmax, ymin, ymax, nb_points_x,
                                nb_points_y, PointsX, PointsY, Hx_grid, Hy_grid, Hz_grid);
            
            FindPolynomialCoef(Hz_grid, PointsX, PointsY, L, matHi(num_l), sign_Hi);
          }
        
        // interpolation par spline cubique
        if (interp_spline_L)
          {
            Vector<Matrix<Complex_wp> > matHdiv(matHi);
            matHi.Reallocate(listeL.GetM());
            for (int n = 0; n < listeL.GetM(); n++)
              matHi(n).Reallocate(matHdiv(0).GetM(), matHdiv(0).GetN());
            
            for (int j = 0; j < matHdiv(0).GetM(); j++)
              for (int k = 0; k < matHdiv(0).GetN(); k++)
                {
                  SplineInterpolation<Complex_wp> spline;
                  VectComplex_wp ydiv(Ns);
                  for (int i = 0; i < Ns; i++)
                    ydiv(i) = matHdiv(i)(j, k);
                  
                  spline.Init(Lcomp, ydiv);
                  for (int n = 0; n < listeL.GetM(); n++)
                    matHi(n)(j, k) = spline.Evaluate(listeL(n));
                }
          }

        if (type_algo_min == TARGET_POINT)
          {
            VectR2 vec_target(1); vec_target(0) = target;
            VectComplex_wp Hz;
            ComputeInterpolationField(liste_wire, current_wire, numL, matHi,
                                      sign_Hi, vec_target, Hz);

            Real_wp min_val = abs(Hz(0));
            cout << "Valeur minimale de Hz en A/m = " << min_val << endl;
            cout << "Atteinte au point " << target << endl;
            
            Real_wp I = level_field_zeroDB / min_val * pow(10.0, level_dB/20.0);
            cout << "Pour atteindre " << level_dB << " decibels, il faut multiplier I par " << I << endl;            
          }

        // calcul final
        VectR2 Points2D; VectReal_wp TetaPoints;
        for (int num = 0; num < var_grid.GetM(); num++)
          {
            Points2D.Clear();
            var_grid(num).GenerateGridPoints(Points2D, TetaPoints);
            
            VectComplex_wp Hz;
            ComputeInterpolationField(liste_wire, current_wire, numL, matHi,
                                      sign_Hi, Points2D, Hz);

            VectReal_wp val(Points2D.GetM());
	    for (int k = 0; k < Points2D.GetM(); k++)
              val(k) = 20.0*log10(abs(Hz(k)/level_field_zeroDB)+epsilon);
            
	    WriteMatlab(val, var_grid(num), "FieldH_Grid" +to_str(num) + "_Z.dat", precision_output_file);
          }
        
      }
    else if (type_algo_fourier == TEST_ERROR)
      {
        Real_wp lmin = -2.0, lmax = 4.0, dl = (lmax - lmin) / 100;
        Vector<Matrix<Complex_wp> > matHi(101);
        Real_wp err_max = 0;
        for (int num_l = 0; num_l <= 100; num_l++)
          {
            Real_wp L = exp(lmin + num_l*dl);
            Vector<VectR2> elem_wire(1); elem_wire(0).Reallocate(2);
            VectComplex_wp I0_wire(1); I0_wire(0) = 1.0;
            elem_wire(0)(0).Init(-L, 0);
            elem_wire(0)(1).Init(L, 0);
            //DISP(L); DISP(liste_wire(0)(0)); DISP(liste_wire(0)(1));
            
            int nb_points_x = nb_points_fourier, nb_points_y = nb_points_fourier;
            Real_wp xmin(xmin_fourier), xmax(xmax_fourier), ymin(ymin_fourier), ymax(ymax_fourier);
            xmin -= L; xmax += L;
            
            VectReal_wp PointsX, PointsY; VectComplex_wp Hx_grid, Hy_grid, Hz_grid;
            ComputeFieldFourier(elem_wire, I0_wire, xmin, xmax, ymin, ymax, nb_points_x,
                                nb_points_y, PointsX, PointsY, Hx_grid, Hy_grid, Hz_grid);
            
            TinyVector<bool, 2> sign_Hi;
            Real_wp err = FindPolynomialCoef(Hz_grid, PointsX, PointsY, L, matHi(num_l), sign_Hi);        
            err_max = max(err_max, err);
          }
        
        DISP(err_max);
        
        // interpolation par spline cubique
        int Ns = nb_pts_interpolation(3); int dNs = 100/(Ns-1);
        VectReal_wp xdiv(Ns); VectComplex_wp ydiv(Ns);
        for (int j = 0; j < matHi(0).GetM(); j++)
          for (int k = 0; k < matHi(0).GetN(); k++)
            {
              SplineInterpolation<Complex_wp> spline;
              for (int i = 0; i < Ns; i++)
                {
                  int num_l = dNs*i;
                  if (i == Ns-1)
                    num_l = 100;
                  
                  xdiv(i) = lmin + num_l*dl;
                  ydiv(i) = matHi(num_l)(j, k);
                }
              
              spline.Init(xdiv, ydiv); err_max = 0;
              for (int i = 0; i <= 100; i++)
                {
                  Complex_wp coef = spline.Evaluate(lmin + i*dl);
                  err_max = max(err_max, abs(coef - matHi(i)(j, k))/abs(coef));
                }
              
              DISP(j); DISP(k); DISP(err_max);
            }
      }
  }
  
  void RunAll()
  {
    if (!smart_loop)
      {
        DetectVertexWall();

        // calcul auracast
        // loop on output grids
        VectR2 Points2D; VectReal_wp TetaPoints;
        for (int num = 0; num < var_grid.GetM(); num++)
          {
            Points2D.Clear();
            var_grid(num).GenerateGridPoints(Points2D, TetaPoints);
            
            VectReal_wp val(Points2D.GetM());
            val.Zero();
            
            ComputeFieldAura(Points2D, val);
            
            WriteMatlab(val, var_grid(num), "FieldA_Grid" +to_str(num) + ".dat", precision_output_file);
          }
        
        return;
      }
    
    // writing the lines in a text file
    for (int n = 0; n < liste_wire.GetM(); n++)
      {
	string nom_fichier = "poly" + to_str(n) + ".dat";
	ofstream file_out(nom_fichier.data());
	file_out.precision(15);
	for (int i = 0; i < liste_wire(n).GetM(); i++)
	  file_out << liste_wire(n)(i)(0) << " " << liste_wire(n)(i)(1) << '\n';

	file_out.close();
      }

    if (rho_metal > 0)
      {
        TestFourier();
        return;
      }
    
    // loop on output grids
    VectR2 Points2D; VectReal_wp TetaPoints;
    for (int num = 0; num < var_grid.GetM(); num++)
      {
	Points2D.Clear();
	var_grid(num).GenerateGridPoints(Points2D, TetaPoints);
	
	VectReal_wp val(Points2D.GetM());
	val.Zero();
	
	Vector<TinyVector<Complex_wp, 3> > fieldH(Points2D.GetM());	
	ComputeFieldH(Points2D, liste_wire, current_wire, h, fieldH);

	if (choice_field_output == -1)
	  {
	    for (int k = 0; k < Points2D.GetM(); k++)
	      {
                Real_wp Ht = sqrt(absSquare(fieldH(k)(0)) + absSquare(fieldH(k)(1)));
                val(k) = 20.0*log10(Ht/level_field_zeroDB+epsilon);
              }
            
	    WriteMatlab(val, var_grid(num), "FieldH_Grid" +to_str(num) + "_XY.dat", precision_output_file);
            
	    /*for (int k = 0; k < Points2D.GetM(); k++)
	      val(k) = 20.0*log10(abs(fieldH(k)(1)/level_field_zeroDB)+epsilon);
	    
              WriteMatlab(val, var_grid(num), "FieldH_Grid" +to_str(num) + "_Y.dat", precision_output_file); */
	    
	    for (int k = 0; k < Points2D.GetM(); k++)
	      val(k) = 20.0*log10(abs(fieldH(k)(2)/level_field_zeroDB)+epsilon);
	    
	    WriteMatlab(val, var_grid(num), "FieldH_Grid" +to_str(num) + "_Z.dat", precision_output_file);
	  }
        else if (choice_field_output == -2)
	  {
	    for (int k = 0; k < Points2D.GetM(); k++)
	      val(k) = realpart(fieldH(k)(2));
            
	    WriteMatlab(val, var_grid(num), "FieldH_Grid" +to_str(num) + ".dat", precision_output_file);
	  }
	else
	  {
	    for (int k = 0; k < Points2D.GetM(); k++)
	      val(k) = 20.0*log10(abs(fieldH(k)(choice_field_output)/level_field_zeroDB)+epsilon);

	    WriteMatlab(val, var_grid(num), "FieldH_Grid" +to_str(num) + ".dat", precision_output_file);
	  }
      }

    // finding the minimum value of Hz
    Points2D.Clear();
    if (type_algo_min == AUTO)
      {
	// mean value of points in liste_wire
	int num = 0; R2 ptCenter;
	for (int n = 0; n < liste_wire.GetM(); n++)
	  for (int i = 0; i < liste_wire(n).GetM()-1; i++)
	    {
	      ptCenter += liste_wire(n)(i);
	      num++;
	    }

	ptCenter *= Real_wp(1) / num;
	cout << "Centre salle = " << ptCenter << endl;
	Points2D.PushBack(ptCenter);
      }
    else if (type_algo_min == TARGET_POINT)
      {
	Points2D.PushBack(target);
      }
    else if (type_algo_min == MINIMUM_BOX)
      {
	Points2D.Reallocate(nb_points_x*nb_points_y);
	Real_wp dx = (maxBox(0) - minBox(0)) / (nb_points_x-1);
	Real_wp dy = (maxBox(1) - minBox(1)) / (nb_points_y-1);
	for (int i = 0; i < nb_points_x; i++)
	  for (int j = 0; j < nb_points_y; j++)
	    Points2D(i*nb_points_y + j).Init(minBox(0) + i*dx, minBox(1) + j*dy);
      }
    
    Vector<TinyVector<Complex_wp, 3> > fieldH(Points2D.GetM());	
    ComputeFieldH(Points2D, liste_wire, current_wire, h, fieldH);

    Real_wp min_val(1e300); int arg_min = -1;
    for (int i = 0; i < fieldH.GetM(); i++)
      {
	if (abs(fieldH(i)(2)) < min_val)
	  {
	    min_val = abs(fieldH(i)(2));
	    arg_min = i;
	  }
      }

    cout << "Valeur minimale de Hz en A/m = " << min_val << endl;
    cout << "Atteinte au point " << Points2D(arg_min) << endl;

    Real_wp I = level_field_zeroDB / min_val * pow(10.0, level_dB/20.0);
    cout << "Pour atteindre " << level_dB << " decibels, il faut multiplier I par " << I << endl;
  }
  
};
  
int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  if (argc != 2)
    {
      cout << "Donnez un fichier de donnees" << endl;
      abort();
    }

  string input_file(argv[1]);
  
  ComputationFieldH var;

  // on lit le fichier de donnees
  ReadInputFile(input_file, var);

  // on calcule H sur les points de sortie demandes
  var.RunAll();

  return FinalizeMontjoie();
}
