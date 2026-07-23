#define MONTJOIE_WITH_THREE_DIM
#define MONTJOIE_WITH_NODAL_H1

#include "Quadrature/MontjoieQuadrature.hxx"

using namespace Montjoie;

void CheckTriangleQuadrature(int p, const VectR2& Points2D, const VectReal_wp& Weights2D)
{
  // testing integral of x^i y^j
  for (int i = 0; i <= p; i++)
    for (int j = 0; j <= p-i; j++)
      {
        double val_exact = tgamma(i+1)*tgamma(j+1)/tgamma(i+j+3);
        double val = 0;
        for (int m = 0; m < Points2D.GetM(); m++)
          val += Weights2D(m)*pow(Points2D(m)(0), i)*pow(Points2D(m)(1), j);
        
        if (abs(val-val_exact)/val > 1e-13)
          {
            cout << " Order " << i+j << " not exactly integrated for order p = " << p <<endl;
            DISP(i); DISP(j); DISP(val); DISP(val_exact);
          }
      }  
}

void GenerateTriangleFormulas()
{
  ofstream file_out("triangle.cxx");
  file_out.precision(16);
  file_out<<"switch(p) \n { \n ";
  for (int p = 1; p <= 19; p++)
    {
      int nb_points = 0; double x, y, t, omega;
      Vector<double> Weights2D; Vector<R2> Points2D;
      
      string name_file("src/Quadrature/Table_Quadrature.txt");
      string title = string("Gauss  quadrature  points  and  weights  on  the  ")
			      + "reference  triangle  order  p=";
      title += to_str(p);
      
      // the file where the formulas are stored
      ifstream file_in(name_file.data());
      string ligne;
      
      if (!file_in.is_open())
	{
	  cout<<" File containing quadrature formulas not found"<<endl;
	  abort();
	}
      
      while (!file_in.eof())
	{
	  getline(file_in, ligne);
	  DeleteSpaceAtExtremityOfString(ligne);
	  
	  if (!ligne.compare(title))
	    {
	      // line has been found
	      file_in >> nb_points;
	      Points2D.Reallocate(nb_points);
	      Weights2D.Reallocate(nb_points);
	      for (int j = 0; j < nb_points; j++)
		{
		  file_in >> x >> y >> omega;
		  Points2D(j).Init(0.5*x + 0.5, 0.5*y + 0.5);
		  Weights2D(j) = 0.25*omega;
		}
	    }
	}
      
      CheckTriangleQuadrature(p, Points2D, Weights2D);
      
      // finding couples
      int nb_couples = 0;
      int nb_remaining_points = nb_points;
      Vector<bool> PointUsed(nb_points); PointUsed.Fill(false);
      Vector<R2> points2d(nb_points);
      Vector<double> weights2d(nb_points);
      TinyVector<R2, 6> sym_points;
      IVect NbPointsCouple(nb_points); NbPointsCouple.Fill(0);
      int n1 = 0, n3 = 0, n6 = 0;
      while (nb_remaining_points > 0)
	{
	  int k = 0;
	  // searching an unused point
	  for (int j = 0; j < nb_points; j++)
	    if (!PointUsed(j))
	      {
		PointUsed(j) = true;
		nb_remaining_points--;
		k++;
		// generating other points by symmetry
		x = Points2D(j)(0);
		y = Points2D(j)(1);
		t = 1.0 - x - y;
		sym_points(0).Init(x, y); 
		sym_points(1).Init(x, t);
		sym_points(2).Init(y, x); 
		sym_points(3).Init(y, t);
		sym_points(4).Init(t, x); 
		sym_points(5).Init(t, y);
		
		// storing this point
		if (abs(x-y) < 1e-13)
		  points2d(nb_couples).Init(x, t);
		else if (abs(x-t) < 1e-13)
		  points2d(nb_couples).Init(x, y);
		else
		  points2d(nb_couples).Init(y, x);
		
		weights2d(nb_couples) = Weights2D(j);
		break;
	      }
	  
	  for (int j = 0; j < nb_points; j++)
	    for (int m = 0; m < 6; m++)
	      if (!PointUsed(j))
		if (Points2D(j) == sym_points(m))
		  {
		    PointUsed(j) = true;
		    nb_remaining_points--;
		    k++;
		  }
	  
	  NbPointsCouple(nb_couples) = k;
	  
	  if (k == 1)
	    {
	      if (n1 == 1)
		{
		  abort();
		}
	      n1++;
	    }
	  else if (k == 3)
	    n3++;
	  else if (k == 6)
	    n6++;
	  
	  nb_couples++;	  	  
	}
      
      bool cpp_code = false;
      if (cpp_code)
	{
	  file_out << "case " << p << " : " << endl;
	  file_out <<" { \n n1 = " << n1 << "; n3 = "<<n3<<"; n6 = "<<n6<<";\n";
	  file_out << "nb_points = n1 + 3*n3 + 6*n6;" << endl;
	  file_out << "points2d.Reallocate(nb_points); \n weights2d.Reallocate(nb_points); \n";
	  int ind = 0;
	  for (int j = 0; j < nb_couples; j++)
	    if (NbPointsCouple(j) == 1)
	      {
		file_out << "points2d(0).Init(1.0/3, 1.0/3);\n weights2d(0) = " << weights2d(j) << ";\n";
		ind++;
	      }
	  
	  for (int j = 0; j < nb_couples; j++)
	    if (NbPointsCouple(j) == 3)
	      {
		file_out << "points2d("<<ind<<").Init("<< points2d(j)(0) <<", "<<points2d(j)(1);
		file_out << "); \n weights2d(" << ind << ") = " << weights2d(j) << ";\n";
		ind += 3;
	      }
	  
	  for (int j = 0; j < nb_couples; j++)
	    if (NbPointsCouple(j) == 6)
	      {
		file_out << "points2d("<<ind<<").Init("<< points2d(j)(0) <<", "<<points2d(j)(1);
		file_out << "); \n weights2d(" << ind << ") = " << weights2d(j) << ";\n";
		ind += 6;
	      }
	  
	  file_out << "} \n break; \n ";
	}
      else
	{
	  file_out << "case(" << p << ") " << endl;
	  file_out <<"\n n1 = " << n1 << "; n3 = "<<n3<<"; n6 = "<<n6<<";\n";
	  file_out << "nb_points = n1 + 3*n3 + 6*n6" << endl;
	  file_out << "call resize(points2d, nb_points) \n call resize(weights2d, nb_points) \n";
	  int ind = 1;
	  for (int j = 0; j < nb_couples; j++)
	    if (NbPointsCouple(j) == 1)
	      {
		file_out << "points2d(:, 1) = (/ 1.0_wp/3.0_wp, 1.0_wp/3.0_wp /) \n weights2d(1) = " << weights2d(j) << "_wp\n";
		ind++;
	      }
	  
	  for (int j = 0; j < nb_couples; j++)
	    if (NbPointsCouple(j) == 3)
	      {
		file_out << "points2d(:, "<<ind<<") = (/"<< points2d(j)(0) <<"_wp, "<<points2d(j)(1)<<"_wp /)\n";
		file_out << "weights2d(" << ind << ") = " << weights2d(j) << "_wp;\n";
		ind += 3;
	      }
	  
	  for (int j = 0; j < nb_couples; j++)
	    if (NbPointsCouple(j) == 6)
	      {
		file_out << "points2d(:, "<<ind<<") = (/"<< points2d(j)(0) <<"_wp, "<<points2d(j)(1)<<"_wp /)\n";
		file_out << "weights2d(" << ind << ") = " << weights2d(j) << "_wp;\n";
		ind += 6;
	      }
	  
	  file_out << "\n";
	}
    }
  
  file_out.close();

}


void CheckTetrahedronQuadrature(int p, const VectR3& Points3D, VectReal_wp& Weights3D)
{
  int nb_points = Points3D.GetM();
  // testing integral of x^i y^j z^k
  for (int i = 0; i <= p; i++)
    for (int j = 0; j <= p-i; j++)
      for (int k = 0; k <= p-i-j; k++)
        {
          double val_exact = tgamma(i+1)*tgamma(j+1)*tgamma(k+1)/tgamma(i+j+k+4);
          double val = 0;
          for (int m = 0; m < nb_points; m++)
            val += Weights3D(m)*pow(Points3D(m)(0), i)*pow(Points3D(m)(1), j)*pow(Points3D(m)(2), k);
          
          if (abs(val-val_exact)/val > 1e-10)
            {
              cout << " Order " << i+j+k << " not exactly integrated for order p = " << p <<endl;
              DISP(i); DISP(j); DISP(k); DISP(val); DISP(val_exact);
            }
        }
}

void GenerateTetrahedronFormulas()
{
  ofstream file_out("tetra.cxx");
  file_out.precision(16);
  file_out<<"switch(p) \n { \n ";
  for (int p = 1; p <= 21; p++)
    {
      int nb_points = 0; double x, y, z, t, omega;
      Vector<double> Weights3D; Vector<R3> Points3D;
      
      string name_file("src/Quadrature/Table_Quadrature.txt");
      string title = string("Gauss  quadrature  constants  for  the  ")
			      + "reference  tetrahedron  order  p=";
      title += to_str(p);
      
      // the file where the formulas are stored
      ifstream file_in(name_file.data());
      string ligne;
      
      if (!file_in.is_open())
	{
	  cout<<" File containing quadrature formulas not found"<<endl;
	  abort();
	}
      
      while (!file_in.eof())
	{
	  getline(file_in, ligne);
	  DeleteSpaceAtExtremityOfString(ligne);
	  
	  if (!ligne.compare(title))
	    {
	      // line has been found
	      file_in >> nb_points;
	      Points3D.Reallocate(nb_points);
	      Weights3D.Reallocate(nb_points);
	      for (int j = 0; j < nb_points; j++)
		{
		  file_in >> x >> y >> z >> omega;
		  Points3D(j).Init(0.5*x + 0.5, 0.5*y + 0.5, 0.5*z + 0.5);
		  Weights3D(j) = 0.125*omega;
		}
	    }
	}
      
      if (nb_points > 0)
	{
          CheckTetrahedronQuadrature(p, Points3D, Weights3D);

	  // finding couples
	  int nb_couples = 0;
	  int nb_remaining_points = nb_points;
	  Vector<bool> PointUsed(nb_points); PointUsed.Fill(false);
	  Vector<R3> points3d(nb_points);
	  Vector<double> weights3d(nb_points);
	  TinyVector<R3, 24> sym_points;
	  Vector<R3> other_points(30);
	  IVect num_points(30); 
	  IVect NbPointsCouple(nb_points); NbPointsCouple.Fill(0);
	  int n1 = 0, n4 = 0, n6 = 0, n12 = 0, n24 = 0;
	  while (nb_remaining_points > 0)
	    {
	      int k = 0, expected_orbit = 24;
	      // searching an unused point
	      for (int j = 0; j < nb_points; j++)
		if (!PointUsed(j))
		  {
		    PointUsed(j) = true;
		    nb_remaining_points--;
		    num_points(k) = j;
		    other_points(k) = Points3D(j);
		    k++;
		    // generating other points by symmetry
		    x = Points3D(j)(0);
		    y = Points3D(j)(1);
		    z = Points3D(j)(2);
		    t = 1.0 - x - y - z;
		    sym_points(0).Init(x, y, z); 
		    sym_points(1).Init(x, t, z);
		    sym_points(2).Init(y, x, z); 
		    sym_points(3).Init(y, t, z);
		    sym_points(4).Init(t, x, z); 
		    sym_points(5).Init(t, y, z);

		    sym_points(6).Init(x, y, t); 
		    sym_points(7).Init(x, z, t);
		    sym_points(8).Init(y, x, t); 
		    sym_points(9).Init(y, z, t);
		    sym_points(10).Init(z, x, t); 
		    sym_points(11).Init(z, y, t);

		    sym_points(12).Init(z, y, x); 
		    sym_points(13).Init(z, t, x);
		    sym_points(14).Init(y, z, x); 
		    sym_points(15).Init(y, t, x);
		    sym_points(16).Init(t, z, x); 
		    sym_points(17).Init(t, y, x);

		    sym_points(18).Init(x, z, y); 
		    sym_points(19).Init(x, t, y);
		    sym_points(20).Init(z, x, y); 
		    sym_points(21).Init(z, t, y);
		    sym_points(22).Init(t, x, y); 
		    sym_points(23).Init(t, z, y);
		 
		    // storing this point
		    if (abs(x-y) < 1e-13)
		      {
			if (abs(x-z) < 1e-13)
			  {
			    points3d(nb_couples).Init(x, x, t);
			    if (abs(x-t) < 1e-13)
			      expected_orbit = 1;
			    else
			      expected_orbit = 4;
			  }
			else
			  {
			    points3d(nb_couples).Init(x, x, z);
			    if (abs(x-t) < 1e-13)
			      expected_orbit = 4;
			    else if (abs(z-t) < 1e-13)
			      expected_orbit = 6;
			    else
			      expected_orbit = 12;
			  }
		      }
		    else if (abs(x-t) < 1e-13)
		      {
			if (abs(x-z) < 1e-13)
			  {
			    points3d(nb_couples).Init(x, x, y);
			    expected_orbit = 4;
			  }
			else
			  {
			    points3d(nb_couples).Init(x, x, z);
			    if (abs(y-z) < 1e-13)
			      expected_orbit = 6;
			    else
			      expected_orbit = 12;
			  }
		      }
		    else if (abs(x-z) < 1e-13)
		      {
			points3d(nb_couples).Init(x, x, y);
			if (abs(y-t) < 1e-13)
			  expected_orbit = 6;
			else
			  expected_orbit = 12;
		      }
		    else if (abs(y-z) < 1e-13)
		      {
			points3d(nb_couples).Init(y, y, x);
			if (abs(y-t) < 1e-13)
			  expected_orbit = 4;
			else
			  expected_orbit = 12;
		      }
		    else if (abs(y-t) < 1e-13)
		      {
			points3d(nb_couples).Init(y, y, x);
			expected_orbit = 12;
		      }
		    else if (abs(z-t) < 1e-13)
		      {
			points3d(nb_couples).Init(z, z, x);
			expected_orbit = 12;
		      }
		    else
		      points3d(nb_couples).Init(x, y, z);
		    
		    weights3d(nb_couples) = Weights3D(j);
		    break;
		  }
	      
	      for (int j = 0; j < nb_points; j++)
		for (int m = 0; m < 24; m++)
		  if (!PointUsed(j))
		    if (Points3D(j) == sym_points(m))
		      {
			PointUsed(j) = true;
			nb_remaining_points--;
			num_points(k) = j;
			other_points(k) = Points3D(j);
			//DISP(j); DISP(sym_points(m));
			k++;
		      }
	      
	      NbPointsCouple(nb_couples) = k;
	      
	      if (expected_orbit != k)
		{
		  abort();
		}
	      
	      if (k == 1)
		n1++;
	      else if (k == 4)
		n4++;
	      else if (k == 6)
		n6++;
	      else if (k == 12)
		n12++;
	      else if (k == 24)
		n24++;
	      
	      nb_couples++;	  	  
	    }
	  
	  DISP(n1);
	  bool cpp_code = false;
	  if (cpp_code)
	    {
	      file_out << "case " << p << " : " << endl;
	      file_out <<" { \n n1 = " << n1 << "; n4 = "<<n4<<"; n6 = "<<n6<<"; n12 = "<<n12<<"; n24 = "<<n24<<";\n";
	      file_out << "nb_points = n1 + 4*n4 + 6*n6 + 12*n12 + 24*n24;" << endl;
	      file_out << "points3d.Reallocate(nb_points); \n weights3d.Reallocate(nb_points); \n";
	      int ind = 0;
	      for (int j = 0; j < nb_couples; j++)
		if (NbPointsCouple(j) == 1)
		  {
		    file_out << "points3d(0).Init(0.25, 0.25, 0.25);\n weights3d(0) = " << weights3d(j) << ";\n";
		    ind++;
		  }
	      
	      for (int j = 0; j < nb_couples; j++)
		if (NbPointsCouple(j) == 4)
		  {
		    file_out << "points3d("<<ind<<").Init("<< points3d(j)(0) <<", "<<points3d(j)(1)<<", "<<points3d(j)(2)<<");\n";
		    file_out << "weights3d(" << ind << ") = " << weights3d(j) << ";\n";
		    ind += 4;
		  }
	  	  
	      for (int j = 0; j < nb_couples; j++)
		if (NbPointsCouple(j) == 6)
		  {
		    file_out << "points3d("<<ind<<").Init("<< points3d(j)(0) <<", "<<points3d(j)(1)<<", "<<points3d(j)(2)<<");\n";
		    file_out << "weights3d(" << ind << ") = " << weights3d(j) << ";\n";
		    ind += 6;
		  }

	      for (int j = 0; j < nb_couples; j++)
		if (NbPointsCouple(j) == 12)
		  {
		    file_out << "points3d("<<ind<<").Init("<< points3d(j)(0) <<", "<<points3d(j)(1)<<", "<<points3d(j)(2)<<");\n";
		    file_out << "weights3d(" << ind << ") = " << weights3d(j) << ";\n";
		    ind += 12;
		  }

	      for (int j = 0; j < nb_couples; j++)
		if (NbPointsCouple(j) == 24)
		  {
		    file_out << "points3d("<<ind<<").Init("<< points3d(j)(0) <<", "<<points3d(j)(1)<<", "<<points3d(j)(2)<<");\n";
		    file_out << "weights3d(" << ind << ") = " << weights3d(j) << ";\n";
		    ind += 24;
		  }
	      
	      file_out << "} \n break; \n ";
	    }
	  else
	    {
	      file_out << "case(" << p << ") " << endl;
	      file_out <<"\n n1 = " << n1 << "; n4 = "<<n4<<"; n6 = "<<n6<<"; n12 = "<<n12<<"; n24 = "<<n24<<";\n";
	      file_out << "nb_points = n1 + 4*n4 + 6*n6 + 12*n12 + 24*n24" << endl;
	      file_out << "call resize(points3d, nb_points) \n call resize(weights3d, nb_points) \n";
	      int ind = 1;
	      for (int j = 0; j < nb_couples; j++)
		if (NbPointsCouple(j) == 1)
		  {
		    file_out << "points3d(:, 1) = (/ 0.25_wp, 0.25_wp, 0.25_wp /) \n weights3d(1) = " << weights3d(j) << "_wp\n";
		    ind++;
		  }
	      
	      for (int j = 0; j < nb_couples; j++)
		if (NbPointsCouple(j) == 4)
		  {
		    file_out << "points3d(:, "<<ind<<") = (/"<< points3d(j)(0) <<"_wp, "<<points3d(j)(1)<<"_wp, "<<points3d(j)(2)<<"_wp /)\n";
		    file_out << "weights3d(" << ind << ") = " << weights3d(j) << "_wp;\n";
		    ind += 4;
		  }
	      
	      for (int j = 0; j < nb_couples; j++)
		if (NbPointsCouple(j) == 6)
		  {
		    file_out << "points3d(:, "<<ind<<") = (/"<< points3d(j)(0) <<"_wp, "<<points3d(j)(1)<<"_wp, "<<points3d(j)(2)<<"_wp /)\n";
		    file_out << "weights3d(" << ind << ") = " << weights3d(j) << "_wp;\n";
		    ind += 6;
		  }

	      for (int j = 0; j < nb_couples; j++)
		if (NbPointsCouple(j) == 12)
		  {
		    file_out << "points3d(:, "<<ind<<") = (/"<< points3d(j)(0) <<"_wp, "<<points3d(j)(1)<<"_wp, "<<points3d(j)(2)<<"_wp /)\n";
		    file_out << "weights3d(" << ind << ") = " << weights3d(j) << "_wp;\n";
		    ind += 12;
		  }

	      for (int j = 0; j < nb_couples; j++)
		if (NbPointsCouple(j) == 24)
		  {
		    file_out << "points3d(:, "<<ind<<") = (/"<< points3d(j)(0) <<"_wp, "<<points3d(j)(1)<<"_wp, "<<points3d(j)(2)<<"_wp /)\n";
		    file_out << "weights3d(" << ind << ") = " << weights3d(j) << "_wp;\n";
		    ind += 24;
		  }
	      
	      file_out << "\n";
	    }
	}
    }
  
  file_out.close();
  
}


void CheckPyramidQuadrature(int p, const VectR3& Points3D, const VectReal_wp& Weights3D)
{
  // integration de C_p
  bool test_integral = true;
  for (int k = 0; k <= p; k++)
    for (int i = 0; i <= k; i++)
      for (int j = 0; j <= k-i; j++)
        {
          Real_wp val_ex = 4.0/Real_wp(k+3);
          if (i%2 == 1)
            val_ex = 0.0;
          else
            val_ex *= 1.0/Real_wp(i+1);
          
          if (j%2 == 1)
            val_ex = 0.0;
          else
            val_ex *= 1.0/Real_wp(j+1);
          
          Real_wp val_num = 0;
          for (int n = 0; n < Points3D.GetM(); n++)
            val_num += Weights3D(n)*pow(Points3D(n)(0), i)*pow(Points3D(n)(1), j)*pow(1.0-Points3D(n)(2), k-i-j);
          
          if (abs(val_num - val_ex) > 1e-12)
            {
              cout << "Integrale de x^" << i << " y^" << j << " (1-z)^" << k << " differe" << endl;
              DISP(val_ex); DISP(val_num);
              test_integral = false;
            }
        }
  
  if (test_integral)
    cout << "Rule exact for C_" << p << endl;
}

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  //VectReal_wp xi, omega;
  //int r = to_num<int>(argv[1]);
  //ComputeGaussJacobi(xi, omega, r, 2.0, 0.0);
  //DISP(xi); DISP(omega);

  // PyramidClassical tri;
  //TetrahedronClassical tri;
  //TriangleClassical tri;
  // tri.ConstructFiniteElement(r);
  
  
  R2::threshold = 1e-12;
  R3::threshold = 1e-12;
  // GenerateTriangleFormulas();
  //VectR2 Points2D; VectReal_wp Weights2D;
  VectR3 Points3D; VectReal_wp Weights3D;
  VectReal_wp points, weights;
  int p = atoi(argv[1]);
  //for (int p = 1; p <= 10; p++)
    {
      PyramidQuadrature::ConstructPolynomialRule(p, Points3D, Weights3D);
      //PyramidQuadrature::ConstructQuadrature(p/2, Points3D, Weights3D, points, weights, PyramidQuadrature::QUADRATURE_JACOBI2);
      CheckPyramidQuadrature(p, Points3D, Weights3D);
      
      //TriangleQuadrature::ConstructQuadrature(p, Points2D, Weights2D, TriangleQuadrature::QUADRATURE_GAUSS);
      //CheckTriangleQuadrature(p, Points2D, Weights2D);

      //TetrahedronQuadrature::ConstructQuadrature(p, Points3D, Weights3D, TetrahedronQuadrature::QUADRATURE_QUASI_LUMPED);
      //CheckTetrahedronQuadrature(p, Points3D, Weights3D);
    }
  
  // GenerateTetrahedronFormulas();
  
  
  return FinalizeMontjoie();
}
