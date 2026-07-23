#ifndef MONTJOIE_FILE_PYRAMID_QUADRATURE_CXX

namespace Montjoie
{
  
  //! constructs quadrature formulas over symmetric pyramid
  /*!
    \param[in] r order of quadrature rule
    \param[out] points3d 3-D integration points \xi_i
    \param[out] weights3d 3-D integration weights \omega_i
    \param[out] points1d_z 1-D points along z axis
    \param[out] weights1d_z 1-D weights along z axis
    \param[in] type_quad type of quadrature
    The quadrature formulas can be used to integrate :
    \int_K f(x, y, z) = \sum \omega_i f(\xi_i)
    where K is the symmetric pyramid defined by :
    -(1-z) <= x, y <= (1-z)   and  0 <= z <= 1
   */
  void PyramidQuadrature
  ::ConstructQuadrature(int r, VectR3& points3d, VectReal_wp& weights3d,
                        VectReal_wp& points1d_z, VectReal_wp& weights1d_z, int type_quad)
  {
    VectReal_wp points1d, weights1d;
    ComputeGaussLegendre(points1d, weights1d, r);
        
    // quadrature points on the volume
    if (type_quad == QUADRATURE_GAUSS)
      {
        points3d.Reallocate((r+1)*(r+1)*(r+1));
        weights3d.Reallocate(points3d.GetM());
        int node = 0;
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r; j++)
            for (int k = 0; k <= r; k++)
              {
                Real_wp x = points1d(i), y = points1d(j), z = points1d(k);
                points3d(node).Init((2.0*x-1.0)*(1.0-z), (2.0*y-1.0)*(1.0-z), z);
                weights3d(node) = 4.0*(1.0-z)*(1.0-z)*weights1d(i)*weights1d(j)*weights1d(k);
                node++;
              }
      }
    else if (type_quad == QUADRATURE_JACOBI1)
      {
        ComputeGaussJacobi(points1d_z, weights1d_z, r, Real_wp(1), Real_wp(0));
        points3d.Reallocate((r+1)*(r+1)*(r+1));
        weights3d.Reallocate(points3d.GetM());
        int node = 0;
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r; j++)
            for (int k = 0; k <= r; k++)
              {
                Real_wp x = points1d(i), y = points1d(j), z = points1d_z(k);
                points3d(node).Init((2.0*x-1.0)*(1.0-z), (2.0*y-1.0)*(1.0-z), z);
                weights3d(node) = 4.0*(1.0-z)*weights1d(i)*weights1d(j)*weights1d_z(k);
                node++;
              }
      }
    else if (type_quad == QUADRATURE_JACOBI2)
      {
        ComputeGaussJacobi(points1d_z, weights1d_z, r, Real_wp(2), Real_wp(0));
        points3d.Reallocate((r+1)*(r+1)*(r+1));
        weights3d.Reallocate(points3d.GetM());
        int node = 0;
        for (int i = 0; i <= r; i++)
          for (int j = 0; j <= r; j++)
            for (int k = 0; k <= r; k++)
              {
                Real_wp x = points1d(i), y = points1d(j), z = points1d_z(k);
                points3d(node).Init((2.0*x-1.0)*(1.0-z), (2.0*y-1.0)*(1.0-z), z);
                weights3d(node) = 4.0*weights1d(i)*weights1d(j)*weights1d_z(k);
                node++;
              }
      }
    else if (type_quad == QUADRATURE_NON_PRODUCT)
      {
        int p = 2*r+1;
        switch(r)
          {
          case 1 : p = 2; break;
          case 2 : p = 5; break;
          case 3 : p = 7; break;
          case 4 : p = 9; break;
          }
        
        ConstructPolynomialRule(p, points3d, weights3d);
      }
    else
      {
        abort();
      }
  }
  
  
  void PyramidQuadrature::ConstructPolynomialRule(int p, VectR3& points3d, VectReal_wp& weights3d)
  {
    int n4 = 0, n8 = 0;
    switch (p)
      {
      case 0 :
      case 1 :
        //n1 = 1;
        weights3d.Reallocate(1);
        weights3d(0) = 8.0/3.0;
        points3d.Reallocate(1);
        points3d(0).Init(0.0, 0.0, -1.0/2.0);
        break;
      case 2 :
        /*weights3d.Reallocate(4);
        points3d.Reallocate(4);
        weights3d(0) = 0.699453551912568;
        points3d(0).Init(-0.750000000000000, 0.0, -0.5);
        weights3d(1) = 0.699453551912568;
        points3d(1).Init(-0.266666666666667, 0.700991361493770, -0.5);
        weights3d(2) = 0.560442489670798;
        points3d(2).Init(-0.266666666666667, -0.386753854617252, 0.131034482758621);
        weights3d(3) = 0.707317073170732;
        points3d(3).Init(-0.266666666666667, -0.386753854617252, -1.0); */
        
        /*n4 = 1; //n1 = 1;
        weights3d.Reallocate(5);
        weights3d.Fill(1.6/3);
        points3d.Reallocate(5);
        points3d(0).Init(0.5, 0.5, -0.887298334620742);
        points3d(4).Init(0.0, 0.0, 0.274596669241483); */
        
        {
          n4 = 1;
          Real_wp gamma1 = -0.5 + sqrt(Real_wp(3.0)/Real_wp(5.0));
          Real_wp gamma2 = -0.5 - 0.25*sqrt(Real_wp(3.0)/Real_wp(5.0));
          weights3d.Reallocate(5); points3d.Reallocate(5);
          weights3d(0) = 8.0/15.0;
          points3d(0).Init(0.5, 0.5, gamma2);
          weights3d(4) = 8.0/15.0;
          points3d(4).Init(0.0, 0.0, gamma1);
        }
        break;
      case 3 :
        n4 = 1; //n1 = 2; 
        weights3d.Reallocate(6);
        weights3d(0) = Real_wp(36)/85;
        weights3d(4) = Real_wp(3)/10;
        weights3d(5) = Real_wp(343)/510;
        points3d.Reallocate(6);
        points3d(0).Init(0.793492047615872, 0.0, -2.0/3);
        points3d(4).Init(0.0, 0.0, -1.0);
        points3d(5).Init(0.0, 0.0, 0.142857142857143);
        break;
      case 4 :
        n4 = 2; //n1 = 2;
        weights3d.Reallocate(10);
        points3d.Reallocate(10);
        weights3d(0) = 0.169384241783417;
        points3d(0).Init(0.657966997126504, 0.657966997126504, -0.921503432203399);
        weights3d(4) = 0.283532234381299;
        points3d(4).Init(0.650581556392544, 0.0, -0.355231700855205);
        weights3d(8) = 0.551689073546826;
        points3d(8).Init(0.0, 0.0, -0.749726093801093);
        weights3d(9) = 0.303311688460986;
        points3d(9).Init(0.0, 0.0, 0.354465577768317);
        break;
      case 5 :
        /*n4 = 3; //n1 = 3;
        weights3d.Reallocate(15);
        points3d.Reallocate(15);
        weights3d(0) = 0.097062415074088;
        points3d(0).Init(0.719461795119535, 0.719461795119535, -0.882163616032127);
        weights3d(4) = 0.153089998418834;
        points3d(4).Init(0.440029973561290, 0.440029973561290, -0.151149945443633);
        weights3d(8) = 0.224973598340578;
        points3d(8).Init(0.689317760758028, 0.0, -0.750000000000000);
        weights3d(12) = 0.140916095906755;
        points3d(12).Init(0.0, 0.0, -1.014030732756705);
        weights3d(13) = 0.443318442055681;
        points3d(13).Init(0.0, 0.0, -0.370504617476054);
        weights3d(14) = 0.181928081370229;
        points3d(14).Init(0.0, 0.0, 0.460567421625817);*/
        
        n4 = 4; // n1 = 1;
        weights3d.Reallocate(17);
        points3d.Reallocate(17);
        weights3d(0) = 0.072421368643694;
        points3d(0).Init(0.798041533155477, 0.798041533155477, -0.809167735391489);
        weights3d(4) = 0.061389172107864;
        points3d(4).Init(0.443342939439341, 0.443342939439341, -0.017170927262009);
        weights3d(8) = 0.230560253161839;
        points3d(8).Init(0.639204459383276, 0.0, -0.873027592022377);
        weights3d(12) = 0.255653204690824;
        points3d(12).Init(0.468380753435554, 0.0, -0.365143572456810);
        weights3d(16) = 0.186570672249784;
        points3d(16).Init(0.0, 0.0, 0.449318156750961);
        break;
      case 6 :
        weights3d.Reallocate(23);
        points3d.Reallocate(23);
        weights3d(0) = 0.061960020365082;
        weights3d(1) = 0.061960020365082;
        weights3d(2) = 0.061960020365082;
        weights3d(3) = 0.061960020365082;
        points3d(0).Init(-0.827892022220389, 0.228553292960554, -0.975368124041314);
        points3d(1).Init(0.827892022220389, -0.228553292960554, -0.975368124041314);
        points3d(2).Init(-0.228553292960554, 0.827892022220389, -0.975368124041314);
        points3d(3).Init(0.228553292960554, -0.827892022220389, -0.975368124041314);
        weights3d(4) = 0.114181837798273;
        weights3d(5) = 0.114181837798273;
        weights3d(6) = 0.114181837798273;
        weights3d(7) = 0.114181837798273;
        points3d(4).Init(-0.733314084871403, -0.246957217097095, -0.534303519773502);
        points3d(5).Init(0.733314084871403, 0.246957217097095, -0.534303519773502);
        points3d(6).Init(-0.246957217097095, -0.733314084871403, -0.534303519773502);
        points3d(7).Init(0.246957217097095, 0.733314084871403, -0.534303519773502);
        weights3d(8) = 0.064683549806387;
        weights3d(9) = 0.064683549806387;
        points3d(8).Init(-0.759322249982095, -0.759322249982095, -0.901297748120975);
        points3d(9).Init(0.759322249982095, 0.759322249982095, -0.901297748120975);
        weights3d(10) = 0.216648112769220;
        weights3d(11) = 0.216648112769220;
        points3d(10).Init(-0.287601758293552, -0.287601758293552, -0.853725850581071);
        points3d(11).Init(0.287601758293552, 0.287601758293552, -0.853725850581071);
        weights3d(12) = 0.052280462560358;
        weights3d(13) = 0.052280462560358;
        points3d(12).Init(0.781129332284802, -0.781129332284802, -0.774925175941892);
        points3d(13).Init(-0.781129332284802, 0.781129332284802, -0.774925175941892);
        weights3d(14) = 0.258213931202755;
        weights3d(15) = 0.258213931202755;
        points3d(14).Init(0.396094511344884, -0.396094511344884, -0.586090047053785);
        points3d(15).Init(-0.396094511344884, 0.396094511344884, -0.586090047053785);
        weights3d(16) = 0.199865246633000;
        weights3d(17) = 0.199865246633000;
        points3d(16).Init(0.169033458451920, 0.169033458451920, -0.131363109889133);
        points3d(17).Init(-0.169033458451920, -0.169033458451920, -0.131363109889133);
        weights3d(18) = 0.090143958829821;
        weights3d(19) = 0.090143958829821;
        points3d(18).Init(0.377437091528202, -0.377437091528202, 0.040035206144175);
        points3d(19).Init(-0.377437091528202, 0.377437091528202, 0.040035206144175);
        weights3d(20) = 0.053641092748552;
        weights3d(21) = 0.053641092748552;
        points3d(20).Init(0.345937985742664, 0.345937985742664, 0.146410060502067);
        points3d(21).Init(-0.345937985742664, -0.345937985742664, 0.146410060502067);
        weights3d(22) = 0.091146524913060;
        points3d(22).Init(0.0, 0.0, 0.573059549806815);
        break;
      case 7 :
        /*n4 = 7; //n1 = 4;
        weights3d.Reallocate(32);
        points3d.Reallocate(32);
        weights3d(0) = 0.100194280213435;
        points3d(0).Init(0.367156142481928, 0.367156142481928, -0.920057020952451);
        weights3d(4) = 0.052887805212567;
        points3d(4).Init(0.766599367404359, 0.766599367404359, -0.904109268863743);
        weights3d(8) = 0.190435579642634;
        points3d(8).Init(0.264977352419382, 0.264977352419382, -0.546968302614360);
        weights3d(12) = 0.077167203781053;
        points3d(12).Init(0.603600802215458, 0.603600802215458, -0.506946742162290);
        weights3d(16) = 0.077949744563810;
        points3d(16).Init(0.279829398847256, 0.279829398847256, 0.179508326447283);
        weights3d(20) = 0.071178024781337;
        points3d(20).Init(0.864098759787715, 0.0, -0.866666666666667);
        weights3d(24) = 0.076562500000028;
        points3d(24).Init(0.617213399848321, 0.0, -0.333333333333333);
        weights3d(28) = -0.183091556415915;
        points3d(28).Init(0.0, 0.0, -0.422227897975400);
        weights3d(29) = 0.222463043407430;
        points3d(29).Init(0.0, 0.0, -0.065238394218177);
        weights3d(30) = 0.049690153330811;
        points3d(30).Init(0.0, 0.0, 0.749146150859227);
        weights3d(31) = -0.007895526435117;
        points3d(31).Init(0.0, 0.0, 0.931878142180172);*/
        
        n4 = 8; // n1 = 1;
        weights3d.Reallocate(33);
        points3d.Reallocate(33);
        weights3d(0) = 0.085676176018999;
        points3d(0).Init(0.368103812374190, 0.368103812374190, -0.935295641606863);
        weights3d(4) = 0.051269218633384;
        points3d(4).Init(0.768857166809385, 0.768857166809385, -0.908041262677320);
        weights3d(8) = 0.159827658208448;
        points3d(8).Init(0.306891446901222, 0.306891446901222, -0.611361877687655);
        weights3d(12) = 0.075035550101029;
        points3d(12).Init(0.609244741163901, 0.609244741163901, -0.511276520841257);
        weights3d(16) = 0.090156104778395;
        points3d(16).Init(0.193453966706627, 0.193453966706627, -0.025433079112450);
        weights3d(20) = 0.042264816233125;
        points3d(20).Init(0.306935338280579, 0.306935338280579, 0.228313791540606);
        weights3d(24) = 0.071178024781071;
        points3d(24).Init(0.864098759788004, 0.0, -0.866666666667237);
        weights3d(28) = 0.076562500000108;
        points3d(28).Init(0.617213399849506, 0.0, -0.333333333336100);
        weights3d(32) = 0.058786471648435;
        points3d(32).Init(0.0, 0.0, 0.631573924240567);
        break;
      case 8 :
        /*n8 = 1; n4 = 9; // n1 = 2;
        weights3d.Reallocate(46);
        points3d.Reallocate(46);
        weights3d(0) = 0.023838342903583;
        points3d(0).Init(0.889748462878860, 0.524423401771527, -0.958535174038917);
        weights3d(8) = 0.054117091870623;
        points3d(8).Init(0.300849716344273, 0.300849716344273, -0.938570100408951);
        weights3d(12) = 0.023835574737103;
        points3d(12).Init(0.809860521960780, 0.809860521960780, -0.782213283414361);
        weights3d(16) = 0.087532227329390;
        points3d(16).Init(0.476421500507226, 0.476421500507226, -0.705910687424969);
        weights3d(20) = 0.096769980749473;
        points3d(20).Init(0.481571723769423, 0.481571723769423, -0.426083866123888);
        weights3d(24) = 0.009918864214805;
        points3d(24).Init(0.552525436533552, 0.552525436533552, -0.038433293418787);
        weights3d(28) = 0.037785656462796;
        points3d(28).Init(0.168708304836178, 0.168708304836178, 0.428377493198300);
        weights3d(32) = 0.072301468548313;
        points3d(32).Init(0.698562689379681, 0.0, -0.870630768348530);
        weights3d(36) = 0.047029494133604;
        points3d(36).Init(0.791343124611894, 0.0, -0.584001111150429);
        weights3d(40) = 0.122694030065541;
        points3d(40).Init(0.374288082186939, 0.0, -0.056825101556676);
        weights3d(44) = 0.258379638668609;
        points3d(44).Init(0.0, 0.0, -0.550796514382408);
        weights3d(45) = 0.009642732322797;
        points3d(45).Init(0.0, 0.0, 0.823127024581939);*/
        
        n8 = 2; n4 = 7; // n1 = 4;
        weights3d.Reallocate(48);
        points3d.Reallocate(48);
        weights3d(0) = 0.025969596552691;
        points3d(0).Init(0.874496930876371, 0.467687186260233, -0.963816320572130);
        weights3d(8) = 0.034434101526449;
        points3d(8).Init(0.788619323446728, 0.106279251382570, -0.676881943842562);
        weights3d(16) = 0.021100483265270;
        points3d(16).Init(0.829065607340491, 0.829065607340491, -0.838053411841242);
        weights3d(20) = 0.113756971169142;
        points3d(20).Init(0.510947174328476, 0.510947174328476, -0.680630507039655);
        weights3d(24) = 0.015734010695881;
        points3d(24).Init(0.656721871362728, 0.656721871362728, -0.354301071531389);
        weights3d(28) = 0.096871999024165;
        points3d(28).Init(0.368483094260005, 0.368483094260005, -0.272073996132499);
        weights3d(32) = 0.037948587645290;
        points3d(32).Init(0.246229581685580, 0.246229581685580, 0.329810776523420);
        weights3d(36) = 0.098799327331805;
        points3d(36).Init(0.491778676738584, 0.0, -0.913630070654998);
        weights3d(40) = 0.059982554235940;
        points3d(40).Init(0.503396929834625, 0.0, -0.132398768099718);
        weights3d(44) = 0.246909857116978;
        points3d(44).Init(0.0, 0.0, -0.542336234274656);
        weights3d(45) = 0.120187722312389;
        points3d(45).Init(0.0, 0.0, 0.114808588761911);
        weights3d(46) = 0.011710066663757;
        points3d(46).Init(0.0, 0.0, 0.379206250032304);
        weights3d(47) = 0.027853702470455;
        points3d(47).Init(0.0, 0.0, 0.716879813872603);
        break;        
      case 9 :
        n8 = 7; n4 = 5; // n1 = 2;
        weights3d.Reallocate(78);
        points3d.Reallocate(78);
        weights3d(0) = 0.001339095382101;
        points3d(0).Init(0.533017207297908, 0.172596065690238, -1.138765456262421);
        weights3d(8) = 0.025194086642736;
        points3d(8).Init(0.886106901903220, 0.333803313872032, -0.923676148202636);
        weights3d(16) = 0.064930141290614;
        points3d(16).Init(0.478020191825880, 0.142511506898247, -0.875253195996511);
        weights3d(24) = 0.035391605594706;
        points3d(24).Init(0.749736308663014, 0.271722978493158, -0.612963089790895);
        weights3d(32) = 0.065039916059217;
        points3d(32).Init(0.382904862544845, 0.108978010806538, -0.502821761881133);
        weights3d(40) = 0.055612199255489;
        points3d(40).Init(0.245462634502461, 0.468333316506249, -0.120081804327846);
        weights3d(48) = 0.016672665315589;
        points3d(48).Init(0.254855721354154, 0.143006403348963, 0.395085664029986);
        weights3d(56) = 0.033703591459062;
        points3d(56).Init(0.766249052769552, 0.766249052769552, -0.925989117611179);
        weights3d(60) = 0.000169600348221;
        points3d(60).Init(1.254395263858391, 1.254395263858391, -0.870855710453712);
        weights3d(64) = 0.055530265982886;
        points3d(64).Init(0.640261671994421, 0.640261671994421, -0.624603536324118);
        weights3d(68) = 0.003410468611445;
        points3d(68).Init(0.664149290834120, 0.664149290834120, -0.237688208519910);
        weights3d(72) = 0.019357771948949;
        points3d(72).Init(0.054282651302894, 0.054282651302894, 0.194831715390382);
        weights3d(76) = 0.085548572634762;
        points3d(76).Init(0.0, 0.0, -0.074744356612811);
        weights3d(77) = 0.018993624306037;
        points3d(77).Init(0.0, 0.0, 0.752811544075902);
        break;
      default :
	{
	  cout << "formula not available" << endl;
	  abort();
	}
      }
    
    // generating other points by symmetry
    int nb = 0;
    for (int i = 0; i < n8; i++)
      {
        Real_wp x = points3d(nb)(0), y = points3d(nb)(1), z = points3d(nb)(2);
        weights3d(nb+1) = weights3d(nb);
        weights3d(nb+2) = weights3d(nb);
        weights3d(nb+3) = weights3d(nb);
        weights3d(nb+4) = weights3d(nb);
        weights3d(nb+5) = weights3d(nb);
        weights3d(nb+6) = weights3d(nb);
        weights3d(nb+7) = weights3d(nb);
        points3d(nb+1).Init(x, -y, z);
        points3d(nb+2).Init(-x, -y, z);
        points3d(nb+3).Init(-x, y, z);
        points3d(nb+4).Init(y, x, z);
        points3d(nb+5).Init(y, -x, z);
        points3d(nb+6).Init(-y, x, z);
        points3d(nb+7).Init(-y, -x, z);
        nb += 8;
      }
    
    for (int i = 0; i < n4; i++)
      {
        Real_wp x = points3d(nb)(0), y = points3d(nb)(1), z = points3d(nb)(2);
        weights3d(nb+1) = weights3d(nb);
        weights3d(nb+2) = weights3d(nb);
        weights3d(nb+3) = weights3d(nb);
        if (y == 0.0)
          {
            points3d(nb+1).Init(-x, 0.0, z);
            points3d(nb+2).Init(0.0, x, z);
            points3d(nb+3).Init(0.0, -x, z);
          }
        else
          {
            points3d(nb+1).Init(x, -x, z);
            points3d(nb+2).Init(-x, x, z);
            points3d(nb+3).Init(-x, -x, z);
          }
        
        nb += 4;
      }
    
    //DISP(points3d); DISP(weights3d);
    // computation of points on the symmetric pyramid
    Real_wp poids = 1.0/2;
    for (int i = 0; i < weights3d.GetM(); i++)
      {
        weights3d(i) *= poids;
        points3d(i)(2) = 0.5*points3d(i)(2) + 0.5;
      }
  }
  
}

#define MONTJOIE_PYRAMID_QUADRATURE_CXX
#endif
