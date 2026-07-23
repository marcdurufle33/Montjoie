#ifndef MONTJOIE_FILE_TRIANGLE_QUADRATURE_CXX

namespace Montjoie
{
  
  //! constructs quadrature formulas over unit triangle
  /*!
    \param[in] p polynomials of degree at most p must be integrated exactly
    \param[out] points2d 2-D integration points \xi_i
    \param[out] weights2d 2-D integration weights \omega_i
    \param[in] type_quad type of quadrature
    The quadrature formulas can be used to integrate :
    \int_K f(x, y, z) = \sum \omega_i f(\xi_i)
    where K is the unit triangle defined by :
    x, y >= 0  and x+y <= 1
   */
  void TriangleQuadrature
  ::ConstructQuadrature(int p, VectR2& points2d, VectReal_wp& weights2d, int type_quad)
  {
    // quadrature formula over a triangle
    // Gauss-like quadrature -> points are strictely inside the triangle
    int n1, n3, n6, nb_points = 0;
    bool formula_not_available = false;
    
    if (type_quad == QUADRATURE_GAUSS)
      {
        switch (p)
          {
          case 0:
          case 1:
            { 
              n1 = 1; n3 = 0; n6 = 0;
              nb_points = n1 + 3*n3 + 6*n6;
              points2d.Reallocate(nb_points); 
              weights2d.Reallocate(nb_points); 
              points2d(0).Init(1.0/3, 1.0/3);
              weights2d(0) = 0.5;
            } 
            break; 
          case 2 : 
            { 
              n1 = 0; n3 = 1; n6 = 0;
              nb_points = n1 + 3*n3 + 6*n6;
              points2d.Reallocate(nb_points); 
              weights2d.Reallocate(nb_points); 
              points2d(0).Init(0.1666666666666665, 0.666666666666667); 
              weights2d(0) = 0.1666666666666667;
            } 
            break;
          case 3 :
            { 
              n1 = 1; n3 = 1; n6 = 0;
              nb_points = n1 + 3*n3 + 6*n6;
              points2d.Reallocate(nb_points); 
              weights2d.Reallocate(nb_points); 
              points2d(0).Init(1.0/3, 1.0/3);
              weights2d(0) = -0.28125;
              points2d(1).Init(0.2, 0.6); 
              weights2d(1) = 0.2604166666666667;
            } 
            break; 
          case 4 : 
            { 
              n1 = 0; n3 = 2; n6 = 0;
              nb_points = n1 + 3*n3 + 6*n6;
              points2d.Reallocate(nb_points); 
              weights2d.Reallocate(nb_points);
              //points2d(0).Init(0, 0);
              points2d(0).Init(0.445948490915965, 0.10810301816807); 
              weights2d(0) = 0.1116907948390055;
              //points2d(3).Init(0.5, 0.0);
              points2d(3).Init(0.09157621350977102, 0.816847572980458); 
              weights2d(3) = 0.054975871827661;
            } 
            break; 
          case 5 : 
            { 
              n1 = 1; n3 = 2; n6 = 0;
              nb_points = n1 + 3*n3 + 6*n6;
              points2d.Reallocate(nb_points); 
              weights2d.Reallocate(nb_points); 
              points2d(0).Init(1.0/3, 1.0/3);
              weights2d(0) = 0.1125;
              points2d(1).Init(0.470142064105115, 0.05971587178977); 
              weights2d(1) = 0.066197076394253;
              points2d(4).Init(0.101286507323456, 0.797426985353088); 
              weights2d(4) = 0.0629695902724135;
            } 
            break; 
          case 6 : 
            { 
              n1 = 0; n3 = 2; n6 = 1;
              nb_points = n1 + 3*n3 + 6*n6;
              points2d.Reallocate(nb_points); 
              weights2d.Reallocate(nb_points); 
              points2d(0).Init(0.24928674517091, 0.50142650965818); 
              weights2d(0) = 0.0583931378631895;
              points2d(3).Init(0.063089014491502, 0.873821971016996); 
              weights2d(3) = 0.0254224531851035;
              points2d(6).Init(0.636502499121399, 0.310352451033784); 
              weights2d(6) = 0.041425537809187;
            } 
            break; 
          case 7 : 
            { 
              n1 = 1; n3 = 2; n6 = 1;
              nb_points = n1 + 3*n3 + 6*n6;
              points2d.Reallocate(nb_points); 
              weights2d.Reallocate(nb_points); 
              points2d(0).Init(1.0/3, 1.0/3);
              weights2d(0) = -0.07478502223384099;
              points2d(1).Init(0.26034596607904, 0.4793080678419199); 
              weights2d(1) = 0.087807628716604;
              points2d(4).Init(0.06513010290221599, 0.869739794195568); 
              weights2d(4) = 0.026673617804419;
              points2d(7).Init(0.63844418856981, 0.312865496004874); 
              weights2d(7) = 0.0385568804451285;
            } 
            break; 
          case 8 : 
            { 
              n1 = 1; n3 = 3; n6 = 1;
              nb_points = n1 + 3*n3 + 6*n6;
              points2d.Reallocate(nb_points); 
              weights2d.Reallocate(nb_points); 
              points2d(0).Init(1.0/3, 1.0/3);
              weights2d(0) = 0.0721578038388935;
              points2d(1).Init(0.459292588292723, 0.08141482341455397); 
              weights2d(1) = 0.0475458171336425;
              points2d(4).Init(0.17056930775176, 0.65886138449648); 
              weights2d(4) = 0.051608685267359;
              points2d(7).Init(0.05054722831703101, 0.898905543365938); 
              weights2d(7) = 0.016229248811599;
              points2d(10).Init(0.728492392955404, 0.263112829634638); 
              weights2d(10) = 0.0136151570872175;
            } 
            break; 
          case 9 : 
            { 
              n1 = 1; n3 = 4; n6 = 1;
              nb_points = n1 + 3*n3 + 6*n6;
              points2d.Reallocate(nb_points); 
              weights2d.Reallocate(nb_points); 
              points2d(0).Init(1.0/3, 1.0/3);
              weights2d(0) = 0.0485678981413995;
              points2d(1).Init(0.489682519198738, 0.02063496160252398); 
              weights2d(1) = 0.0156673501135695;
              points2d(4).Init(0.437089591492937, 0.125820817014126); 
              weights2d(4) = 0.038913770502387;
              points2d(7).Init(0.188203535619033, 0.623592928761934); 
              weights2d(7) = 0.039823869463605;
              points2d(10).Init(0.04472951339445302, 0.910540973211094); 
              weights2d(10) = 0.012788837829349;
              points2d(13).Init(0.741198598784498, 0.221962989160766); 
              weights2d(13) = 0.0216417696886445;
            } 
            break; 
          case 10 : 
            { 
              n1 = 1; n3 = 2; n6 = 3;
              nb_points = n1 + 3*n3 + 6*n6;
              points2d.Reallocate(nb_points); 
              weights2d.Reallocate(nb_points); 
              points2d(0).Init(1.0/3, 1.0/3);
              weights2d(0) = 0.045408995191377;
              points2d(1).Init(0.485577633383657, 0.02884473323268599); 
              weights2d(1) = 0.0183629788782335;
              points2d(4).Init(0.109481575485037, 0.781036849029926); 
              weights2d(4) = 0.022660529717764;
              points2d(7).Init(0.550352941820999, 0.307939838764121); 
              weights2d(7) = 0.03637895842271;
              points2d(13).Init(0.728323904597411, 0.246672560639903); 
              weights2d(13) = 0.0141636212655285;
              points2d(19).Init(0.9236559335875, 0.06680325101220003); 
              weights2d(19) = 0.0047108334818665;
            } 
            break; 
          case 11 : 
            { 
              n1 = 0; n3 = 5; n6 = 2;
              nb_points = n1 + 3*n3 + 6*n6;
              points2d.Reallocate(nb_points); 
              weights2d.Reallocate(nb_points); 
              points2d(0).Init(0.534611048270758, -0.06922209654151601); 
              weights2d(0) = 0.0004635031644805;
              points2d(3).Init(0.398969302965855, 0.20206139406829); 
              weights2d(3) = 0.0385747674574065;
              points2d(6).Init(0.203309900431282, 0.593380199137436); 
              weights2d(6) = 0.029661488690387;
              points2d(9).Init(0.119350912282581, 0.761298175434838); 
              weights2d(9) = 0.018092270251709;
              points2d(12).Init(0.03236494811127599, 0.935270103777448); 
              weights2d(12) = 0.006829865501339;
              points2d(15).Init(0.593201213428213, 0.356620648261293); 
              weights2d(15) = 0.026168555981102;
              points2d(21).Init(0.807489003159792, 0.171488980304042); 
              weights2d(21) = 0.0103538298195705;
            } 
            break; 
          case 12 : 
            { 
              n1 = 0; n3 = 5; n6 = 3;
              nb_points = n1 + 3*n3 + 6*n6;
              points2d.Reallocate(nb_points); 
              weights2d.Reallocate(nb_points); 
              points2d(0).Init(0.488217389773805, 0.02356522045238996); 
              weights2d(0) = 0.0128655332202275;
              points2d(3).Init(0.43972439229446, 0.12055121541108); 
              weights2d(3) = 0.021846272269019;
              points2d(6).Init(0.271210385012116, 0.4575792299757679); 
              weights2d(6) = 0.0314291121089425;
              points2d(9).Init(0.127576145541586, 0.7448477089168281); 
              weights2d(9) = 0.0173980564653545;
              points2d(12).Init(0.02131735045320998, 0.95736529909358); 
              weights2d(12) = 0.0030831305257795;
              points2d(15).Init(0.608943235779788, 0.275713269685514); 
              weights2d(15) = 0.0201857788831905;
              points2d(21).Init(0.6958360867878031, 0.28132558098994); 
              weights2d(21) = 0.0111783866011515;
              points2d(27).Init(0.858014033544073, 0.116251915907597); 
              weights2d(27) = 0.0086581155543295;
            } 
            break; 
          case 13 : 
            { 
              n1 = 1; n3 = 6; n6 = 3;
              nb_points = n1 + 3*n3 + 6*n6;
              points2d.Reallocate(nb_points); 
              weights2d.Reallocate(nb_points); 
              points2d(0).Init(1.0/3, 1.0/3);
              weights2d(0) = 0.026260461700401;
              points2d(1).Init(0.495048184939705, 0.009903630120590035); 
              weights2d(1) = 0.005640072604665;
              points2d(4).Init(0.468716635109574, 0.06256672978085198); 
              weights2d(4) = 0.015711759181227;
              points2d(7).Init(0.414521336801277, 0.170957326397446); 
              weights2d(7) = 0.023536251252097;
              points2d(10).Init(0.229399572042831, 0.541200855914338); 
              weights2d(10) = 0.0236817932681775;
              points2d(13).Init(0.11442449519633, 0.77115100960734); 
              weights2d(13) = 0.015583764522897;
              points2d(16).Init(0.024811391363459, 0.950377217273082); 
              weights2d(16) = 0.003987885732537;
              points2d(19).Init(0.63635117456166, 0.268794997058761); 
              weights2d(19) = 0.018424201364366;
              points2d(25).Init(0.690169159986905, 0.291730066734288); 
              weights2d(25) = 0.008700731651911001;
              points2d(31).Init(0.851409537834241, 0.126357385491669); 
              weights2d(31) = 0.0077608934195225;
            } 
            break; 
          case 14 : 
            { 
              n1 = 0; n3 = 6; n6 = 4;
              nb_points = n1 + 3*n3 + 6*n6;
              points2d.Reallocate(nb_points); 
              weights2d.Reallocate(nb_points); 
              points2d(0).Init(0.488963910362179, 0.02207217927564198); 
              weights2d(0) = 0.0109417906847145;
              points2d(3).Init(0.417644719340454, 0.164710561319092); 
              weights2d(3) = 0.0163941767720625;
              points2d(6).Init(0.273477528308839, 0.453044943382322); 
              weights2d(6) = 0.025887052253646;
              points2d(9).Init(0.177205532412543, 0.645588935174914); 
              weights2d(9) = 0.0210812943684965;
              points2d(12).Init(0.06179988309087298, 0.876400233818254); 
              weights2d(12) = 0.0072168498348885;
              points2d(15).Init(0.01939096124870099, 0.961218077502598); 
              weights2d(15) = 0.0024617018012;
              points2d(18).Init(0.770608554774996, 0.172266687821356); 
              weights2d(18) = 0.012332876606282;
              points2d(24).Init(0.570222290846683, 0.336861459796345); 
              weights2d(24) = 0.0192857553935305;
              points2d(30).Init(0.686980167808088, 0.298372882136258); 
              weights2d(30) = 0.007218154056767;
              points2d(36).Init(0.8797571713701711, 0.118974497696957); 
              weights2d(36) = 0.0025051144192505;
            } 
            break; 
          case 15 : 
            { 
              n1 = 0; n3 = 6; n6 = 5;
              nb_points = n1 + 3*n3 + 6*n6;
              points2d.Reallocate(nb_points); 
              weights2d.Reallocate(nb_points); 
              points2d(0).Init(0.506972916858243, -0.01394583371648594); 
              weights2d(0) = 0.0009584378214245;
              points2d(3).Init(0.431406354283023, 0.137187291433954); 
              weights2d(3) = 0.0221245136355725;
              points2d(6).Init(0.277693644847144, 0.444612710305712); 
              weights2d(6) = 0.025593274359426;
              points2d(9).Init(0.126464891041254, 0.747070217917492); 
              weights2d(9) = 0.011843867935344;
              points2d(12).Init(0.07080838597468597, 0.8583832280506281); 
              weights2d(12) = 0.0066448878450105;
              points2d(15).Init(0.01896517024107303, 0.9620696595178539); 
              weights2d(15) = 0.002374458304096;
              points2d(18).Init(0.604954466893291, 0.261311371140087); 
              weights2d(18) = 0.0192750362997965;
              points2d(24).Init(0.575586555512814, 0.388046767090269); 
              weights2d(24) = 0.013607907160312;
              points2d(30).Init(0.724462663076655, 0.285712220049916); 
              weights2d(30) = 0.0010910386833985;
              points2d(36).Init(0.747556466051838, 0.215599664072284); 
              weights2d(36) = 0.0107526599238655;
              points2d(42).Init(0.883964574092416, 0.103575616576386); 
              weights2d(42) = 0.0038369713155245;
            } 
            break; 
          case 16 : 
            { 
              n1 = 1; n3 = 7; n6 = 5;
              nb_points = n1 + 3*n3 + 6*n6;
              points2d.Reallocate(nb_points); 
              weights2d.Reallocate(nb_points); 
              points2d(0).Init(1.0/3, 1.0/3);
              weights2d(0) = 0.023437848713821;
              points2d(1).Init(0.497380541948438, 0.005238916103123992); 
              weights2d(1) = 0.0032029392892925;
              points2d(4).Init(0.413469438549352, 0.1730611229012959); 
              weights2d(4) = 0.0208551483696935;
              points2d(7).Init(0.470458599066991, 0.059082801866018); 
              weights2d(7) = 0.013445742125032;
              points2d(10).Init(0.240553749969521, 0.518892500060958); 
              weights2d(10) = 0.021066261380825;
              points2d(13).Init(0.147965794222573, 0.704068411554854); 
              weights2d(13) = 0.0150001334213865;
              points2d(16).Init(0.07546518765747401, 0.849069624685052); 
              weights2d(16) = 0.007100049462512;
              points2d(19).Init(0.01659640262302498, 0.96680719475395); 
              weights2d(19) = 0.0017912311756365;
              points2d(22).Init(0.599868711174861, 0.296555596579887); 
              weights2d(22) = 0.0163865737303135;
              points2d(28).Init(0.6421935249415049, 0.337723063403079); 
              weights2d(28) = 0.0076491531242205;
              points2d(34).Init(0.7995927209713269, 0.204748281642812); 
              weights2d(34) = 0.0011931220964195;
              points2d(40).Init(0.768699721401368, 0.189358492130623); 
              weights2d(40) = 0.0095423963779495;
              points2d(46).Init(0.900399064086661, 0.08528361568265702); 
              weights2d(46) = 0.003425027273271;
            } 
            break; 
          case 17 : 
            { 
              n1 = 1; n3 = 8; n6 = 6;
              nb_points = n1 + 3*n3 + 6*n6;
              points2d.Reallocate(nb_points); 
              weights2d.Reallocate(nb_points); 
              points2d(0).Init(1.0/3, 1.0/3);
              weights2d(0) = 0.0167185996454015;
              points2d(1).Init(0.497170540556774, 0.00565891888645198); 
              weights2d(1) = 0.0025467077202535;
              points2d(4).Init(0.482176322624625, 0.03564735475074998); 
              weights2d(4) = 0.007335432263819;
              points2d(7).Init(0.450239969020782, 0.09952006195843599); 
              weights2d(7) = 0.012175439176836;
              points2d(10).Init(0.400266239377397, 0.1994675212452059); 
              weights2d(10) = 0.0155537754344845;
              points2d(13).Init(0.252141267970953, 0.495717464058094); 
              weights2d(13) = 0.01562855560931;
              points2d(16).Init(0.162047004658461, 0.675905990683078); 
              weights2d(16) = 0.0124078271698325;
              points2d(19).Init(0.07587588226074599, 0.848248235478508); 
              weights2d(19) = 0.0070280365352785;
              points2d(22).Init(0.015654726967822, 0.968690546064356); 
              weights2d(22) = 0.0015973380868895;
              points2d(25).Init(0.655493203809423, 0.334319867363658); 
              weights2d(25) = 0.0040598276594965;
              points2d(31).Init(0.57233759053202, 0.292221537796944); 
              weights2d(31) = 0.0134028711415815;
              points2d(37).Init(0.626001190286228, 0.31957488542319); 
              weights2d(37) = 0.009229996605411;
              points2d(43).Init(0.796427214974071, 0.190704224192292); 
              weights2d(43) = 0.004238434267164;
              points2d(49).Init(0.752351005937729, 0.180483211648746); 
              weights2d(49) = 0.009146398385012499;
              points2d(55).Init(0.9046255040956079, 0.08071131367956402); 
              weights2d(55) = 0.0033328160020825;
            } 
            break; 
          case 18 : 
            { 
              n1 = 1; n3 = 9; n6 = 7;
              nb_points = n1 + 3*n3 + 6*n6;
              points2d.Reallocate(nb_points); 
              weights2d.Reallocate(nb_points); 
              points2d(0).Init(1.0/3, 1.0/3);
              weights2d(0) = 0.0154049699688235;
              points2d(1).Init(0.493344808630921, 0.01331038273815799); 
              weights2d(1) = 0.004536218339702;
              points2d(4).Init(0.469210594241957, 0.06157881151608602); 
              weights2d(4) = 0.009380658469797;
              points2d(7).Init(0.436281395887006, 0.127437208225988); 
              weights2d(7) = 0.0097205489927385;
              points2d(10).Init(0.394846170673416, 0.210307658653168); 
              weights2d(10) = 0.013876974305405;
              points2d(13).Init(0.249794568803157, 0.500410862393686); 
              weights2d(13) = 0.0161281126757285;
              points2d(16).Init(0.161432193743843, 0.677135612512314); 
              weights2d(16) = 0.012537016308461;
              points2d(19).Init(0.07659822748537098, 0.846803545029258); 
              weights2d(19) = 0.007635963985916;
              points2d(22).Init(0.02425243935345001, 0.9514951212931); 
              weights2d(22) = 0.0033969610114815;
              points2d(25).Init(0.04314636721696502, 0.91370726556607); 
              weights2d(25) = -0.00111154936496;
              points2d(28).Init(0.6326579688566361, 0.358911494940944); 
              weights2d(28) = 0.003165957038203;
              points2d(34).Init(0.574410971510855, 0.294402476751957); 
              weights2d(34) = 0.013628769024569;
              points2d(40).Init(0.6247790467925121, 0.325017801641814); 
              weights2d(40) = 0.0088383928247325;
              points2d(46).Init(0.748933176523037, 0.184737559666046); 
              weights2d(46) = 0.009189742319035001;
              points2d(52).Init(0.769207005420443, 0.218796800013321); 
              weights2d(52) = 0.004052366404096;
              points2d(58).Init(0.8839623022734671, 0.101179597136408); 
              weights2d(58) = 0.0038170645353625;
              points2d(64).Init(1.014347260005363, 0.02087475528258598); 
              weights2d(64) = 2.3093830397e-05;
            } 
            break; 
          case 19 : 
            { 
              n1 = 1; n3 = 8; n6 = 8;
              nb_points = n1 + 3*n3 + 6*n6;
              points2d.Reallocate(nb_points); 
              weights2d.Reallocate(nb_points); 
              points2d(0).Init(1.0/3, 1.0/3);
              weights2d(0) = 0.0164531656944595;
              points2d(1).Init(0.489609987073006, 0.02078002585398797); 
              weights2d(1) = 0.005165365945636;
              points2d(4).Init(0.454536892697893, 0.09092621460421402); 
              weights2d(4) = 0.011193623631508;
              points2d(7).Init(0.401416680649431, 0.1971666387011379); 
              weights2d(7) = 0.015133062934734;
              points2d(10).Init(0.255551654403098, 0.488896691193804); 
              weights2d(10) = 0.015245483901099;
              points2d(13).Init(0.17707794215213, 0.64584411569574); 
              weights2d(13) = 0.0120796063708205;
              points2d(16).Init(0.110061053227952, 0.779877893544096); 
              weights2d(16) = 0.0080254017934005;
              points2d(19).Init(0.05552862425184002, 0.88894275149632); 
              weights2d(19) = 0.004042290130892;
              points2d(22).Init(0.01262186377722901, 0.974756272445542); 
              weights2d(22) = 0.0010396810137425;
              points2d(25).Init(0.600633794794645, 0.395754787356943); 
              weights2d(25) = 0.0019424384524905;
              points2d(31).Init(0.557603261588784, 0.307929983880436); 
              weights2d(31) = 0.012787080306011;
              points2d(37).Init(0.720987025817365, 0.26456694840652); 
              weights2d(37) = 0.004440451786669;
              points2d(43).Init(0.594527068955871, 0.358539352205951); 
              weights2d(43) = 0.008062273380865501;
              points2d(49).Init(0.839331473680839, 0.157807405968595); 
              weights2d(49) = 0.0012459709087455;
              points2d(55).Init(0.701087978926173, 0.07505059697591099); 
              weights2d(55) = 0.009121420059475501;
              points2d(61).Init(0.822931324069857, 0.142421601113383); 
              weights2d(61) = 0.0051292818680995;
              points2d(67).Init(0.924344252620784, 0.06549462808293799); 
              weights2d(67) = 0.001899964427651;
            } 
            break; 
          case 20 : 
            { 
              n1 = 1; n3 = 8; n6 = 10;
              //n1 = 1; n3 = 2; n6 = 0;
              nb_points = n1 + 3*n3 + 6*n6;
              points2d.Reallocate(nb_points); 
              weights2d.Reallocate(nb_points);
              points2d.Reallocate(nb_points); 
              weights2d.Reallocate(nb_points); 
              FillZero(points2d); FillZero(weights2d);
              points2d(0).Init(0.0, 0.0);
              weights2d(0) = 0.02761042699769952;
              points2d(1).Init(-0.4977490260133565, 0.0);
              weights2d(1) = 0.001779029547326740;
              points2d(4).Init(-0.3587903720915737, 0.0);
              weights2d(4) = 0.02011239811396117;
              points2d(7).Init(-0.1932918138657104, 0.0);
              weights2d(7) = 0.02681784725933157;
              points2d(10).Init(0.2064993924016380, 0.0);
              weights2d(10) = 0.02452313380150201;
              points2d(13).Init(0.3669431077237697, 0.0);
              weights2d(13) = 0.01639457841069539;
              points2d(16).Init(0.6767931784861860, 0.0);       
              weights2d(16) = 0.01479590739864960;
              points2d(19).Init(0.8827927364865920, 0.0);
              weights2d(19) = 0.004579282277704251;
              points2d(22).Init(0.9664768608120111, 0.0);
              weights2d(22) = 0.001651826515576217;
              points2d(25).Init(-0.4919755727189941, -0.7513212483763635);
              weights2d(25) = 0.002349170908575584;
              points2d(31).Init(-0.4880677744007016, -0.5870191642967427);
              weights2d(31) = 0.004465925754181793;
              points2d(37).Init(-0.4843664025781043, -0.1717270984114328);
              weights2d(37) = 0.006099566807907972;
              points2d(43).Init(-0.4835533778058150, -0.3833898305784408);
              weights2d(43) = 0.006891081327188203;
              points2d(49).Init(-0.4421499318718065, -0.6563281974461070);
              weights2d(49) = 0.007997475072478163;
              points2d(55).Init(-0.4466292382741727, -0.06157647932662624);
              weights2d(55) = 0.007386134285336024;
              points2d(61).Init(-0.4254937754558538, -0.4783124082660027);
              weights2d(61) = 0.01279933187864826;
              points2d(67).Init(-0.4122204123735024, -0.2537089901614676);
              weights2d(67) = 0.01725807117569655;
              points2d(73).Init(-0.3177533194934086, -0.3996183176834929);
              weights2d(73) = 0.01867294590293547;
              points2d(79).Init(-0.2889337325840919, -0.1844183967233982);
              weights2d(79) = 0.02281822405839526;
              Real_wp coefx = 1.0/3, coefy = 1.0/sqrt(Real_wp(3));
              for (int i = 0; i < nb_points; i++)
                {
                  Real_wp x = points2d(i)(0);
                  Real_wp y = points2d(i)(1);
                  points2d(i).Init(coefx*(1.0-x) + coefy*y,
                                   coefx*(1.0-x) - coefy*y);
                  
                  weights2d(i) *= 0.5;
                }
            } 
            break; 	    
          default :
            formula_not_available = true;
          }
      }
    else if (type_quad == QUADRATURE_MASS_LUMPED)
      {
        switch (p)
          {
          case 0 :
            n1 = 0; n3 = 1; n6 = 0;
            nb_points = n1 + 3*n3 + 6*n6;
            points2d.Reallocate(nb_points); 
            weights2d.Reallocate(nb_points); 
            points2d(0).Init(0.0, 0.0);
            weights2d(0) = 1.0/6;
            break;
          case 2 :
            n1 = 1; n3 = 2; n6 = 0;
            nb_points = n1 + 3*n3 + 6*n6;
            points2d.Reallocate(nb_points); 
            weights2d.Reallocate(nb_points);             
            points2d(0).Init(1.0/3, 1.0/3);
            weights2d(0) = 0.225;
            points2d(1).Init(0.0, 0.0);
            weights2d(1) = 0.025;
            points2d(4).Init(0.5, 0.0);
            weights2d(4) = 1.0/15;
            break;
          default:
            {
              cout << "order not implemented " << endl;
              abort();
            }
          }
      }
    else if (type_quad == QUADRATURE_QUASI_LUMPED)
      {
        switch (p)
          {
          case 0:
            n1 = 0; n3 = 1; n6 = 0;
            nb_points = n1 + 3*n3 + 6*n6;
            points2d.Reallocate(nb_points); 
            weights2d.Reallocate(nb_points); 
            points2d(0).Init(0.0, 0.0);
            weights2d(0) = 1.0/6;
            break;
          case 2 :
            n1 = 1; n3 = 2; n6 = 0;
            nb_points = n1 + 3*n3 + 6*n6;
            points2d.Reallocate(nb_points); 
            weights2d.Reallocate(nb_points);             
            points2d(0).Init(1.0/3, 1.0/3);
            weights2d(0) = 0.225;
            points2d(1).Init(0.0, 0.0);
            weights2d(1) = 0.025;
            points2d(4).Init(0.5, 0.0);
            weights2d(4) = 1.0/15;
            break;
          case 4 :
            n1 = 1; n3 = 2; n6 = 1;
            nb_points = n1 + 3*n3 + 6*n6;
            points2d.Reallocate(nb_points); 
            weights2d.Reallocate(nb_points);             
            points2d(0).Init(1.0/3, 1.0/3);
            weights2d(0) = 0.267141408471108;
            points2d(1).Init(0.0, 0.0);
            weights2d(1) = 0.00757317419721512;
            points2d(4).Init(0.17, 0.66);
            weights2d(4) = 0.118355535056377;
            points2d(7).Init(0.0, 0.276393202250021);
            weights2d(7) = 0.0591787439613527;
            Mlt(Real_wp(0.5), weights2d);
            break;
          default:
            {
              cout << "order not implemented " << endl;
              abort();
            }
          }
      }
    else
      formula_not_available = true;
    
    if (formula_not_available)
      {
        int r = p/2;
        GetTensorizedQuadratureRule(r, points2d, weights2d, type_quad);
      }
    else
      {
        int nb_couples = n1 + n3 + n6;
        IVect NbPoints_couple(nb_couples);
        int ind = 0;
        for (int i = 0; i < n1; i++)
          NbPoints_couple(ind++) = 1;
        
        for (int i = 0; i < n3; i++)
          NbPoints_couple(ind++) = 3;
        
        for (int i = 0; i < n6; i++)
          NbPoints_couple(ind++) = 6;
        
        CompleteTrianglePointsWithSymmetry(points2d, weights2d, nb_couples, NbPoints_couple);
      }
  }
  
  
  void TriangleQuadrature
  ::GetTensorizedQuadratureRule(int r, VectR2& points2d, VectReal_wp& weights2d, int type_quad)
  {
    // other orders : we use Gauss-Radau formulas
    Globatto<Real_wp> gauss_x, gauss_y;
    gauss_x.ConstructQuadrature(r, gauss_x.QUADRATURE_GAUSS);
    gauss_y.ConstructQuadrature(r, gauss_x.QUADRATURE_JACOBI, 1.0, 0.0);
    
    int nb_points_quadrature = (r+1)*(r+1);
    points2d.Reallocate(nb_points_quadrature);
    weights2d.Reallocate(nb_points_quadrature);
    int nb = 0; Real_wp x, y, one(1);
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
        {
          x = gauss_x.Points(i);
          y = gauss_y.Points(j);
          points2d(nb)(0) = x*(one-y);
          points2d(nb)(1) = y;
          weights2d(nb) = gauss_x.Weights(i)*gauss_y.Weights(j);
          nb++;
	}
  }
  
  
  //! we assume that points2d and weights2d are partially filled with original points
  void TriangleQuadrature::
  CompleteTrianglePointsWithSymmetry(VectR2& points2d, VectReal_wp& weights2d,
                                     int nb_couples, const IVect& NbPoints_couple)
  {    
    Real_wp one(1), two(2);
    int ind = 0;
    // we generate the other points by symmetry
    for (int i = 0; i < nb_couples; i++)
      {
	int np = NbPoints_couple(i);
	
	if (np == 3)
	  {
            Real_wp x = points2d(ind)(0);
            Real_wp y = points2d(ind)(1);
            if (x != y)
	      {
		if (abs(two*x + y - one) > R2::threshold)
		  y = x;
              }
	    
            y = one - two*x;
            
            points2d(ind).Init(x, x);
	    points2d(ind+1).Init(y, x);
	    points2d(ind+2).Init(x, y);
	  }
	else if (np == 6)
	  {
            Real_wp x = points2d(ind)(0);
            Real_wp y = points2d(ind)(1);
	    Real_wp z = one - x - y;
            points2d(ind).Init(x, y);
	    points2d(ind+1).Init(y, x);
	    points2d(ind+2).Init(y, z);
	    points2d(ind+3).Init(z, y);
	    points2d(ind+4).Init(x, z);
	    points2d(ind+5).Init(z, x);
	  }
	
	for (int j = 1; j < np; j++)
	  weights2d(ind+j) = weights2d(ind);
	
	ind += np;
      }
  }
  
}

#define MONTJOIE_FILE_TRIANGLE_QUADRATURE_CXX
#endif
  
