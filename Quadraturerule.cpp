#include <vector>
#include <iostream>
#include <algorithm>

#include "Quadraturerule.H"

// return gauss quadrature information for line integral, first index in the gauss quadrature coordinate number, second index runs from 0 to 1. 0 is xi coordinate. 1 is weight.
std::vector<std::vector<double>> gauss_line_integral(const int& IntOrder){

  // Order 1 Legendre-Gauss Points
  int n1 = 1;
  double x1[] = {
    0.500000000000000
  };
  double w1[] = {
    1.000000000000000
  };
  
  // Order 3 Legendre-Gauss points
  int n3 = 2;
  double x3[] = {
    0.211324865405187, 0.788675134594813
  };
  double w3[] = {
    0.500000000000000, 0.500000000000000
  };

  // Order 5 Legendre-Gauss points
  int n5 = 3;
  double x5[] = {
    0.112701665379258, 0.500000000000000, 0.887298334620742
  };
  double w5[] = {
    0.277777777777778, 0.444444444444444, 0.277777777777778
  };

  // Order 7 Legendre-Gauss points
  int n7 = 4;
  double x7[] = {
    0.069431844202974, 0.330009478207572, 0.669990521792428, 0.930568155797026
  };
  double w7[] = {
    0.173927422568727, 0.326072577431273, 0.326072577431273, 0.173927422568727
  };

  // Order 9 Legendre-Gauss points
  int n9 = 5;
  double x9[] = {
    0.046910077030668, 0.230765344947158, 0.500000000000000, 0.769234655052841,
    0.953089922969332
  };
  double w9[] = {
    0.118463442528095, 0.239314335249683, 0.284444444444444, 0.239314335249683,
    0.118463442528095
  };

  // Order 11 Legendre-Gauss points
  int n11 = 6;
  double x11[] = {
    0.033765242898424, 0.169395306766868, 0.380690406958402, 0.619309593041598,
    0.830604693233132, 0.966234757101576
  };
  double w11[] = {
    0.085662246189585, 0.180380786524069, 0.233956967286345, 0.233956967286345,
    0.180380786524069, 0.085662246189585
  };

  // Order 13 Legendre-Gauss points
  int n13 = 7;
  double x13[] = {
    0.025446043828621, 0.129234407200303, 0.297077424311301, 0.500000000000000,
    0.702922575688699, 0.870765592799697, 0.974553956171379
  };
  double w13[] = {
    0.064742483084435, 0.139852695744638, 0.190915025252560, 0.208979591836735,
    0.190915025252560, 0.139852695744638, 0.064742483084435
  };

  // Order 15 Legendre-Gauss points
  int n15 = 8;
  double x15[] = {
    0.019855071751232, 0.101666761293187, 0.237233795041836, 0.408282678752175,
    0.591717321247825, 0.762766204958164, 0.898333238706813, 0.980144928248768
  };
  double w15[] = {
    0.050614268145188, 0.111190517226687, 0.156853322938944, 0.181341891689181,
    0.181341891689181, 0.156853322938944, 0.111190517226687, 0.050614268145188
  };

  // Order 17 Legendre-Gauss points
  int n17 = 9;
  double x17[] = {
    0.015919880246187, 0.081984446336682, 0.193314283649705, 0.337873288298096,
    0.500000000000000, 0.662126711701905, 0.806685716350295, 0.918015553663318,
    0.984080119753813
  };
  double w17[] = {
    0.040637194180787, 0.090324080347429, 0.130305348201468, 0.156173538520001,
    0.165119677500630, 0.156173538520001, 0.130305348201468, 0.090324080347429,
    0.040637194180787
  };

  // Order 19 Legendre-Gauss points
  int n19 = 10;
  double x19[] = {
    0.013046735741414, 0.067468316655508, 0.160295215850488, 0.283302302935376,
    0.425562830509184, 0.574437169490816, 0.716697697064624, 0.839704784149512,
    0.932531683344492, 0.986953264258586
  };
  double w19[] = {
    0.033335672154344, 0.074725674575290, 0.109543181257991, 0.134633359654998,
    0.147762112357376, 0.147762112357376, 0.134633359654998, 0.109543181257991,
    0.074725674575290, 0.033335672154344
  };

  // with N points, we can integrate accurately up to order 2N-1
  // (start with IntOrder = 0)
  int IntOrder2Rule[] = {1,1,3,3,5,5,7,7,9,9,11,11,13,13,15,15,17,17,19,19};

    int n;
    double x[IntOrder2Rule[IntOrder]];
    double w[IntOrder2Rule[IntOrder]];

  // choose necessary quadrature rule
  switch(IntOrder2Rule[IntOrder]) {
  case 1:
    n = n1; std::copy(x1, x1 + n1, x); std::copy(w1, w1 + n1, w);
    break;
  case 3:
    n = n3; std::copy(x3, x3 + n3, x); std::copy(w3, w3 + n3, w);
    break;
  case 5:
    n = n5; std::copy(x5, x5 + n5, x); std::copy(w5, w5 + n5, w);
    break;
  case 7:
    n = n7; std::copy(x7, x7 + n7, x); std::copy(w7, w7 + n7, w);
    break;
  case 9:
    n = n9; std::copy(x9, x9 + n9, x); std::copy(w9, w9 + n9, w);
    break;
  case 11:
    n = n11; std::copy(x11, x11 + n11, x); std::copy(w11, w11 + n11, w);
    break;
  case 13:
    n = n13; std::copy(x13, x13 + n13, x); std::copy(w13, w13 + n13, w);
    break;
  case 15:
    n = n15; std::copy(x15, x15 + n15, x); std::copy(w15, w15 + n15, w);
    break;
  case 17:
    n = n17; std::copy(x17, x17 + n17, x); std::copy(w17, w17 + n17, w);
    break;
  case 19:
    n = n19; std::copy(x19, x19 + n19, x); std::copy(w19, w19 + n19, w);
    break;
  default:
    printf("ERROR: Unsupported quadrature order \n");
    exit(EXIT_FAILURE);
    break;
  }
    
    // return gauss quadrature information for line integral, first index in the gauss quadrature coordinate number, second index runs from 0 to 1. 0 is xi coordinate. 1 is weight.
    std::vector<std::vector<double>> gauss_quadrature_weights(n, std::vector<double>(2));

    for (int i = 0; i < n ; ++i) {
        gauss_quadrature_weights[i][0] = x[ i ]; // xi coordinates of gauss quadrature
        gauss_quadrature_weights[i][1] = w[ i ]; // weight
    }

    // return gauss cuadrature information
    return gauss_quadrature_weights;

}

std::vector<std::vector<double>> gauss_area_integral(const int& IntOrder){

 /*
    Dunavant points generated with .m code written by John Burkard

    http://people.scs.fsu.edu/~burkardt/f_src/dunavant/dunavant.html
    
   1. David Dunavant,
      High Degree Efficient Symmetrical Gaussian Quadrature Rules for the Triangle,
      International Journal for Numerical Methods in Engineering,
      Volume 21, 1985, pages 1129-1148.
   2. James Lyness, Dennis Jespersen,
      Moderate Degree Symmetric Quadrature Rules for the Triangle,
      Journal of the Institute of Mathematics and its Applications,
      Volume 15, Number 1, February 1975, pages 19-32.

   Note, the coordinates in the x[] vectors are stored unrolled in
   sequential x,y pairs, e.g. [x1 y1  x2 y2  x3 y3 ... ]
   
   This is given by Prof. Fidkowski.

  */

//   Order 1 Dunavant Points
  int n1 = 1;
  double x1[] = {
    0.333333333333333, 0.333333333333333
  };
  double w1[] = {
    0.500000000000000
  };

  // Order 2 Dunavant Points
  int n2 = 3;
  double x2[] = {
    0.666666666666667, 0.166666666666667, 0.166666666666667, 0.166666666666667,
    0.166666666666667, 0.666666666666667
  };
  double w2[] = {
    0.166666666666666, 0.166666666666666, 0.166666666666666
  };

  // Order 3 Dunavant Points
  int n3 = 4;
  double x3[] = {
    0.333333333333333, 0.333333333333333, 0.600000000000000, 0.200000000000000,
    0.200000000000000, 0.200000000000000, 0.200000000000000, 0.600000000000000
  };
  double w3[] = {
    -0.281250000000000, 0.260416666666667, 0.260416666666667, 0.260416666666667
  };

  // Order 4 Dunavant Points
  int n4 = 6;
  double x4[] = {
    0.108103018168070, 0.445948490915965, 0.445948490915965, 0.445948490915965,
    0.445948490915965, 0.108103018168070, 0.816847572980459, 0.091576213509771,
    0.091576213509771, 0.091576213509771, 0.091576213509771, 0.816847572980459
  };
  double w4[] = {
    0.111690794839005, 0.111690794839005, 0.111690794839005, 0.054975871827661,
    0.054975871827661, 0.054975871827661
  };

  // Order 5 Dunavant Points
  int n5 = 7;
  double x5[] = {
    0.333333333333333, 0.333333333333333, 0.059715871789770, 0.470142064105115,
    0.470142064105115, 0.470142064105115, 0.470142064105115, 0.059715871789770,
    0.797426985353087, 0.101286507323456, 0.101286507323456, 0.101286507323456,
    0.101286507323456, 0.797426985353087
  };
  double w5[] = {
    0.112500000000000, 0.066197076394253, 0.066197076394253, 0.066197076394253,
    0.062969590272414, 0.062969590272414, 0.062969590272414
  };

  // Order 6 Dunavant Points
  int n6 = 12;
  double x6[] = {
    0.501426509658179, 0.249286745170910, 0.249286745170910, 0.249286745170910,
    0.249286745170910, 0.501426509658179, 0.873821971016996, 0.063089014491502,
    0.063089014491502, 0.063089014491502, 0.063089014491502, 0.873821971016996,
    0.053145049844817, 0.310352451033784, 0.310352451033784, 0.636502499121399,
    0.636502499121399, 0.053145049844817, 0.310352451033784, 0.053145049844817,
    0.636502499121399, 0.310352451033784, 0.053145049844817, 0.636502499121399
  };
  double w6[] = {
    0.058393137863189, 0.058393137863189, 0.058393137863189, 0.025422453185103,
    0.025422453185103, 0.025422453185103, 0.041425537809187, 0.041425537809187,
    0.041425537809187, 0.041425537809187, 0.041425537809187, 0.041425537809187
  };

  // Order 7 Dunavant Points
  int n7 = 13;
  double x7[] = {
    0.333333333333333, 0.333333333333333, 0.479308067841920, 0.260345966079040,
    0.260345966079040, 0.260345966079040, 0.260345966079040, 0.479308067841920,
    0.869739794195568, 0.065130102902216, 0.065130102902216, 0.065130102902216,
    0.065130102902216, 0.869739794195568, 0.048690315425316, 0.312865496004874,
    0.312865496004874, 0.638444188569810, 0.638444188569810, 0.048690315425316,
    0.312865496004874, 0.048690315425316, 0.638444188569810, 0.312865496004874,
    0.048690315425316, 0.638444188569810
  };
  double w7[] = {
    -0.074785022233841, 0.087807628716604, 0.087807628716604, 0.087807628716604,
    0.026673617804419, 0.026673617804419, 0.026673617804419, 0.038556880445128,
    0.038556880445128, 0.038556880445128, 0.038556880445128, 0.038556880445128,
    0.038556880445128
  };

  // Order 8 Dunavant Points
  int n8 = 16;
  double x8[] = {
    0.333333333333333, 0.333333333333333, 0.081414823414554, 0.459292588292723,
    0.459292588292723, 0.459292588292723, 0.459292588292723, 0.081414823414554,
    0.658861384496480, 0.170569307751760, 0.170569307751760, 0.170569307751760,
    0.170569307751760, 0.658861384496480, 0.898905543365938, 0.050547228317031,
    0.050547228317031, 0.050547228317031, 0.050547228317031, 0.898905543365938,
    0.008394777409958, 0.263112829634638, 0.263112829634638, 0.728492392955404,
    0.728492392955404, 0.008394777409958, 0.263112829634638, 0.008394777409958,
    0.728492392955404, 0.263112829634638, 0.008394777409958, 0.728492392955404
  };
  double w8[] = {
    0.072157803838894, 0.047545817133642, 0.047545817133642, 0.047545817133642,
    0.051608685267359, 0.051608685267359, 0.051608685267359, 0.016229248811599,
    0.016229248811599, 0.016229248811599, 0.013615157087217, 0.013615157087217,
    0.013615157087217, 0.013615157087217, 0.013615157087217, 0.013615157087217
  };

  // Order 9 Dunavant Points
  int n9 = 19;
  double x9[] = {
    0.333333333333333, 0.333333333333333, 0.020634961602525, 0.489682519198738,
    0.489682519198738, 0.489682519198738, 0.489682519198738, 0.020634961602525,
    0.125820817014127, 0.437089591492937, 0.437089591492937, 0.437089591492937,
    0.437089591492937, 0.125820817014127, 0.623592928761935, 0.188203535619033,
    0.188203535619033, 0.188203535619033, 0.188203535619033, 0.623592928761935,
    0.910540973211095, 0.044729513394453, 0.044729513394453, 0.044729513394453,
    0.044729513394453, 0.910540973211095, 0.036838412054736, 0.221962989160766,
    0.221962989160766, 0.741198598784498, 0.741198598784498, 0.036838412054736,
    0.221962989160766, 0.036838412054736, 0.741198598784498, 0.221962989160766,
    0.036838412054736, 0.741198598784498
  };
  double w9[] = {
    0.048567898141400, 0.015667350113570, 0.015667350113570, 0.015667350113570,
    0.038913770502387, 0.038913770502387, 0.038913770502387, 0.039823869463605,
    0.039823869463605, 0.039823869463605, 0.012788837829349, 0.012788837829349,
    0.012788837829349, 0.021641769688645, 0.021641769688645, 0.021641769688645,
    0.021641769688645, 0.021641769688645, 0.021641769688645
  };

  // Order 10 Dunavant Points
  int n10 = 25;
  double x10[] = {
    0.333333333333333, 0.333333333333333, 0.028844733232685, 0.485577633383657,
    0.485577633383657, 0.485577633383657, 0.485577633383657, 0.028844733232685,
    0.781036849029926, 0.109481575485037, 0.109481575485037, 0.109481575485037,
    0.109481575485037, 0.781036849029926, 0.141707219414880, 0.307939838764121,
    0.307939838764121, 0.550352941820999, 0.550352941820999, 0.141707219414880,
    0.307939838764121, 0.141707219414880, 0.550352941820999, 0.307939838764121,
    0.141707219414880, 0.550352941820999, 0.025003534762686, 0.246672560639903,
    0.246672560639903, 0.728323904597411, 0.728323904597411, 0.025003534762686,
    0.246672560639903, 0.025003534762686, 0.728323904597411, 0.246672560639903,
    0.025003534762686, 0.728323904597411, 0.009540815400299, 0.066803251012200,
    0.066803251012200, 0.923655933587500, 0.923655933587500, 0.009540815400299,
    0.066803251012200, 0.009540815400299, 0.923655933587500, 0.066803251012200,
    0.009540815400299, 0.923655933587500
  };
  double w10[] = {
    0.045408995191377, 0.018362978878233, 0.018362978878233, 0.018362978878233,
    0.022660529717764, 0.022660529717764, 0.022660529717764, 0.036378958422710,
    0.036378958422710, 0.036378958422710, 0.036378958422710, 0.036378958422710,
    0.036378958422710, 0.014163621265528, 0.014163621265528, 0.014163621265528,
    0.014163621265528, 0.014163621265528, 0.014163621265528, 0.004710833481867,
    0.004710833481867, 0.004710833481867, 0.004710833481867, 0.004710833481867,
    0.004710833481867
  };

  // Order 12 Dunavant Points
  int n12 = 33;
  double x12[] = {
    0.023565220452390, 0.488217389773805, 0.488217389773805, 0.488217389773805,
    0.488217389773805, 0.023565220452390, 0.120551215411079, 0.439724392294460,
    0.439724392294460, 0.439724392294460, 0.439724392294460, 0.120551215411079,
    0.457579229975768, 0.271210385012116, 0.271210385012116, 0.271210385012116,
    0.271210385012116, 0.457579229975768, 0.744847708916828, 0.127576145541586,
    0.127576145541586, 0.127576145541586, 0.127576145541586, 0.744847708916828,
    0.957365299093579, 0.021317350453210, 0.021317350453210, 0.021317350453210,
    0.021317350453210, 0.957365299093579, 0.115343494534698, 0.275713269685514,
    0.275713269685514, 0.608943235779788, 0.608943235779788, 0.115343494534698,
    0.275713269685514, 0.115343494534698, 0.608943235779788, 0.275713269685514,
    0.115343494534698, 0.608943235779788, 0.022838332222257, 0.281325580989940,
    0.281325580989940, 0.695836086787803, 0.695836086787803, 0.022838332222257,
    0.281325580989940, 0.022838332222257, 0.695836086787803, 0.281325580989940,
    0.022838332222257, 0.695836086787803, 0.025734050548330, 0.116251915907597,
    0.116251915907597, 0.858014033544073, 0.858014033544073, 0.025734050548330,
    0.116251915907597, 0.025734050548330, 0.858014033544073, 0.116251915907597,
    0.025734050548330, 0.858014033544073
  };
  double w12[] = {
    0.012865533220227, 0.012865533220227, 0.012865533220227, 0.021846272269019,
    0.021846272269019, 0.021846272269019, 0.031429112108943, 0.031429112108943,
    0.031429112108943, 0.017398056465355, 0.017398056465355, 0.017398056465355,
    0.003083130525780, 0.003083130525780, 0.003083130525780, 0.020185778883191,
    0.020185778883191, 0.020185778883191, 0.020185778883191, 0.020185778883191,
    0.020185778883191, 0.011178386601152, 0.011178386601152, 0.011178386601152,
    0.011178386601152, 0.011178386601152, 0.011178386601152, 0.008658115554329,
    0.008658115554329, 0.008658115554329, 0.008658115554329, 0.008658115554329,
    0.008658115554329
  };

  // Order 13 Dunavant Points
  int n13 = 37;
  double x13[] = {
    0.333333333333333, 0.333333333333333, 0.009903630120591, 0.495048184939705,
    0.495048184939705, 0.495048184939705, 0.495048184939705, 0.009903630120591,
    0.062566729780852, 0.468716635109574, 0.468716635109574, 0.468716635109574,
    0.468716635109574, 0.062566729780852, 0.170957326397447, 0.414521336801277,
    0.414521336801277, 0.414521336801277, 0.414521336801277, 0.170957326397447,
    0.541200855914337, 0.229399572042831, 0.229399572042831, 0.229399572042831,
    0.229399572042831, 0.541200855914337, 0.771151009607340, 0.114424495196330,
    0.114424495196330, 0.114424495196330, 0.114424495196330, 0.771151009607340,
    0.950377217273082, 0.024811391363459, 0.024811391363459, 0.024811391363459,
    0.024811391363459, 0.950377217273082, 0.094853828379579, 0.268794997058761,
    0.268794997058761, 0.636351174561660, 0.636351174561660, 0.094853828379579,
    0.268794997058761, 0.094853828379579, 0.636351174561660, 0.268794997058761,
    0.094853828379579, 0.636351174561660, 0.018100773278807, 0.291730066734288,
    0.291730066734288, 0.690169159986905, 0.690169159986905, 0.018100773278807,
    0.291730066734288, 0.018100773278807, 0.690169159986905, 0.291730066734288,
    0.018100773278807, 0.690169159986905, 0.022233076674090, 0.126357385491669,
    0.126357385491669, 0.851409537834241, 0.851409537834241, 0.022233076674090,
    0.126357385491669, 0.022233076674090, 0.851409537834241, 0.126357385491669,
    0.022233076674090, 0.851409537834241
  };
  double w13[] = {
    0.026260461700401, 0.005640072604665, 0.005640072604665, 0.005640072604665,
    0.015711759181227, 0.015711759181227, 0.015711759181227, 0.023536251252097,
    0.023536251252097, 0.023536251252097, 0.023681793268178, 0.023681793268178,
    0.023681793268178, 0.015583764522897, 0.015583764522897, 0.015583764522897,
    0.003987885732537, 0.003987885732537, 0.003987885732537, 0.018424201364366,
    0.018424201364366, 0.018424201364366, 0.018424201364366, 0.018424201364366,
    0.018424201364366, 0.008700731651911, 0.008700731651911, 0.008700731651911,
    0.008700731651911, 0.008700731651911, 0.008700731651911, 0.007760893419522,
    0.007760893419522, 0.007760893419522, 0.007760893419522, 0.007760893419522,
    0.007760893419522
  };

  // Order 14 Dunavant Points
  int n14 = 42;
  double x14[] = {
    0.022072179275643, 0.488963910362179, 0.488963910362179, 0.488963910362179,
    0.488963910362179, 0.022072179275643, 0.164710561319092, 0.417644719340454,
    0.417644719340454, 0.417644719340454, 0.417644719340454, 0.164710561319092,
    0.453044943382323, 0.273477528308839, 0.273477528308839, 0.273477528308839,
    0.273477528308839, 0.453044943382323, 0.645588935174913, 0.177205532412543,
    0.177205532412543, 0.177205532412543, 0.177205532412543, 0.645588935174913,
    0.876400233818255, 0.061799883090873, 0.061799883090873, 0.061799883090873,
    0.061799883090873, 0.876400233818255, 0.961218077502598, 0.019390961248701,
    0.019390961248701, 0.019390961248701, 0.019390961248701, 0.961218077502598,
    0.057124757403648, 0.172266687821356, 0.172266687821356, 0.770608554774996,
    0.770608554774996, 0.057124757403648, 0.172266687821356, 0.057124757403648,
    0.770608554774996, 0.172266687821356, 0.057124757403648, 0.770608554774996,
    0.092916249356972, 0.336861459796345, 0.336861459796345, 0.570222290846683,
    0.570222290846683, 0.092916249356972, 0.336861459796345, 0.092916249356972,
    0.570222290846683, 0.336861459796345, 0.092916249356972, 0.570222290846683,
    0.014646950055654, 0.298372882136258, 0.298372882136258, 0.686980167808088,
    0.686980167808088, 0.014646950055654, 0.298372882136258, 0.014646950055654,
    0.686980167808088, 0.298372882136258, 0.014646950055654, 0.686980167808088,
    0.001268330932872, 0.118974497696957, 0.118974497696957, 0.879757171370171,
    0.879757171370171, 0.001268330932872, 0.118974497696957, 0.001268330932872,
    0.879757171370171, 0.118974497696957, 0.001268330932872, 0.879757171370171
  };
  double w14[] = {
    0.010941790684715, 0.010941790684715, 0.010941790684715, 0.016394176772063,
    0.016394176772063, 0.016394176772063, 0.025887052253646, 0.025887052253646,
    0.025887052253646, 0.021081294368497, 0.021081294368497, 0.021081294368497,
    0.007216849834889, 0.007216849834889, 0.007216849834889, 0.002461701801200,
    0.002461701801200, 0.002461701801200, 0.012332876606282, 0.012332876606282,
    0.012332876606282, 0.012332876606282, 0.012332876606282, 0.012332876606282,
    0.019285755393531, 0.019285755393531, 0.019285755393531, 0.019285755393531,
    0.019285755393531, 0.019285755393531, 0.007218154056767, 0.007218154056767,
    0.007218154056767, 0.007218154056767, 0.007218154056767, 0.007218154056767,
    0.002505114419250, 0.002505114419250, 0.002505114419250, 0.002505114419250,
    0.002505114419250, 0.002505114419250
  };

  // Order 17 Dunavant Points
  int n17 = 61;
  double x17[] = {
    0.333333333333333, 0.333333333333333, 0.005658918886452, 0.497170540556774,
    0.497170540556774, 0.497170540556774, 0.497170540556774, 0.005658918886452,
    0.035647354750751, 0.482176322624625, 0.482176322624625, 0.482176322624625,
    0.482176322624625, 0.035647354750751, 0.099520061958437, 0.450239969020782,
    0.450239969020782, 0.450239969020782, 0.450239969020782, 0.099520061958437,
    0.199467521245206, 0.400266239377397, 0.400266239377397, 0.400266239377397,
    0.400266239377397, 0.199467521245206, 0.495717464058095, 0.252141267970953,
    0.252141267970953, 0.252141267970953, 0.252141267970953, 0.495717464058095,
    0.675905990683077, 0.162047004658461, 0.162047004658461, 0.162047004658461,
    0.162047004658461, 0.675905990683077, 0.848248235478508, 0.075875882260746,
    0.075875882260746, 0.075875882260746, 0.075875882260746, 0.848248235478508,
    0.968690546064356, 0.015654726967822, 0.015654726967822, 0.015654726967822,
    0.015654726967822, 0.968690546064356, 0.010186928826919, 0.334319867363658,
    0.334319867363658, 0.655493203809423, 0.655493203809423, 0.010186928826919,
    0.334319867363658, 0.010186928826919, 0.655493203809423, 0.334319867363658,
    0.010186928826919, 0.655493203809423, 0.135440871671036, 0.292221537796944,
    0.292221537796944, 0.572337590532020, 0.572337590532020, 0.135440871671036,
    0.292221537796944, 0.135440871671036, 0.572337590532020, 0.292221537796944,
    0.135440871671036, 0.572337590532020, 0.054423924290583, 0.319574885423190,
    0.319574885423190, 0.626001190286228, 0.626001190286228, 0.054423924290583,
    0.319574885423190, 0.054423924290583, 0.626001190286228, 0.319574885423190,
    0.054423924290583, 0.626001190286228, 0.012868560833637, 0.190704224192292,
    0.190704224192292, 0.796427214974071, 0.796427214974071, 0.012868560833637,
    0.190704224192292, 0.012868560833637, 0.796427214974071, 0.190704224192292,
    0.012868560833637, 0.796427214974071, 0.067165782413524, 0.180483211648746,
    0.180483211648746, 0.752351005937729, 0.752351005937729, 0.067165782413524,
    0.180483211648746, 0.067165782413524, 0.752351005937729, 0.180483211648746,
    0.067165782413524, 0.752351005937729, 0.014663182224828, 0.080711313679564,
    0.080711313679564, 0.904625504095608, 0.904625504095608, 0.014663182224828,
    0.080711313679564, 0.014663182224828, 0.904625504095608, 0.080711313679564,
    0.014663182224828, 0.904625504095608
  };
  double w17[] = {
    0.016718599645402, 0.002546707720254, 0.002546707720254, 0.002546707720254,
    0.007335432263819, 0.007335432263819, 0.007335432263819, 0.012175439176836,
    0.012175439176836, 0.012175439176836, 0.015553775434484, 0.015553775434484,
    0.015553775434484, 0.015628555609310, 0.015628555609310, 0.015628555609310,
    0.012407827169833, 0.012407827169833, 0.012407827169833, 0.007028036535279,
    0.007028036535279, 0.007028036535279, 0.001597338086889, 0.001597338086889,
    0.001597338086889, 0.004059827659497, 0.004059827659497, 0.004059827659497,
    0.004059827659497, 0.004059827659497, 0.004059827659497, 0.013402871141582,
    0.013402871141582, 0.013402871141582, 0.013402871141582, 0.013402871141582,
    0.013402871141582, 0.009229996605411, 0.009229996605411, 0.009229996605411,
    0.009229996605411, 0.009229996605411, 0.009229996605411, 0.004238434267164,
    0.004238434267164, 0.004238434267164, 0.004238434267164, 0.004238434267164,
    0.004238434267164, 0.009146398385012, 0.009146398385012, 0.009146398385012,
    0.009146398385012, 0.009146398385012, 0.009146398385012, 0.003332816002083,
    0.003332816002083, 0.003332816002083, 0.003332816002083, 0.003332816002083,
    0.003332816002083
  };

  // Order 19 Dunavant Points
  int n19 = 73;
  double x19[] = {
    0.333333333333333, 0.333333333333333, 0.020780025853987, 0.489609987073006,
    0.489609987073006, 0.489609987073006, 0.489609987073006, 0.020780025853987,
    0.090926214604215, 0.454536892697893, 0.454536892697893, 0.454536892697893,
    0.454536892697893, 0.090926214604215, 0.197166638701138, 0.401416680649431,
    0.401416680649431, 0.401416680649431, 0.401416680649431, 0.197166638701138,
    0.488896691193805, 0.255551654403098, 0.255551654403098, 0.255551654403098,
    0.255551654403098, 0.488896691193805, 0.645844115695741, 0.177077942152130,
    0.177077942152130, 0.177077942152130, 0.177077942152130, 0.645844115695741,
    0.779877893544096, 0.110061053227952, 0.110061053227952, 0.110061053227952,
    0.110061053227952, 0.779877893544096, 0.888942751496321, 0.055528624251840,
    0.055528624251840, 0.055528624251840, 0.055528624251840, 0.888942751496321,
    0.974756272445543, 0.012621863777229, 0.012621863777229, 0.012621863777229,
    0.012621863777229, 0.974756272445543, 0.003611417848412, 0.395754787356943,
    0.395754787356943, 0.600633794794645, 0.600633794794645, 0.003611417848412,
    0.395754787356943, 0.003611417848412, 0.600633794794645, 0.395754787356943,
    0.003611417848412, 0.600633794794645, 0.134466754530780, 0.307929983880436,
    0.307929983880436, 0.557603261588784, 0.557603261588784, 0.134466754530780,
    0.307929983880436, 0.134466754530780, 0.557603261588784, 0.307929983880436,
    0.134466754530780, 0.557603261588784, 0.014446025776115, 0.264566948406520,
    0.264566948406520, 0.720987025817365, 0.720987025817365, 0.014446025776115,
    0.264566948406520, 0.014446025776115, 0.720987025817365, 0.264566948406520,
    0.014446025776115, 0.720987025817365, 0.046933578838178, 0.358539352205951,
    0.358539352205951, 0.594527068955871, 0.594527068955871, 0.046933578838178,
    0.358539352205951, 0.046933578838178, 0.594527068955871, 0.358539352205951,
    0.046933578838178, 0.594527068955871, 0.002861120350567, 0.157807405968595,
    0.157807405968595, 0.839331473680839, 0.839331473680839, 0.002861120350567,
    0.157807405968595, 0.002861120350567, 0.839331473680839, 0.157807405968595,
    0.002861120350567, 0.839331473680839, 0.223861424097916, 0.075050596975911,
    0.075050596975911, 0.701087978926173, 0.701087978926173, 0.223861424097916,
    0.075050596975911, 0.223861424097916, 0.701087978926173, 0.075050596975911,
    0.223861424097916, 0.701087978926173, 0.034647074816760, 0.142421601113383,
    0.142421601113383, 0.822931324069857, 0.822931324069857, 0.034647074816760,
    0.142421601113383, 0.034647074816760, 0.822931324069857, 0.142421601113383,
    0.034647074816760, 0.822931324069857, 0.010161119296278, 0.065494628082938,
    0.065494628082938, 0.924344252620784, 0.924344252620784, 0.010161119296278,
    0.065494628082938, 0.010161119296278, 0.924344252620784, 0.065494628082938,
    0.010161119296278, 0.924344252620784
  };
  double w19[] = {
    0.016453165694459, 0.005165365945636, 0.005165365945636, 0.005165365945636,
    0.011193623631508, 0.011193623631508, 0.011193623631508, 0.015133062934734,
    0.015133062934734, 0.015133062934734, 0.015245483901099, 0.015245483901099,
    0.015245483901099, 0.012079606370821, 0.012079606370821, 0.012079606370821,
    0.008025401793400, 0.008025401793400, 0.008025401793400, 0.004042290130892,
    0.004042290130892, 0.004042290130892, 0.001039681013742, 0.001039681013742,
    0.001039681013742, 0.001942438452491, 0.001942438452491, 0.001942438452491,
    0.001942438452491, 0.001942438452491, 0.001942438452491, 0.012787080306011,
    0.012787080306011, 0.012787080306011, 0.012787080306011, 0.012787080306011,
    0.012787080306011, 0.004440451786669, 0.004440451786669, 0.004440451786669,
    0.004440451786669, 0.004440451786669, 0.004440451786669, 0.008062273380866,
    0.008062273380866, 0.008062273380866, 0.008062273380866, 0.008062273380866,
    0.008062273380866, 0.001245970908745, 0.001245970908745, 0.001245970908745,
    0.001245970908745, 0.001245970908745, 0.001245970908745, 0.009121420059476,
    0.009121420059476, 0.009121420059476, 0.009121420059476, 0.009121420059476,
    0.009121420059476, 0.005129281868099, 0.005129281868099, 0.005129281868099,
    0.005129281868099, 0.005129281868099, 0.005129281868099, 0.001899964427651,
    0.001899964427651, 0.001899964427651, 0.001899964427651, 0.001899964427651,
    0.001899964427651
  };
  
  // (start with IntOrder = 0)
  int IntOrder2Rule[] = {1,1,2,3,4,5,6,7,8,9,10,12,12,13,14,17,17,17,19,19};

    int n;
    double x[2 * IntOrder2Rule[IntOrder]];
    double w[IntOrder2Rule[IntOrder]];

  // choose necessary quadrature rule
  switch(IntOrder2Rule[IntOrder]) {
  case 1:
    n = n1; std::copy(x1, x1 + n1, x); std::copy(w1, w1 + n1, w);
    break;
  case 2:
    n = n2; std::copy(x2, x2 + n2, x); std::copy(w2, w2 + n2, w);
    break;
  case 3:
    n = n3; std::copy(x3, x3 + n3, x); std::copy(w3, w3 + n3, w);
    break;
  case 4:
    n = n4; std::copy(x4, x4 + n4, x); std::copy(w4, w4 + n4, w);
    break;
  case 5:
    n = n5; std::copy(x5, x5 + n5, x); std::copy(w5, w5 + n5, w);
    break;
  case 6:
    n = n6; std::copy(x6, x6 + n6, x); std::copy(w6, w6 + n6, w);
    break;
  case 7:
    n = n7; std::copy(x7, x7 + n7, x); std::copy(w7, w7 + n7, w);
    break;
  case 8:
    n = n8; std::copy(x8, x8 + n8, x); std::copy(w8, w8 + n8, w);
    break;
  case 9:
    n = n9; std::copy(x9, x9 + n9, x); std::copy(w9, w9 + n9, w);
    break;
  case 10:
    n = n10; std::copy(x10, x10 + n10, x); std::copy(w10, w10 + n10, w);
    break;
  case 12:
    n = n12; std::copy(x12, x12 + n12, x); std::copy(w12, w12 + n12, w);
    break;
  case 13:
    n = n13; std::copy(x13, x13 + n13, x); std::copy(w13, w13 + n13, w);
    break;
  case 14:
    n = n14; std::copy(x14, x14 + n14, x); std::copy(w14, w14 + n14, w);
    break;
  case 17:
    n = n17; std::copy(x17, x17 + n17, x); std::copy(w17, w17 + n17, w);
    break;
  case 19:
    n = n19; std::copy(x19, x19 + n19, x); std::copy(w19, w19 + n19, w);
    break;
  default:
    printf("ERROR: Unsupported quadrature order \n");
    exit(EXIT_FAILURE);
    break;
  }

    // save gauss quadrature information for area integral, first index in the gauss quadrature coordinate number, second index runs from 0 to 2. 0 is xi coordinate. 1 is eta coordinate. 2 is weight.
    std::vector<std::vector<double>> gauss_quadrature_weights(n, std::vector<double>(3));

    for (int i = 0; i < n ; ++i) {
        gauss_quadrature_weights[i][0] = x[ 2 * i ]; // xi coordinates of gauss quadrature
        gauss_quadrature_weights[i][1] = x[ 2 * i + 1]; // eta coordinates of gauss quadrature
        gauss_quadrature_weights[i][2] = w[ 1 ]; // weight
    }

    // return gauss cuadrature information
    return gauss_quadrature_weights;

}