#include "CMGTools/RootTools/interface/RochCor2012.h"
#include <TLorentzVector.h>


const double RochCor2012::pi = 3.14159265358979323846;
  
const float RochCor2012::genm_smr = 9.09915e+01; //gen mass peak with eta dependent gaussian smearing => better match in Z mass profile vs. eta/phi
const float RochCor2012::genm = 91.06; //gen mass peak without smearing => Z mass profile vs. eta/phi in CMS note
  
const float RochCor2012::mrecm = 90.8177; //rec mass peak in MC (2011A)
const float RochCor2012::drecm = 90.5332; //rec mass peak in data (2011A)
const float RochCor2012::mgscl_stat = 0.0001; //stat. error of global factor for mass peak in MC (2011A)  
const float RochCor2012::mgscl_syst = 0.0006; //syst. error of global factor for mass peak in MC (2011A)  
const float RochCor2012::dgscl_stat = 0.0001; //stat. error of global factor for mass peak in data (2011A)
const float RochCor2012::dgscl_syst = 0.0008; //syst. error of global factor for mass peak in data (2011A)
  
  //iteration2 after FSR : after Z Pt correction
const float RochCor2012::delta = -5.48477e-06;
const float RochCor2012::delta_stat = 4.38582e-07;
const float RochCor2012::delta_syst = 6.992e-07;
  
const float RochCor2012::sf = 33.4956;
const float RochCor2012::sf_stat = 0.312614;
const float RochCor2012::sf_syst = 9.29;
  
const float RochCor2012::apar = 1.0; //+- 0.002
const float RochCor2012::bpar = -5.03313e-06; //+- 1.57968e-06
const float RochCor2012::cpar = -4.41463e-05; //+- 1.92775e-06
const float RochCor2012::d0par = -0.000148871; //+- 3.16301e-06
const float RochCor2012::e0par = 1.59501; //+- 0.0249021
const float RochCor2012::d1par = 7.95495e-05; //+- 1.12386e-05
const float RochCor2012::e1par = -0.364823; //+- 0.17896
const float RochCor2012::d2par = 0.000152032; //+- 5.68386e-06
const float RochCor2012::e2par = 0.410195; //+- 0.0431732

const float RochCor2012::netabin[9] = {-2.4,-2.1,-1.4,-0.7,0.0,0.7,1.4,2.1,2.4};


const float RochCor2012::dcor_bf[8][8]={{0.000107019,0.000110872,0.000104914,0.000088356,0.000098934,0.000097173,0.000076036,0.000038352},
				    {0.000027861,0.000082929,0.000079044,0.000090323,0.000096823,0.000077726,0.000089548,0.000094251},
				    {0.000096623,0.000080357,0.000116928,0.000064621,0.000089653,0.000104048,0.000075808,0.000104033},
				    {0.000056079,0.000074857,0.000091264,0.000094404,0.000089990,0.000076868,0.000063443,0.000058654},
				    {0.000020316,0.000039838,0.000085370,0.000099329,0.000089701,0.000080920,0.000061223,0.000036706},
				    {0.000079484,0.000058457,0.000115420,0.000113329,0.000085512,0.000123082,0.000063315,0.000046771},
				    {-0.000014849,0.000063752,0.000076973,0.000115789,0.000099257,0.000090677,0.000079501,0.000033056},
				    {0.000013455,0.000084340,0.000099403,0.000102317,0.000104791,0.000125368,0.000116473,0.000060607}};
const float RochCor2012::dcor_ma[8][8]={{0.000100487,0.000017322,-0.000130583,-0.000151975,-0.000160698,-0.000163992,0.000058854,0.000537338},
				    {0.000797318,0.000133807,-0.000012265,-0.000039851,-0.000014145,-0.000020469,-0.000038474,-0.000089369},
				    {0.000127640,0.000079630,0.000096678,0.000101597,0.000103066,0.000051955,-0.000052329,-0.000232106},
				    {-0.000822448,-0.000093050,0.000224437,0.000199975,0.000208382,0.000161882,0.000008895,-0.000055717},
				    {-0.000659872,-0.000009989,0.000204988,0.000187290,0.000147524,0.000176117,0.000146629,-0.000072046},
				    {-0.000232191,0.000003074,0.000063732,0.000058353,0.000024090,0.000055329,0.000196989,0.000152461},
				    {0.000150637,-0.000069255,-0.000120556,-0.000129774,-0.000181534,-0.000185724,-0.000138647,0.000415873},
				    {-0.000069055,-0.000059101,-0.000171624,-0.000183451,-0.000229660,-0.000318313,-0.000075301,0.000616515}};
const float RochCor2012::mcor_bf[8][8]={{0.000051096,0.000053354,0.000088067,0.000069625,0.000111574,0.000054434,0.000049471,0.000029711},
				    {-0.000011598,0.000000375,0.000060415,0.000029346,0.000033469,0.000016833,0.000015773,0.000072907},
				    {-0.000092219,0.000008603,0.000056443,0.000072303,0.000068745,0.000029366,0.000000943,0.000025621},
				    {-0.000046305,0.000026115,0.000040818,0.000041803,0.000029672,0.000028866,-0.000037959,0.000082943},
				    {-0.000008281,0.000014787,0.000049244,0.000048117,0.000030259,0.000017593,-0.000052198,0.000045353},
				    {0.000055014,0.000010201,0.000047900,0.000049045,0.000057969,0.000055183,0.000046076,0.000138166},
				    {-0.000013461,0.000002464,0.000078624,0.000083703,0.000082662,0.000059725,0.000023476,0.000066592},
				    {0.000015052,0.000014736,0.000068918,0.000076644,0.000089550,0.000082794,0.000028474,0.000108749}};
const float RochCor2012::mcor_ma[8][8]={{0.000385236,0.000085244,-0.000048248,-0.000046507,-0.000039849,-0.000074739,-0.000027738,0.000147033},
				    {0.000117999,-0.000031813,-0.000046653,-0.000014969,-0.000052586,-0.000031938,0.000064566,0.000151759},
				    {-0.000049906,-0.000028434,0.000014103,-0.000003392,0.000000074,0.000017360,0.000041907,0.000002294},
				    {0.000203101,0.000176454,0.000070179,0.000045394,0.000059526,0.000098039,0.000121830,0.000117598},
				    {0.000226645,0.000158641,0.000060774,0.000062295,0.000058320,0.000068351,0.000102766,0.000224776},
				    {-0.000416113,-0.000065963,0.000036343,0.000031159,0.000037239,0.000038536,0.000012286,0.000196513},
				    {-0.000465413,-0.000067328,-0.000001936,0.000020537,-0.000037656,-0.000004558,0.000001761,-0.000050296},
				    {0.000273417,0.000001643,-0.000056352,-0.000055155,-0.000073923,-0.000052147,-0.000106017,-0.000373331}};


const float RochCor2012::dcor_bfer[8][8]={{0.000060607,0.000033327,0.000028317,0.000026616,0.000026635,0.000028568,0.000033216,0.000059162},
				      {0.000058346,0.000032683,0.000028301,0.000026814,0.000026649,0.000028376,0.000032958,0.000058995},
				      {0.000059405,0.000033098,0.000028297,0.000026759,0.000026544,0.000028262,0.000032849,0.000059237},
				      {0.000058096,0.000032952,0.000028116,0.000026624,0.000026511,0.000028165,0.000032528,0.000058668},
				      {0.000059777,0.000033148,0.000028356,0.000026645,0.000026687,0.000028264,0.000033157,0.000059752},
				      {0.000058886,0.000032805,0.000028226,0.000026821,0.000026998,0.000028379,0.000032977,0.000057995},
				      {0.000058887,0.000033128,0.000028232,0.000026494,0.000026707,0.000028133,0.000032632,0.000058308},
				      {0.000058322,0.000033103,0.000028333,0.000026601,0.000026502,0.000028277,0.000032777,0.000058169}};
const float RochCor2012::dcor_maer[8][8]={{0.000060607,0.000033327,0.000028317,0.000026616,0.000026635,0.000028568,0.000033216,0.000059162},
				      {0.000058346,0.000032683,0.000028301,0.000026814,0.000026649,0.000028376,0.000032958,0.000058995},
				      {0.000059405,0.000033098,0.000028297,0.000026759,0.000026544,0.000028262,0.000032849,0.000059237},
				      {0.000058096,0.000032952,0.000028116,0.000026624,0.000026511,0.000028165,0.000032528,0.000058668},
				      {0.000059777,0.000033148,0.000028356,0.000026645,0.000026687,0.000028264,0.000033157,0.000059752},
				      {0.000058886,0.000032805,0.000028226,0.000026821,0.000026998,0.000028379,0.000032977,0.000057995},
				      {0.000058887,0.000033128,0.000028232,0.000026494,0.000026707,0.000028133,0.000032632,0.000058308},
				      {0.000058322,0.000033103,0.000028333,0.000026601,0.000026502,0.000028277,0.000032777,0.000058169}};
const float RochCor2012::mcor_bfer[8][8]={{0.000057442,0.000032512,0.000028476,0.000026441,0.000026600,0.000028609,0.000033152,0.000058520},
				      {0.000059126,0.000033023,0.000028420,0.000026487,0.000026433,0.000028348,0.000033027,0.000058155},
				      {0.000058834,0.000033739,0.000028770,0.000026521,0.000026452,0.000028520,0.000033074,0.000058708},
				      {0.000057004,0.000033420,0.000028307,0.000026526,0.000026567,0.000028409,0.000032725,0.000059488},
				      {0.000059638,0.000032488,0.000028209,0.000026119,0.000026508,0.000028089,0.000033404,0.000059850},
				      {0.000058238,0.000033400,0.000028385,0.000026755,0.000026870,0.000028312,0.000032999,0.000056991},
				      {0.000059148,0.000034134,0.000028568,0.000026543,0.000026583,0.000027788,0.000032530,0.000059904},
				      {0.000060323,0.000033258,0.000028902,0.000026263,0.000026405,0.000028263,0.000032459,0.000059398}};
const float RochCor2012::mcor_maer[8][8]={{0.000057442,0.000032512,0.000028476,0.000026441,0.000026600,0.000028609,0.000033152,0.000058520},
				      {0.000059126,0.000033023,0.000028420,0.000026487,0.000026433,0.000028348,0.000033027,0.000058155},
				      {0.000058834,0.000033739,0.000028770,0.000026521,0.000026452,0.000028520,0.000033074,0.000058708},
				      {0.000057004,0.000033420,0.000028307,0.000026526,0.000026567,0.000028409,0.000032725,0.000059488},
				      {0.000059638,0.000032488,0.000028209,0.000026119,0.000026508,0.000028089,0.000033404,0.000059850},
				      {0.000058238,0.000033400,0.000028385,0.000026755,0.000026870,0.000028312,0.000032999,0.000056991},
				      {0.000059148,0.000034134,0.000028568,0.000026543,0.000026583,0.000027788,0.000032530,0.000059904},
				      {0.000060323,0.000033258,0.000028902,0.000026263,0.000026405,0.000028263,0.000032459,0.000059398}};

//=======================================================================================================

const float RochCor2012::dmavg[8][8]={{0.025806983,0.025159891,0.024975842,0.025455723,0.025425207,0.024926903,0.025207309,0.026048885},
				  {0.025750965,0.025097309,0.024989121,0.025452482,0.025527396,0.024992650,0.025032483,0.025745209},
				  {0.025804636,0.025158395,0.025016371,0.025488043,0.025484602,0.025030160,0.025056485,0.025666000},
				  {0.025856440,0.025213197,0.025026234,0.025434566,0.025508311,0.025036940,0.025187392,0.025890424},
				  {0.025876249,0.025281013,0.024931971,0.025371491,0.025499941,0.025069176,0.025232100,0.025805794},
				  {0.025650327,0.025166171,0.025005627,0.025408096,0.025460445,0.025019452,0.025223815,0.025741512},
				  {0.025876350,0.025175626,0.024978362,0.025447193,0.025421566,0.024991240,0.025167436,0.025891024},
				  {0.025871826,0.025199998,0.024895251,0.025438601,0.025413080,0.024985093,0.025161586,0.025951909}};
const float RochCor2012::dpavg[8][8]={{0.025863485,0.025162835,0.025030071,0.025502551,0.025451110,0.024937023,0.025201191,0.025910362},
				  {0.025625774,0.024974553,0.025011365,0.025516130,0.025527164,0.024980470,0.025081821,0.025774996},
				  {0.025843368,0.025180360,0.024976660,0.025487988,0.025465345,0.025063817,0.025222588,0.025848188},
				  {0.026043606,0.025279730,0.024916318,0.025446864,0.025459185,0.024968978,0.025218397,0.026001596},
				  {0.026105632,0.025315669,0.024938580,0.025393917,0.025404370,0.024991038,0.025172806,0.025925590},
				  {0.025834291,0.025225312,0.024968093,0.025365828,0.025470435,0.024970995,0.025149210,0.025854499},
				  {0.025907967,0.025225916,0.025005253,0.025466169,0.025463152,0.024998250,0.025108125,0.025770905},
				  {0.025902320,0.025180646,0.024998309,0.025482915,0.025518436,0.025018455,0.025214326,0.025771985}};
const float RochCor2012::mmavg[8][8]={{0.026128584,0.025420381,0.025425881,0.025276864,0.025621923,0.024943013,0.024997684,0.024951389},
				  {0.025892397,0.025324140,0.024767618,0.025466916,0.025547697,0.024931444,0.024865503,0.024965166},
				  {0.025884594,0.025299910,0.024911579,0.025545151,0.025692597,0.025396692,0.025361045,0.025631407},
				  {0.026016195,0.024945150,0.024826044,0.025509113,0.025358243,0.025056858,0.025043673,0.025082216},
				  {0.025525420,0.025148992,0.025057081,0.025180734,0.025324788,0.024978238,0.025439115,0.025149794},
				  {0.025679313,0.024968559,0.024960063,0.025743545,0.025488564,0.024816573,0.025175969,0.025249409},
				  {0.025616512,0.024970570,0.024939109,0.025640487,0.025595783,0.024913486,0.025093671,0.025319956},
				  {0.025649297,0.025295899,0.025035573,0.025319031,0.025256918,0.025451105,0.024914639,0.025190688}};
const float RochCor2012::mpavg[8][8]={{0.025604088,0.025265064,0.024723310,0.025622864,0.025484807,0.025519480,0.025055743,0.024885352},
				  {0.025536602,0.024972232,0.024960433,0.025341841,0.025646744,0.024954528,0.025367321,0.024594933},
				  {0.025733960,0.024972564,0.025092482,0.025435571,0.025381427,0.024627063,0.025213528,0.026130881},
				  {0.025760733,0.025261410,0.024756177,0.025349453,0.025471219,0.024913743,0.025367404,0.025344302},
				  {0.025564542,0.025115323,0.025066901,0.025519108,0.025278455,0.025177200,0.025202904,0.025660497},
				  {0.026146946,0.025127524,0.025185574,0.025478394,0.025289463,0.025376077,0.025378334,0.025742637},
				  {0.026086491,0.025343798,0.025090656,0.025447269,0.025322028,0.024798416,0.024731958,0.024964551},
				  {0.025714620,0.024896867,0.024996409,0.025739719,0.025418277,0.024759519,0.025287654,0.024960841}};

//=======================================================================================================

//===============================================================================================
//parameters for Z pt correction

const float RochCor2012::ptlow[85] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5,
				  6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5,
				  10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5,
				  15.0, 15.5, 16.0, 16.5, 17.0, 17.5, 18.0, 18.5, 19.0, 19.5,
				  20.0, 20.5, 21.0, 21.5, 22.0, 22.5, 23.0, 23.5, 24.0, 24.5,
				  25.0, 26.0, 27.0, 28.0, 29.0,
				  30.0, 32.0, 34.0, 36.0, 38.0,
				  40.0, 44.0, 48.0, 52.0, 56.0,
				  60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0,
				  100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 175.0,
				  200.0, 250.0, 350.0, 500.0, 1000.0};

//int nptbins( sizeof(ptlow)/sizeof(float) - 1 );
  
const float RochCor2012::zptscl[84] = {1.49177,1.45654,1.36283,1.28569,1.2418,1.12336,1.10416,1.08731,0.994051,0.96532,
				   0.94427,0.932725,0.918082,0.899665,0.898398,0.927687,0.908047,0.892392,0.924027,0.945895,
				   0.937149,0.923983,0.923387,0.955362,0.947812,0.962943,0.948781,0.961555,0.95222,0.999207,
				   0.973884,0.993013,0.953487,0.951402,0.985583,0.986603,0.981388,1.00022,1.0294,0.964748,
				   0.974592,1.01546,0.992343,1.00101,0.990866,0.98982,1.02924,1.02265,0.967695,1.02411,
				   0.97331,1.01052,1.01561,0.992594,0.976504,1.01205,0.981111,1.00078,1.02078,1.00719,
				   1.0099,1.02865,1.03845,1.03254,1.09815,1.10263,1.06302,1.0725,1.14703,1.10574,
				   1.13911,1.16947,1.1709,1.11413,1.28793,1.18953,1.20212,1.18112,1.25471,1.15329,
				   1.14276,1.17223,1.09173,2.00229};
  
const float RochCor2012::zptscler[84] = {0.0270027,0.0154334,0.0115338,0.00958085,0.0084683,0.00736665,0.0069567,0.00671434,
				     0.00617693,0.00601943,0.00594735,0.00594569,0.00594903,0.00595495,0.00608115,0.00633704,
				     0.0063916,0.0064468,0.00678106,0.00706769,0.00717517,0.00727958,0.00747182,0.00785544,
				     0.00798754,0.00828787,0.00839147,0.00865826,0.00876775,0.00933276,0.00935768,0.0097289,
				     0.00962058,0.00983828,0.0103044,0.0104871,0.0106575,0.0110388,0.0114986,0.0111494,
				     0.0115202,0.0121059,0.0121345,0.0124923,0.0125972,0.0128401,0.0134519,0.0136279,
				     0.0133414,0.014186,0.00992195,0.0105984,0.0109484,0.0111756,0.0114579,0.00870013,
				     0.00904749,0.00970734,0.0104583,0.0109818,0.00837852,0.00939894,0.010415,0.0113433,
				     0.013007,0.0128788,0.0140174,0.0156993,0.0181717,0.019765,0.0222326,0.0249408,
				     0.0272806,0.0211706,0.0278087,0.0306654,0.0361387,0.041327,0.0341513,0.0440116,
				     0.0473006,0.0680212,0.149162,0.56279};

RochCor2012::~RochCor2012(){
}

RochCor2012::RochCor2012(){
  
  eran.SetSeed(123456);
  sran.SetSeed(1234);
  
  gscler_mc_dev=0;
  gscler_da_dev=0;

  for(int i=0; i<8; ++i){
      for(int j=0; j<8; ++j){
          mptsys_mc_dm[i][j]=0;
          mptsys_mc_da[i][j]=0;
          mptsys_da_dm[i][j]=0;
          mptsys_da_da[i][j]=0;
      }
  }

}

RochCor2012::RochCor2012(int seed){
  eran.SetSeed(123456);
  sran.SetSeed(seed);

  gscler_mc_dev=sran.Gaus(0.0, 1.0);
  gscler_da_dev=sran.Gaus(0.0, 1.0);

  for(int i=0; i<8; ++i){
      for(int j=0; j<8; ++j){
          mptsys_mc_dm[i][j]=sran.Gaus(0.0, 1.0);
          mptsys_mc_da[i][j]=sran.Gaus(0.0, 1.0);
          mptsys_da_dm[i][j]=sran.Gaus(0.0, 1.0);
          mptsys_da_da[i][j]=sran.Gaus(0.0, 1.0);
      }
  }
}

void RochCor2012::momcor_mc( TLorentzVector& mu, float charge, float sysdev, int runopt){
  
  //sysdev == num : deviation = num

  float ptmu = mu.Pt();
  float muphi = mu.Phi();
  float mueta = mu.Eta(); // same with mu.Eta() in Root

  float px = mu.Px();
  float py = mu.Py();
  float pz = mu.Pz();
  float e = mu.E();

  int mu_phibin = phibin(muphi);
  int mu_etabin = etabin(mueta);

  //float mptsys = sran.Gaus(0.0,sysdev);
  
  float dm = (mcor_bf[mu_phibin][mu_etabin] + mptsys_mc_dm[mu_phibin][mu_etabin]*mcor_bfer[mu_phibin][mu_etabin])/mmavg[mu_phibin][mu_etabin];
  float da = mcor_ma[mu_phibin][mu_etabin] + mptsys_mc_da[mu_phibin][mu_etabin]*mcor_maer[mu_phibin][mu_etabin];
  
  float cor = 1.0/(1.0 + dm + charge*da*ptmu);

  //for the momentum tuning - eta,phi,Q correction
  px *= cor;
  py *= cor;
  pz *= cor;
  e  *= cor;
  
  float gscler = 0.0;
  float deltaer = 0.0;
  float sfer = 0.0;
  
  gscler = TMath::Sqrt( TMath::Power(mgscl_stat,2) + TMath::Power(mgscl_syst,2) );
  deltaer = TMath::Sqrt( TMath::Power(delta_stat,2) + TMath::Power(delta_syst,2) );
  sfer = TMath::Sqrt( TMath::Power(sf_stat,2) + TMath::Power(sf_syst,2) );
  
  float tune = 1.0/(1.0 + (delta + sysdev*deltaer)*sqrt(px*px + py*py)*eran.Gaus(1.0,(sf + sysdev*sfer)));
  
  px *= (tune); 
  py *= (tune);  
  pz *= (tune);  
  e  *= (tune);   
  
  float gscl = (genm_smr/mrecm);
  
  px *= (gscl + gscler_mc_dev*gscler);
  py *= (gscl + gscler_mc_dev*gscler);
  pz *= (gscl + gscler_mc_dev*gscler);
  e  *= (gscl + gscler_mc_dev*gscler);

  mu.SetPxPyPzE(px,py,pz,e);
  
}


void RochCor2012::momcor_data( TLorentzVector& mu, float charge, float sysdev, int runopt){
  
  float ptmu = mu.Pt();

  float muphi = mu.Phi();
  float mueta = mu.Eta(); // same with mu.Eta() in Root

  float px = mu.Px();
  float py = mu.Py();
  float pz = mu.Pz();
  float e = mu.E();
  
  int mu_phibin = phibin(muphi);
  int mu_etabin = etabin(mueta);
  
  //float mptsys1 = sran.Gaus(0.0,sysdev);
  
  float dm = (dcor_bf[mu_phibin][mu_etabin] + mptsys_da_dm[mu_phibin][mu_etabin]*dcor_bfer[mu_phibin][mu_etabin])/dmavg[mu_phibin][mu_etabin];
  float da = dcor_ma[mu_phibin][mu_etabin] + mptsys_da_da[mu_phibin][mu_etabin]*dcor_maer[mu_phibin][mu_etabin];
  
  float cor = 1.0/(1.0 + dm + charge*da*ptmu);
  
  px *= cor;
  py *= cor;
  pz *= cor;
  e  *= cor;
  
  //after Z pt correction
  float gscler = 0.0;
  
  gscler = TMath::Sqrt( TMath::Power(dgscl_stat,2) + TMath::Power(dgscl_syst,2) );
  
  float gscl = (genm_smr/drecm);
  
  px *= (gscl + gscler_da_dev*gscler);
  py *= (gscl + gscler_da_dev*gscler);
  pz *= (gscl + gscler_da_dev*gscler);
  e  *= (gscl + gscler_da_dev*gscler);
  
  mu.SetPxPyPzE(px,py,pz,e);
  
}

void RochCor2012::musclefit_data( TLorentzVector& mu, TLorentzVector& mubar){

  float dpar1 = 0.0;
  float dpar2 = 0.0;
  float epar1 = 0.0;
  float epar2 = 0.0;
  
  if(fabs(mu.PseudoRapidity())<=0.9){
    dpar1 = d0par;
    epar1 = e0par;
  }else if(mu.PseudoRapidity()>0.9){
    dpar1 = d1par;
    epar1 = e1par;
  }else if(mu.PseudoRapidity()<-0.9){
    dpar1 = d2par;
    epar1 = e2par;
  }

  if(fabs(mubar.PseudoRapidity())<=0.9){
    dpar2 = d0par;
    epar2 = e0par;
  }else if(mubar.PseudoRapidity()>0.9){
    dpar2 = d1par;
    epar2 = e1par;
  }else if(mubar.PseudoRapidity()<-0.9){
    dpar2 = d2par;
    epar2 = e2par;
  }

  float corr1 = 1.0 + bpar*mu.Pt() + (-1.0)*cpar*mu.Pt()*TMath::Sign(float(1.0),float(mu.PseudoRapidity()))*TMath::Power(mu.PseudoRapidity(),2)
    + (-1.0)*dpar1*mu.Pt()*sin(mu.Phi() + epar1);
  float corr2 = 1.0 + bpar*mubar.Pt() + (1.0)*cpar*mubar.Pt()*TMath::Sign(float(1.0),float(mubar.PseudoRapidity()))*TMath::Power(mubar.PseudoRapidity(),2)
    + (1.0)*dpar2*mubar.Pt()*sin(mubar.Phi() + epar2);
  
  float px1 = mu.Px();
  float py1 = mu.Py();
  float pz1 = mu.Pz();
  float e1 = mu.E();
  
  float px2 = mubar.Px();
  float py2 = mubar.Py();
  float pz2 = mubar.Pz();
  float e2 = mubar.E();

  px1 *= corr1;
  py1 *= corr1;
  pz1 *= corr1;
  e1 *= corr1;
  
  px2 *= corr2;
  py2 *= corr2;
  pz2 *= corr2;
  e2 *= corr2;
  
  mu.SetPxPyPzE(px1,py1,pz1,e1);
  mubar.SetPxPyPzE(px2,py2,pz2,e2);

}

Int_t RochCor2012::phibin(float phi){
  
  int nphibin = -1;
  
  for(int i=0; i<8; i++){
    if(-pi+(2.0*pi/8.0)*i <= phi && -pi+(2.0*pi/8.0)*(i+1) > phi){
      nphibin = i;
      break;
    }
  }
  
  return nphibin;
}

Int_t RochCor2012::etabin(float eta){

  int nbin = -1;
  
  for(int i=0; i<8; i++){
    if(netabin[i] <= eta && netabin[i+1] > eta){
      nbin = i;
      break;
    }
  }
  
  return nbin;
}

float RochCor2012::zptcor(float gzpt) {
  int ibin( 0 );
  
  // mcptscl[] = 84 bins: [0] and [83] are the underflow and overflow
  if ( gzpt > ptlow[nptbins] ) return nptbins-1;
  if ( gzpt < ptlow[0      ] ) return 0;
  
  for ( int i=0; i<nptbins; ++i ) {
    if ( gzpt>=ptlow[i] && gzpt<ptlow[i+1] ) { ibin=i; break; }
  }

  float zptwt = zptscl[ibin];

  return zptwt;
}
