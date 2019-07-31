#include "dncdm_decay_kernel.h"

double dncdm_decay_kernel(double x, int l){
  if (l == 0)
    return 1;
  if (l == 1)
    return x;
  if (l == 2)
    return (x*(-3 + 5*pow(x,2)) + 3*pow(-1 + pow(x,2),2)*atanh(x))/(2.*pow(x,3));
  if (l == 3)
    return (-15*x + 25*pow(x,3) - 8*pow(x,5) + 15*pow(-1 + pow(x,2),2)*atanh(x))/(2.*pow(x,4));
  if (l == 4)
    return -(105*x - 190*pow(x,3) + 81*pow(x,5) + 15*(-7 + pow(x,2))*pow(-1 + pow(x,2),2)*atanh(x))/(4.*pow(x,5));
  if (l == 5)
    return (x*(-315 + 630*pow(x,2) - 343*pow(x,4) + 32*pow(x,6)) - 105*(-3 + pow(x,2))*pow(-1 + pow(x,2),2)*atanh(x))/(4.*pow(x,6));
  if (l == 6)
    return (2*x*(-3465 + 7665*pow(x,2) - 5103*pow(x,4) + 919*pow(x,6)) - 105*pow(-1 + pow(x,2),2)*(33 - 18*pow(x,2) + pow(x,4))*log(1 - x) + 105*pow(-1 + pow(x,2),2)*(33 - 18*pow(x,2) + pow(x,4))*log(1 + x))/(32.*pow(x,7));
  if (l == 7)
    return (-45045*x + 109725*pow(x,3) - 86499*pow(x,5) + 22923*pow(x,7) - 1024*pow(x,9) + 315*pow(-1 + pow(x,2),2)*(143 - 110*pow(x,2) + 15*pow(x,4))*atanh(x))/(80.*pow(x,8));
  if (l == 8)
    return -(45045*x - 120120*pow(x,3) + 109494*pow(x,5) - 38232*pow(x,7) + 3781*pow(x,9) + 315*pow(-1 + pow(x,2),2)*(-143 + 143*pow(x,2) - 33*pow(x,4) + pow(x,6))*atanh(x))/(32.*pow(x,9));
  if (l == 9)
    return (x*(-765765 + 2222220*pow(x,2) - 2300298*pow(x,4) + 995940*pow(x,6) - 155969*pow(x,8) + 4096*pow(x,10)) - 3465*pow(-1 + pow(x,2),2)*(-221 + 273*pow(x,2) - 91*pow(x,4) + 7*pow(x,6))*atanh(x))/(224.*pow(x,10));
  if (l == 10)
    return (x*(-14549535 + 45690645*pow(x,2) - 52954902*pow(x,4) + 27353898*pow(x,6) - 5907275*pow(x,8) + 368961*pow(x,10)) + 3465*pow(-1 + pow(x,2),2)*(4199 - 6188*pow(x,2) + 2730*pow(x,4) - 364*pow(x,6) + 7*pow(x,8))*atanh(x))/(1792.*pow(x,11));
  if (l == 11)
    return (-101846745*x + 344338995*pow(x,3) - 441795354*pow(x,5) + 265086822*pow(x,7) - 73221005*pow(x,9) + 7573735*pow(x,11) - 131072*pow(x,13) + 45045*pow(-1 + pow(x,2),2)*(2261 - 3876*pow(x,2) + 2142*pow(x,4) - 420*pow(x,6) + 21*pow(x,8))*atanh(x))/(5376.*pow(x,12));
  if (l == 12)
    return -(334639305*x - 1212461250*pow(x,3) + 1706175471*pow(x,5) - 1166034012*pow(x,7) + 392983591*pow(x,9) - 57797090*pow(x,11) + 2486305*pow(x,13) + 45045*pow(-1 + pow(x,2),2)*(-7429 + 14535*pow(x,2) - 9690*pow(x,4) + 2550*pow(x,6) - 225*pow(x,8) + 3*pow(x,10))*atanh(x))/(7680.*pow(x,13));
  if (l == 13)
    return (x*(-1673196525 + 6469693230*pow(x,2) - 9908233335*pow(x,4) + 7597351476*pow(x,6) - 3025265243*pow(x,8) + 581869470*pow(x,10) - 42726465*pow(x,12) + 524288*pow(x,14)) - 45045*pow(-1 + pow(x,2),2)*(-37145 + 81719*pow(x,2) - 63954*pow(x,4) + 21318*pow(x,6) - 2805*pow(x,8) + 99*pow(x,10))*atanh(x))/(16896.*pow(x,14));
  if (l == 14)
    return (x*(-15058768725 + 61908271425*pow(x,2) - 102511173765*pow(x,4) + 87144093465*pow(x,6) - 40053761319*pow(x,8) + 9559644731*pow(x,10) - 1020306175*pow(x,12) + 32067947*pow(x,14)) + 45045*pow(-1 + pow(x,2),2)*(334305 - 817190*pow(x,2) + 735471*pow(x,4) - 298452*pow(x,6) + 53295*pow(x,8) - 3366*pow(x,10) + 33*pow(x,12))*atanh(x))/(67584.*pow(x,15));
  if (l == 15)
    return (-436704293025*x + 1902424448925*pow(x,3) - 3386884405905*pow(x,5) + 3161305643925*pow(x,7) - 1647789572715*pow(x,9) + 470846672751*pow(x,11) - 66828525155*pow(x,13) + 3664464223*pow(x,15) - 33554432*pow(x,17) + 765765*pow(-1 + pow(x,2),2)*(570285 - 1533870*pow(x,2) + 1562275*pow(x,4) - 749892*pow(x,6) + 171171*pow(x,8) - 16302*pow(x,10) + 429*pow(x,12))*atanh(x))/(878592.*pow(x,16));
  if (l == 16)
    return -(4512611027925*x - 20767715268300*pow(x,3) + 39558381522660*pow(x,5) - 40218481454580*pow(x,7) + 23446773803310*pow(x,9) - 7816042702260*pow(x,11) + 1394788588452*pow(x,13) - 113028685196*pow(x,15) + 2709067893*pow(x,17) + 765765*pow(-1 + pow(x,2),2)*(-5892945 + 17298645*pow(x,2) + 143*pow(x,4)*(-137655 + 76475*pow(x,2) - 21413*pow(x,4) + 2793*pow(x,6) - 133*pow(x,8) + pow(x,10)))*atanh(x))/(4.100096e6*pow(x,17));
  if (l == 17)
    return (x*(-49638721307175 + 240672588156000*pow(x,2) + pow(x,4)*(-488410081319160 + 537236059365720*pow(x,2) - 346431069866310*pow(x,4) + 132114577760880*pow(x,6) + 7*pow(x,8)*(-4075197330000 + 444930605928*pow(x,2) - 18914963581*pow(x,4) + 134217728*pow(x,6)))) - 14549535*pow(-1 + pow(x,2),2)*(-3411705 + 10855425*pow(x,2) - 13656825*pow(x,4) + 143*pow(x,6)*(60375 - 20125*pow(x,2) + 3381*pow(x,4) - 245*pow(x,6) + 5*pow(x,8)))*atanh(x))/(2.050048e7*pow(x,18));
  if (l == 18)
    return (x*(-347471049150225 + 1770447726622575*pow(x,2) - 3813457159331820*pow(x,4) + 4512261664490580*pow(x,6) - 3189168636701190*pow(x,8) + 1370290746589050*pow(x,10) - 348189480512940*pow(x,12) + 48302582871060*pow(x,14) - 3074378234193*pow(x,16) + 58048958639*pow(x,18)) + 14549535*pow(-1 + pow(x,2),2)*(23881935 - 81880920*pow(x,2) + 112896420*pow(x,4) - 80120040*pow(x,6) + 31081050*pow(x,8) - 6446440*pow(x,10) + 644644*pow(x,12) - 24024*pow(x,14) + 143*pow(x,16))*atanh(x))/(6.5601536e7*pow(x,19));
  if (l == 19)
    return (-1836632688365475*x + 9811920578384925*pow(x,3) - 22357280076751620*pow(x,5) + 28315989541511100*pow(x,7) - 21768906765335730*pow(x,9) + 10410335156462670*pow(x,11) - 3048334822239300*pow(x,13) + 516126281438460*pow(x,15) - 44725843731555*pow(x,17) + 1517387879133*pow(x,19) - 8589934592*pow(x,21) + 14549535*pow(-1 + pow(x,2),2)*(126233085 - 463991880*pow(x,2) + 695987820*pow(x,4) - 548354040*pow(x,6) + 243221550*pow(x,8) - 60386040*pow(x,10) + 7827820*pow(x,12) - 447304*pow(x,14) + 7293*pow(x,16))*atanh(x))/(1.59318016e8*pow(x,20));
  if (l == 20)
    return (-2*x*(723521968143975 - 4044302283471450*pow(x,2) + 9719462192434995*pow(x,4) - 13119291992372280*pow(x,6) + 10900128371081950*pow(x,8) - 5743946177513180*pow(x,10) + 1906766526700750*pow(x,12) - 382692812952120*pow(x,14) + 42502420524795*pow(x,16) - 2181529466970*pow(x,18) + 33287922623*pow(x,20)) + 4849845*pow(-1 + pow(x,2),2)*(-149184555 + 585262485*pow(x,2) - 949074300*pow(x,4) + 822531060*pow(x,6) - 411265530*pow(x,8) + 119399670*pow(x,10) - 19213740*pow(x,12) + 1524900*pow(x,14) - 45747*pow(x,16) + 221*pow(x,18))*log(1 - x) - 4849845*pow(-1 + pow(x,2),2)*(-149184555 + 585262485*pow(x,2) - 949074300*pow(x,4) + 822531060*pow(x,6) - 411265530*pow(x,8) + 119399670*pow(x,10) - 19213740*pow(x,12) + 1524900*pow(x,14) - 45747*pow(x,16) + 221*pow(x,18))*log(1 + x))/(5.7933824e7*pow(x,21));
  if (l > 20)
  {
    printf("l is too large!!");
    exit(0);
  }
  return 0.;
}
