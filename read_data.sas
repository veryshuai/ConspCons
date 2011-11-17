
* CEX Project Main File
* Largely based on NBER CEX extract example program

options nocenter ls=132 ps=84 nonumber;
libname disk '/gpfs/work/dcj138/CEX_code';
run;

/*  Creates quarterly sas datasets and presents quarterly tables
    for 1990 to 1994 Q1.  SEE CBOFAM.OLD for merge of quarterly data for
    annualized tables for NIPA comparison only.
    Note that RENTS are income and RENT is expenditure
*/

%macro one(q);
  filename fam pipe "zcat /gpfs/work/dcj138/NBER_CEX/ffile&q..Z";
  data disk.fam&q;
  infile fam lrecl=1478;
  input
     newid         1-7
     cutenur        10
     repstat        13
     srepstat       14
     totwt       19-29
     adjwt       30-40
     fullyr         41
     numearn     42-45
     numauto     46-49
     vehq        50-53
     famsize     54-57
     membcnt     58-59
                                /* position  */
                                /* X = not included in any summary variable */
/* personal income */
     wages          60 -   69   /*   1      */
     bus            70 -   79   /*   2      */
     farm           80 -   89   /*   3      */
     rents          90 -   99   /*   4      */
     div           100 -  109   /*   5      */
     int           110 -  119   /*   6      */
     pension       120 -  129   /*   7      */
     socsec        130 -  139   /*   8      */
     ssi           140 -  149   /*   9      */
     unemp         150 -  159   /*  10      */
     workcomp      160 -  169   /*  11      */
     welfare       170 -  179   /*  12      */
     scholar       180 -  189   /*  13      */
     foodstmp      190 -  199   /*  14     */

/* pension contributions */
     gvpremia      200 -  209   /*  15     */
     rrpremia      210 -  219   /*  16     */
     sspremia      220 -  229   /*  17     */

/* tax and nontax payments */
     fedtax        230 -  239   /*  18     */
     statax        240 -  249   /*  19     */
     pproptax      250 -  259   /*  20     */
     othtax        260 -  269   /*  21     */
     nontax        270 -  279   /*  22     */

/* food and tobacco */
     foodhome      280 -  289   /*  23     */
     foodout       290 -  299   /*  24     */
     foodwork      300 -  309   /*  25     */
     tobacco       310 -  319   /*  26     */
     alcohol       320 -  329   /*  27     */
     niteclub      330 -  339   /*  28     */

/* clothin, accessories, jewely  */
     clothes       340 -  349   /*  29     */
     tailors       350 -  359   /*  30     */
     jewelry       360 -  369   /*  31     */

/* personal care  */
     toiletry      370 -  379   /*  32     */
     hlthbeau      380 -  389   /*  33     */

/* housing  */
     renthome      390 -  399   /*  34     */
     rentothr      400 -  409   /*  35     */

/* household operation  */
     furnish       410 -  419   /*  36     */
     housuppl      420 -  429   /*  37     */
     elect         430 -  439   /*  38     */
     gas           440 -  449   /*  39     */
     water         450 -  459   /*  40     */
     homefuel      460 -  469   /*  41     */
     telephon      470 -  479   /*  42     */
     servants      480 -  489   /*  43     */

/* medical care  */
     drugs         490 -  499   /*  44     */
     orthopd       500 -  509   /*  45     */
     doctors       510 -  519   /*  46     */
     hospital      520 -  529   /*  47     */
     nurshome      530 -  539   /*  48     */
     helthins      540 -  549   /*  49     */

     busiserv      550 -  559   /*  50     */
     lifeins       560 -  569   /*  51     */

     autos         570 -  579   /*  52     */
     parts         580 -  589   /*  53     */
     carservs      590 -  599   /*  54     */
     gasoline      600 -  609   /*  55     */
     tolls         610 -  619   /*  56     */
     autoins       620 -  629   /*  57     */
     masstran      630 -  639   /*  58     */
     othtrans      640 -  649   /*  59     */
     airfare       650 -  659   /*  60     */

     books         660 -  669   /*  61     */
     pubs          670 -  679   /*  62     */
     recsport      680 -  689   /*  63     */
     othrec        690 -  699   /*  64     */
     gambling      700 -  709   /*  65     */

     highedu       710 -  719   /*  66     */
     lowedu        720 -  729   /*  67     */
     othedu        730 -  739   /*  68     */

     charity       740 -  749   /*  69     */

     insrefnd      750 -  759   /*  70     */
     intauto       760 -  769   /*  71     */
     intoth        770 -  779   /*  72     */
     rentnpay      780 -  789   /*  73     */
     homevalu      790 -  799   /*  74  X  */
     homeval2      800 -  809   /*  75  X  */
     ohint         810 -  819   /*  76     */
     ohtax         820 -  829   /*  77     */
     ohmaint       830 -  839   /*  78     */
     ohprinc       840 -  849   /*  79     */
     ohlump        850 -  859   /*  80     */
     ohsold        860 -  869   /*  81     */
     ohbuy         870 -  879   /*  82     */
     ohmort1       880 -  889   /*  83     */
     ohmort2       890 -  899   /*  84     */
     ohtrans       900 -  909   /*  85     */
     housadd       910 -  919   /*  86     */
     checkng       920 -  929   /*  87  X  */
     dcheckng      930 -  939   /*  88     */
     saving        940 -  949   /*  89  X  */
     dsaving       950 -  959   /*  90     */
     securi        960 -  969   /*  91  X  */
     dsecuri       970 -  979   /*  92     */
     carloan       980 -  989   /*  93     */
     tradein       990 -  999   /*  94     */
     carsold      1000 - 1009   /*  95     */
     carprinc     1010 - 1019   /*  96  X  */
     investb      1020 - 1029   /*  97     */
     pnpremia     1030 - 1039   /*  98     */
     sepremia     1040 - 1049   /*  99     */
     deltaiou     1050 - 1059   /* 100     */
     owe1         1060 - 1069   /* 101  X  */
     owe1q1       1070 - 1079   /* 102  X  */
     owe5         1080 - 1089   /* 103  X  */
     owe5q1       1090 - 1099   /* 104  X  */
     give         1100 - 1109   /* 105     */
     receive      1110 - 1119   /* 106     */
     give2        1120 - 1129   /* 107  X  */
     lumpsums     1130 - 1139   /* 108     */
     othernet     1140 - 1149   /* 109     */
;

/*  additional variables computed from the CBO data  */

families = 1;
   dummy = 1;
 quarter = &q;
 income = wages + bus + farm + rents + div + int + pension + socsec + ssi
           + unemp + workcomp + welfare + foodstmp + scholar + rentnpay;
  premia = gvpremia + rrpremia + sspremia;
 taxplus = fedtax + statax + pproptax + othtax + nontax;
  foodin = foodhome+foodwork;
    food = foodin+foodout;
 nonfood = tobacco+alcohol+niteclub;
clothing = clothes+tailors+jewelry;
perscare = toiletry + hlthbeau;
    rent = renthome + rentothr + rentnpay;
 ownoope = ohint+ohtax+ohmaint+ohprinc+ohlump+ohtrans+housadd;
ownoccup = 12*homeval2;
 housing = rent+ownoccup;
 houseop = furnish+housuppl+elect+gas+water+homefuel+telephon+servants;
 medical = drugs+orthopd+doctors+hospital+nurshome+helthins;
 persbus = busiserv+lifeins;
  transp = autos+parts+carservs+gasoline+tolls+autoins+masstran+
           othtrans+airfare;
recreate = books+pubs+recsport+othrec+gambling;
    educ = highedu+lowedu+othedu;
exptotal = food+nonfood+clothing+perscare+housing+houseop+persbus+
           transp+recreate+educ+charity;

title "Family File Weighted by Adjusted Weight Period&q";
proc means mean sum data=disk.fam&q; weight adjwt;
run;
%mend;


%one(801); %one(802); %one(803); %one(804);
%one(811);%one(812);%one(813);%one(814);
%one(821);%one(822);%one(823);%one(824);
%one(831);%one(832);%one(833);%one(834);
%one(841);%one(842);%one(843);%one(844);
%one(851); %one(852);
%one(861); %one(862); %one(863); %one(864);
%one(871);%one(872);%one(873);%one(874);
%one(881); %one(882); %one(883); %one(884);
%one(891); %one(892); %one(893); %one(894);
%one(901);%one(902);%one(903);%one(904);
%one(911);%one(912);%one(913);%one(914);
%one(921);%one(922);%one(923);%one(924);
%one(931);%one(932);%one(933);%one(934);
%one(941);%one(942); %one(943);%one(944);
%one(951);%one(952);
%one(961);%one(962);%one(963);%one(964);
%one(971);%one(972);%one(973);%one(974);
%one(981);%one(982); %one(983); %one(984);
%one(991); %one(992); %one(993); %one(994);
%one(001);%one(002);%one(003);%one(004);
%one(011);%one(012);%one(013);%one(014);
%one(021);%one(022);%one(023);%one(024);
%one(031);%one(032);
run;

%macro one(q);
  filename memb "/gpfs/work/dcj138/NBER_CEX/mfile&q";
  data disk.memb&q;
  infile memb lrecl=272;
  input
  newid          1-7
  age            8-10
  relation $     11
  educatio $     12-13
  race     $     14
  sex      $     15
  weeksin        16-17
  emplcont $     18
  incoll   $     19
  nonwork  $     20
  marital  $     21
  empstat  $     22
  emptype  $     23
  hrswkd         24-26
  wkswkd         27-28
  occup    $     29-30
  indust   $     31-32
  fedtax         33-42
  gvpremia       43-52
  pripemia       53-62
  rrpremia       63-72
  statax         73-82
  farminc        83-92
  irakeogh       93-102
  fica          103-112
  businc        113-122
  wage          123-132
  socbens       133-142
  ssi           143-152
;
quarter = "&q";
members = 1;

/*  tables  */

   title "member file unweighted for period&q";
   proc means n sum mean;
   run;

%mend;

%one(801); %one(802); %one(803); %one(804);
%one(811);%one(812);%one(813);%one(814);
%one(821);%one(822);%one(823);%one(824);
%one(831);%one(832);%one(833);%one(834);
%one(841);%one(842);%one(843);%one(844);
%one(851); %one(852);
%one(861); %one(862); %one(863); %one(864);
%one(871);%one(872);%one(873);%one(874);
%one(881); %one(882); %one(883); %one(884);
%one(891); %one(892); %one(893); %one(894);
%one(901);%one(902);%one(903);%one(904);
%one(911);%one(912);%one(913);%one(914);
%one(921);%one(922);%one(923);%one(924);
%one(931);%one(932);%one(933);%one(934);
%one(941);%one(942); %one(943);%one(944);
%one(951);%one(952);
%one(961);%one(962);%one(963);%one(964);
%one(971);%one(972);%one(973);%one(974);
%one(981);%one(982); %one(983); %one(984);
%one(991); %one(992); %one(993); %one(994);
%one(001);%one(002);%one(003);%one(004);
%one(011);%one(012);%one(013);%one(014);
%one(021);%one(022);%one(023);%one(024);
%one(031);%one(032);
run;
