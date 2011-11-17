options nocenter ls=132 ps=82 nonumber;
libname disk '/gpfs/work/dcj138/CEX_code/';

/* This program creates merge91.ssd:
   1)  Cocatenates member files centering around 1991 and merges them by newid.
       Creates over65, under18 and marriage flags.  Only flags are kept.
   2)  Cocatenates family files centering around 1991.
   3)  Merges the member and family files
   4)  calculates 1991 taxes and netinc, income net of federal and state income
       taxes and social security payroll tax
   5)  DROPS OBSERVATIONS WITH ZERO CBO WEIGHTS OR ADJWT
   Note that RENTS is income and RENT is expenditure
*/


data temp;
set disk.memb801 disk.memb802 disk.memb803 disk.memb804
    disk.memb811 disk.memb812 disk.memb813 disk.memb814
    disk.memb821 disk.memb822 disk.memb823 disk.memb824
    disk.memb831 disk.memb832 disk.memb833 disk.memb834
    disk.memb841 disk.memb842 disk.memb843 disk.memb844
    disk.memb851 disk.memb852
    disk.memb861 disk.memb862 disk.memb863 disk.memb864
    disk.memb871 disk.memb872 disk.memb873 disk.memb874
    disk.memb881 disk.memb882 disk.memb883 disk.memb884
    disk.memb891 disk.memb892 disk.memb893 disk.memb894
    disk.memb901 disk.memb902 disk.memb903 disk.memb904
    disk.memb911 disk.memb912 disk.memb913 disk.memb914
    disk.memb931 disk.memb932 disk.memb933 disk.memb934
    disk.memb941 disk.memb942 disk.memb943 disk.memb944
    disk.memb951 disk.memb952
    disk.memb961 disk.memb962 disk.memb963 disk.memb964
    disk.memb971 disk.memb972 disk.memb973 disk.memb974
    disk.memb981 disk.memb982 disk.memb983 disk.memb984
    disk.memb991 disk.memb992 disk.memb993 disk.memb994
    disk.memb001 disk.memb002 disk.memb003 disk.memb004
    disk.memb011 disk.memb012 disk.memb013 disk.memb014
    disk.memb021 disk.memb022 disk.memb023 disk.memb024
    disk.memb031 disk.memb032;
proc sort;
by newid;
run;

data memb;
set temp;
by newid;
retain number50 over50 black married male year qt;
if first.newid then do;
    number50 = 0;
    over50 = 0;
    married = 0;
    black = 0;
    male = 0;
end;
year = substr(quarter,1,2);
qt = substr(quarter,3,1);
if age ge 50 then over50 = 1;
if age ge 50 then number50 = number50+1;
if marital = "1" then married = 1;
if race = "2" then black = 1;
if sex = "1" then male = 1;
if last.newid then output;
keep newid number50 over50 black married male year qt;
run;

proc sort;
by newid;
run;

data fam;
set disk.fam801 disk.fam802 disk.fam803 disk.fam804
    disk.fam811 disk.fam812 disk.fam813 disk.fam814
    disk.fam821 disk.fam822 disk.fam823 disk.fam824
    disk.fam831 disk.fam832 disk.fam833 disk.fam834
    disk.fam841 disk.fam842 disk.fam843 disk.fam844
    disk.fam851 disk.fam852
    disk.fam861 disk.fam862 disk.fam863 disk.fam864
    disk.fam871 disk.fam872 disk.fam873 disk.fam874
    disk.fam881 disk.fam882 disk.fam883 disk.fam884
    disk.fam891 disk.fam892 disk.fam893 disk.fam894
    disk.fam901 disk.fam902 disk.fam903 disk.fam904
    disk.fam911 disk.fam912 disk.fam913 disk.fam914
    disk.fam931 disk.fam932 disk.fam933 disk.fam934
    disk.fam941 disk.fam942 disk.fam943 disk.fam944
    disk.fam951 disk.fam952
    disk.fam961 disk.fam962 disk.fam963 disk.fam964
    disk.fam971 disk.fam972 disk.fam973 disk.fam974
    disk.fam981 disk.fam982 disk.fam983 disk.fam984
    disk.fam991 disk.fam992 disk.fam993 disk.fam994
    disk.fam001 disk.fam002 disk.fam003 disk.fam004
    disk.fam011 disk.fam012 disk.fam013 disk.fam014
    disk.fam021 disk.fam022 disk.fam023 disk.fam024
    disk.fam031 disk.fam032;
totwt = totwt/12;
adjwt = adjwt/12;
proc sort;
by newid;
run;

data disk.merged_cex;
merge fam memb;
by newid;
run;

