#include <iostream>
#include "myutils.h"
#include <sys/stat.h>
using namespace std;
#include <sstream>
#include <cstring>
#include "thirdpartprograms.h"

//================================================================================================================
// Some global variables

// name of pdf_viewer defined in thirdpartprograms.h 
string pdfreader=PDF_VIEWER; // you can change using the command set_pdfreader


//====================================================================================
//    ROUTINES TO DEAL WITH PDF


void set_pdfreader(string programname){ pdfreader = programname; }

// To see a pdf file. It uses the pdf reader defined by set_pdfreader.
int view_pdffile(string filename)
{
  char cmd[MAXCOMMANDSIZE];
  sprintf(cmd,"%s %s",pdfreader.c_str(),filename.c_str());
  system(cmd);
  return(0);
}

// Return true if variable x is fractional (within a certain small error).
bool IsFrac(double x) {
  double f;
  f = ceil(x)-x;
  if (f<MY_EPS) return(false);
  if (f>1.0-MY_EPS) return(false);
  return(true);  //    eps  <= ceil(x)-x <= 1-eps
}



//================================================================================================================
//    ROUTINES FOR TYPE CONVERSION

// convert a double value to string
std::string DoubleToString(double val)
{
    std::ostringstream out;
    out << val;
    return out.str();
}

// convert a int value to string
std::string IntToString(int val)
{
    std::ostringstream out;
    out << val;
    return out.str();
}

// convert a string to double
double StringToDouble(string s)
{
  double d;  stringstream(s) >> d;  return(d);
}

// convert a string to int
int StringToInt(string s)
{
  int n;  stringstream(s) >> n;  return(n);
}

//================================================================================================================
//    ROUTINES FOR COLORs


// Given a color code, return its name
string ColorName(int cor)
{
  switch (cor) {
  case WHITE: return("white");
  case BLACK: return("black");
  case RED: return("red");
  case GREEN: return("green");
  case BLUE: return("blue");
  case YELLOW: return("yellow");
  case MAGENTA: return("magenta");
  case CYAN: return("cyan");
  case GRAY: return("gray");
  case ORANGE: return("orange");
  }
  printf("ERROR: Unknown color number %d in routine GetColor.\n",cor);
  exit(1);
}


//================================================================================================================
//    ROUTINES FOR TIME MANIPULATION

 /* return time in seconds since 1/jan/1970*/
long time70(void)
{
  return((long) time(0));
}

// given a time in seconds, print corresponding time in days, hours, minutes and seconds
void printtime(long t)
{
  long dias,horas,minutos,segundos;
  if (t<1) {printf("<1 second");return;}
  minutos = (long) t/60;
  segundos = t - minutos*60;
  horas = (long) minutos/60;
  minutos = minutos-horas*60;
  dias = (long) horas/24;
  horas = horas - dias*24;
  if (dias==1) printf("%ld day",dias);
  if (dias>1) printf("%ld days",dias);
  if (dias && (horas||minutos||segundos)) printf(", ");
  if (horas==1) printf("%ld hour",horas);
  if (horas>1) printf("%ld hours",horas);
  if (horas && (minutos||segundos)) printf(", ");
  if (minutos==1) printf("%ld minute",minutos);
  if (minutos>1) printf("%ld minutes",minutos);
  if (minutos && segundos) printf(", ");
  if (segundos==1) printf("%ld second",segundos);
  if (segundos>1) printf("%ld seconds",segundos);
}

// given a time in seconds, generate a string containing the time in days, hours, minutes and seconds
void sprinttime(char *s,long t)
{
  long dias,horas,minutos,segundos;
  char aux[100];
  if (t<1) {sprintf(s,"%s","<1 second");return;}
  minutos = (long) t/60;
  segundos = t - minutos*60;
  horas = (long) minutos/60;
  minutos = minutos-horas*60;
  dias = (long) horas/24;
  horas = horas - dias*24;
  s[0]='\0';
  if (dias==1) {sprintf(aux,"%ld day",dias);strcat(s,aux);}
  if (dias>1) {sprintf(aux,"%ld days",dias);strcat(s,aux);}
  if (dias && (horas||minutos||segundos)) strcat(s,", ");
  if (horas==1) {sprintf(aux,"%ld hour",horas);strcat(s,aux);}
  if (horas>1) {sprintf(aux,"%ld hours",horas);strcat(s,aux);}
  if (horas && (minutos||segundos)) strcat(s,", ");
  if (minutos==1) {sprintf(aux,"%ld minute",minutos);strcat(s,aux);}
  if (minutos>1) {sprintf(aux,"%ld minutes",minutos);strcat(s,aux);}
  if (minutos && segundos) strcat(s,", ");
  if (segundos==1) {sprintf(aux,"%ld second",segundos);strcat(s,aux);}
  if (segundos>1) {sprintf(aux,"%ld seconds",segundos);strcat(s,aux);}
}


// print the time in the format dd,hh:mm:ss
void shortprinttime(long t)
{
  long dias,horas,minutos,segundos;
  minutos = (long) t/60;
  segundos = t - minutos*60;
  horas = (long) minutos/60;
  minutos = minutos-horas*60;
  dias = (long) horas/24;
  horas = horas - dias*24;
  if (dias>0) printf("%1ldd,%2ldh:%2ldm:%2lds",dias,horas,minutos,segundos);
  else if (horas>0) printf("%1ldh:%2ldm:%2lds",horas,minutos,segundos);
  else if (minutos>0) printf("%2ldm:%2lds",minutos,segundos);
  else printf("%2lds",segundos);
}


