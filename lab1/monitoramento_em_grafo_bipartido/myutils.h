#ifndef MYUTILS_DEFINE
#define MYUTILS_DEFINE
#include <string>
#include <iostream>  
#include <time.h>
#include <math.h>
using namespace std;

//================================================================================================================
// SOME CONSTANTS
//
// maximum number of characters used by a command used by the system routine.
#define MAXCOMMANDSIZE 1000
// used to verify if number are the same, within some error of MY_EPS
#define MY_EPS 0.0001

//================================================================================================================
//    ROUTINES TO DEAL WITH PDF
// set default pdf reader used in some programs
// For example, if it is possible to open a pdf file
// with the command "open file.pdf"
// then, you can use the following command to define
// the "open" as a pdf reader
// set_pdfreader("open");
void set_pdfreader(string programname);

// To see a pdf file. It uses the pdf reader defined by set_pdfreader.
int view_pdffile(string filename);


//================================================================================================================
//    ROUTINES FOR TYPE CONVERSION

// convert a double value to string
std::string DoubleToString(double val);

// convert a string to double
double StringToDouble(string s);

// convert a string to int
int StringToInt(string s);

// convert a int value to string
std::string IntToString(int val);

inline void Pause(void) {cout<<"Pause";std::cin.get();cout<<"\n";}

bool IsFrac(double x);


// colors must match the ColorName in myutils.cpp
typedef enum Color {NOCOLOR,WHITE,BLACK,RED,GREEN,BLUE,YELLOW,MAGENTA,CYAN,GRAY,ORANGE} Color;

// Given a color code, return its name
std::string ColorName(int cor);

//================================================================================================================
//    ROUTINES FOR TIME MANIPULATION

long time70(void);  /* retorna o tempo em segundos desde 1/jan/1970*/
void printtime(long t); /* imprime o tempo em dias, horas, minutos, seg.*/

void sprinttime(char *s,long t); /* imprime o tempo na string s */
void shortprinttime(long t); /* imprime o tempo na string s */

#endif
