#ifndef MYUTILS_DEFINE
#define MYUTILS_DEFINE
#include <string>
#include <iostream>  
#include <time.h>
#include <math.h>
using namespace std;
#include <ostream>
#include <sstream>

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


//====================================================================================================
//     * Type Conversion Routines
// In C++11 we have new convertion functions, but some libraries
// are still not compatible, like the Gurobi Library.

/* float              stof(const string& str, size_t *idx = 0); */
/* double             stod(const string& str, size_t *idx = 0); */
/* long double        stold(const string& str, size_t *idx = 0); */
/* int                stoi(const string& str, size_t *idx = 0, int base = 10); */
/* long               stol(const string& str, size_t *idx = 0, int base = 10); */
/* unsigned long      stoul(const string& str, size_t *idx = 0, int base = 10); */
/* long long          stoll(const string& str, size_t *idx = 0, int base = 10); */
/* unsigned long long stoull(const string& str, size_t *idx = 0, int base = 10); */
/* string to_string(int val); */
/* string to_string(unsigned val); */
/* string to_string(long val); */
/* string to_string(unsigned long val); */
/* string to_string(long long val); */
/* string to_string(unsigned long long val); */
/* string to_string(float val); */
/* string to_string(double val); */
/* string to_string(long double val); */

/* template <typename T> string ToString ( T val ){ */
/*     std::stringstream retval; */
/*     retval << val; */
/*     return retval; */
/* } */



// convert a double value to string
inline std::string DoubleToString(double val) 
{ std::stringstream out; out << val; return out.str(); } 
//{return to_string(val);} // this is only valid in c++11 


// convert a int value to string
inline std::string IntToString(int val)
{ std::stringstream out; out << val; return out.str(); }
//{return to_string(val);}  

// convert a string to double
inline double StringToDouble(string s)
{ double d;  stringstream(s) >> d;  return(d); }
//{return stod(s);} 

// convert a string to int
inline int StringToInt(string s)
{ int n;  stringstream(s) >> n;  return(n); }
//{return stoi(s);} 

inline void Pause(void) {cout<<"Pause";std::cin.get();cout<<"\n";}

//====================================================================================
//     * Functions to test values
bool IsFrac(double x);


//====================================================================================
//     * Working with colors
// colors must match the ColorName in myutils.cpp
typedef enum Color {NOCOLOR,WHITE,BLACK,RED,GREEN,BLUE,YELLOW,MAGENTA,CYAN,GRAY,ORANGE} Color;

// Given a color code, return its name
std::string ColorName(int cor);

//================================================================================================================
//    ROUTINES FOR TIME MANIPULATION

long time70(void);  /* returns the time in seconds since Jan 1, 1970 */
void printtime(long t); /* print the time t in days, hours, minutes and seconds.*/
void sprinttime(char *s,long t); /* prints the time in the string s Example: 1 hour, 2 minutes, 3 seconds*/
void shortprinttime(long t); /* prints the time in the string s. Ex.: 11d,22h:33m:44s   */

#endif
