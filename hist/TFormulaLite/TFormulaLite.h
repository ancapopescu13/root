#ifndef ROOT_TFormulaLite
#define ROOT_TFormulaLite
#ifndef ROOT_TNamed
#include "TNamed.h"
#endif
#include "TInterpreter.h"
#include "TClass.h"
#include <map>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include "Vc/Vc"
#include <Vc/Allocator>
#include <TROOT.h>
#include <list>

template <class T> 
class Initializer
{
	T value;
	TString type;
public:
	Initializer (T arg) { value = arg; }
	TString GetType()
	{
		TString typeName = TString(typeid(value));
		return typeName;
	}
}; 
template <>
class Initializer <double>
{
	double value;
public:
	Initializer (double arg) { value = arg; }
	TString GetTypeDouble() 
	{	
		return "double";
	}
};
template <>
class Initializer <float>
{
	float value;
public:
	Initializer (float arg) { value = arg; }
	TString GetTypeFloat() 
	{
		return "float";
	}
};
template <>
class Initializer <Vc::double_v>
{
	Vc::double_v value;
public:
	Initializer (Vc::double_v arg) { value = arg; }
	TString GetTypeVDouble() 
	{
		return "Vc::double_v";
	}
};
template <>
class Initializer <Vc::float_v>
{
	Vc::float_v value;
public:
	Initializer (Vc::float_v arg) { value = arg; }
	TString GetTypeVFloat() 
	{
		return "Vc::float_v";
	}
};
template <>
class Initializer <std::vector<double>>
{
	std::vector<double> value;
public:
	Initializer (std::vector<double> arg) { value = arg; }
	TString GetTypeVD() 
	{
		return "std::vector<double>";
	}
};
template <>
class Initializer <std::vector<float>>
{
	std::vector<float> value;
public:
	Initializer (std::vector<float> arg) { value = arg; }
	TString GetTypeVF() 
	{
		return "std::vector<float>";
	}
};

template<class T> 
class TFormulaLite
{
private:
   TString                        fName;
   TString                        fForm;
   TString                        fFormula;
   TString                        ind;

   Int_t                          indx;
   Int_t                          pos;
   Int_t                          c;

   std::string                    indexString;

   Bool_t			  fReady = false;
   Bool_t			  fReadyFloat = false;
   Bool_t                         fInjected = false;
   Bool_t                         fPreProcessed = false;
 
   std::vector<double>            paramet;
   std::vector<float>             parametF;
   std::vector<float>             paramValuesF;  
   std::vector<double>            p;
   std::vector<int>               indexing;   
   std::vector<int>               paramPositions;
   std::vector<int>               lengthParam;
   std::vector<int>               lengthIndex;
   std::vector<double>            paramValues;  

public:
   Bool_t                         fReadyToEval = false;
   void*                          fcnptr;

TFormulaLite(const char *name, const char *formula)   :
     fName(name),
     fFormula(formula)
{	
   	if (!fFormula.IsNull()) 
	{ 		
		if( fFormula.Contains("y") || fFormula.Contains("z") || fFormula.Contains("t"))
			fFormula.ReplaceAll("x","x[0]");		
		if( fFormula.Contains("y") )     
			fFormula.ReplaceAll("y","x[1]");		
		if( fFormula.Contains("z") )     
			fFormula.ReplaceAll("z","x[2]");
		if( fFormula.Contains("t") )     
			fFormula.ReplaceAll("t","x[3]");
		fForm = "template<class T> T ";
		fForm += name;
		fForm += "(const T& x) {return " + fFormula + ";}";	
	}
}
const double& GetParameters(const double& params) const
{
	return params;
}
const float& GetParametersFloat(const float& params) const
{
	return params;
}
const char* GetName() const 
{
	return fName;
}
const TString GetSpecialisedName() const 
{
	return fName;
}


void SetParameters(const Int_t& index, const double& params)
{
	for(Int_t i = 0; i < 20; ++i)
	{
		if( i == index )
		{
			c++;
			std::string indexAsString = std::to_string(index);
			std::string indexParameter;
			TString aux = TString::Format("%d", index);
			ind += "p[";
			ind += aux;
			ind += "]";
			lengthIndex.push_back(ind.Length());
			TString par = TString::Format("%f",params);
			lengthParam.push_back(par.Length());
			pos = fFormula.First("p");
			fFormula.Replace(pos, lengthIndex[c-1], par);
			ind.Clear();
		}
	}
}
void SetParametersFloat(const Int_t& index, const float& params)
{
	for(Int_t i = 0; i < 20; ++i)
	{
		if( i == index )
		{
			c++;
			std::string indexAsString = std::to_string(index);
			std::string indexParameter;
			TString aux = TString::Format("%d", index);
			ind += "p[";
			ind += aux;
			ind += "]";
			lengthIndex.push_back(ind.Length());
			TString par = TString::Format("%f",params);
			lengthParam.push_back(par.Length());
			pos = fFormula.First("p");
			fFormula.Replace(pos, lengthIndex[c-1], par);
			fFormula.Replace(pos+lengthParam[c-1]-1, 1, "f");
			ind.Clear();
		}
	}
}


Bool_t ProcessFormula(const std::vector<double>& p, const Int_t dimension, const Int_t nvar = 1)
{
	if(fFormula.Contains("p["))
	{
		c = 0;
		Int_t occ = fFormula.CountChar('p');
		for(Int_t j = 0; j < fFormula.Length(); ++j)
		{
			if( (fFormula[j] == '[') && (fFormula[j-1] == 'p') )
			{
				indexString = fFormula[j+1];
				if ( fFormula[j+2] != ']' )
					indexString += fFormula[j+2];
				indx = std::stoi(indexString);
				indexing.push_back(indx);
			}
		}
		Int_t max = *std::max_element(indexing.begin(), indexing.end());
		if( max >= occ )
			occ = max + 1;
		paramet.resize(occ);
		for(Int_t i = 0; i < fFormula.Length(); ++i)
			if( fFormula[i] == 'p' )
				paramPositions.push_back(i);
		for(Double_t i=0; i < occ; ++i)
		{
			paramet[i] = GetParameters(p[i]);	
			paramValues.push_back(paramet[i]);
		}
		for(Double_t i=0; i < occ; ++i)
			for(Int_t k = 0; k < indexing.size(); ++k)
				if( i == indexing[k] )
					SetParameters(i,paramValues[i]);
	}
	fForm.Clear();
	if( nvar == 1 )
	{
		if( !fPreProcessed )
			if( dimension == 1 )	
			{
				fForm = "template<class T> T ";
				fForm += GetName();
				fForm += "(const T& x) {return " + fFormula + ";}";
				gInterpreter->ProcessLine(fForm);
				fPreProcessed = true;	
			}	
			else
			{
				std::vector<TString> result(dimension);  
				fForm = "template<class T> T ";
				fForm += GetName();
				fForm += "(const T& x) { T result = x; ";
				fForm += "for(Int_t i = 0; i < ";
				fForm += dimension;
				fForm += "; ++i) { result[i] = " ;
				fForm += fFormula;
				fForm += ";} ";
				fForm += "return result;}";
				gInterpreter->ProcessLine(fForm);
				fPreProcessed = true;
			}
	}
	else if( nvar != 1 )
	{
		if( !fPreProcessed )
			if( dimension == 1 )	
			{
				fForm = "template<class T> T ";
				fForm += GetName();
				fForm += "(const T *x) {return " + fFormula + ";}";
				gInterpreter->ProcessLine(fForm);
				fPreProcessed = true;	
			}	
			else
			{
				TString xVariable = "][i]";
				fFormula.ReplaceAll("]", xVariable);
				std::vector<TString> result(dimension);  
				fForm = "template<class T> T ";
				fForm += GetName();
				fForm += "(const T *x) { T result = *x; ";
				fForm += "for(Int_t i = 0; i < ";
				fForm += dimension;
				fForm += "; ++i) { result[i] = " ;
				fForm += fFormula;
				fForm += ";} ";
				fForm += "return result;}";
				gInterpreter->ProcessLine(fForm);
				fPreProcessed = true;
			}
	}	
	fReady = true;
	return fReady;
}
Bool_t ProcessFormulaFloat(const std::vector<float>& p, const Int_t dimension, const Int_t nvar = 1)
{
	if(fFormula.Contains("p["))
	{
		c = 0;
		Int_t occ = fFormula.CountChar('p');
		for(Int_t j = 0; j < fFormula.Length(); ++j)
		{
			if( (fFormula[j] == '[') && (fFormula[j-1] == 'p') )
			{
				indexString = fFormula[j+1];
				if ( fFormula[j+2] != ']' )
					indexString += fFormula[j+2];
				indx = std::stoi(indexString);
				indexing.push_back(indx);
			}
		}
		Int_t max = *std::max_element(indexing.begin(), indexing.end());
		if( max >= occ )
			occ = max + 1;
		parametF.resize(occ);
		for(Int_t i = 0; i < fFormula.Length(); ++i)
			if( fFormula[i] == 'p' )
				paramPositions.push_back(i);
		for(Int_t i=0; i < occ; ++i)
		{
			parametF[i] = GetParametersFloat(p[i]);	
			paramValuesF.push_back(parametF[i]);
		}
		for(Int_t i=0; i < occ; ++i)
			for(Int_t k = 0; k < indexing.size(); ++k)
				if( i == indexing[k] )
					SetParametersFloat(i,paramValuesF[i]);
	}
	fForm.Clear();
	if( nvar == 1 )
	{
		if( !fPreProcessed )
			if( dimension == 1 )	
			{
				fForm = "template<class T> T ";
				fForm += GetName();
				fForm += "(const T& x) {return " + fFormula + ";}";
				gInterpreter->ProcessLine(fForm);
				fPreProcessed = true;	
			}	
			else
			{
				std::vector<TString> result(dimension);  
				fForm = "template<class T> T ";
				fForm += GetName();
				fForm += "(const T& x) { T result = x; ";
				fForm += "for(Int_t i = 0; i < ";
				fForm += dimension;
				fForm += "; ++i) { result[i] = " ;
				fForm += fFormula;
				fForm += ";} ";
				fForm += "return result;}";
				gInterpreter->ProcessLine(fForm);
				fPreProcessed = true;
			}
	}
	else if ( nvar != 1 )
	{
		if( !fPreProcessed )
			if( dimension == 1 )	
			{
				fForm = "template<class T> T ";
				fForm += GetName();
				fForm += "(const T *x) {return " + fFormula + ";}";
				gInterpreter->ProcessLine(fForm);
				fPreProcessed = true;	
			}	
			else
			{
				TString xVariable = "][i]";
				fFormula.ReplaceAll("]", xVariable);
				std::vector<TString> result(dimension);  
				fForm = "template<class T> T ";
				fForm += GetName();
				fForm += "(const T *x) { T result = *x; ";
				fForm += "for(Int_t i = 0; i < ";
				fForm += dimension;
				fForm += "; ++i) { result[i] = " ;
				fForm += fFormula;
				fForm += ";} ";
				fForm += "return result;}";
				gInterpreter->ProcessLine(fForm);
				fPreProcessed = true;
			}
	}
	fReadyFloat = true;
	return fReadyFloat;
}

TString GetTypename(const T& x)
{
	TString type;
	TClass* ptr = TClass::GetClass(typeid(x));
	if(!ptr)
	{
		type = GetPODName(x);
	}
	else
	{
		type = TClass::GetClass(typeid(x))->GetName();
	}

	return type;
}
TString GetPODName(const T& x)
{		
	TString typeName;
	static const std::map<TString,TString> podMap={{"c","char"},
						{"Ds","char16_t"},{"Di","char32_t"},{"w","wchar_t"},{"a","signed char"}, 
						{"s","short"},{"i","int"},{"l","long"},{"x","long long"},
						{"j","unsigned int"},{"h","unsigned char"},{"t","unsigned short int"},{"m","unsigned long 						          int"},{"y","unsigned long long int"},
						{"Ss","string"},{"b","bool"},{"f","float"},{"d","double"},{"e","long double"}
						};	
	try
	{
		typeName = podMap.at(typeid(x).name());
	}
	catch (const std::out_of_range& oor)
	{
		std::cout << "Out of range... " << oor.what() << std::endl;
	}
	return typeName;
}


void* GetSpecialisedFcnPtr(const TString& funcName, const TString& typeName)
{
	const TString textFunction = funcName + "<" + typeName + ">";
	fInjected = true;
	return (void*)gInterpreter->ProcessLine(textFunction);
}


void PreProcessFormula(const T& x, Int_t dimension = 1)
{
	fForm.Clear();
	TString tip = GetTypename(x);
	std::string type;
	type = tip;
	std::string stringSize = "";
	if( dimension != 1 )
		if( type < "std::array" )
		{
			std::string dim = type.substr(type.find(",")+1);
			for(Int_t j = 0; j < dim.length(); ++j)
				if( isdigit(dim[j]) )
					stringSize += dim[j];
			dimension = std::stoi(stringSize);
		}
	if( dimension == 1 )
	{
		if( !fPreProcessed )	
		{
			fForm = "template<class T> T ";
			fForm += GetName();
			fForm += "(const T& x) { return ";
			fForm += fFormula;
			fForm += ";} ";
			if( !fFormula.Contains("p[") )
			{			
				gInterpreter->ProcessLine(fForm);
				fPreProcessed = true;
			}
		}	
	}
	else if( dimension != 1 )	
		if( !fPreProcessed )	
		{
			TString xVariable = "x[i]";
			fFormula.ReplaceAll("x", xVariable);
			std::vector<TString> result(dimension);  
			fForm = "template<class T> T ";
			fForm += GetName();
			fForm += "(const T& x) { T result = x; ";
			fForm += "for(Int_t i = 0; i < ";
			fForm += dimension;
			fForm += "; ++i) { result[i] = " ;
			fForm += fFormula;
			fForm += ";} ";
			fForm += "return result;}";
			if( !fFormula.Contains("p[") )
			{			
				gInterpreter->ProcessLine(fForm);
				fPreProcessed = true;
			}
		}	
}
template <class T2>
class Evaluator
{
	T2 x;
	static void* fcnptr;
public:
	static void* InitializationEval(const T2& x, Bool_t fInjected, TFormulaLite<T2> &fct)
	{
		TString typeName = fct.GetTypename(x);
		TString funcName = fct.GetSpecialisedName();
		return fct.GetSpecialisedFcnPtr(funcName, typeName);
	}
	//Evaluates the formula with a T2 type
	static T2 DoEval(const T2& x, TFormulaLite<T2> &fct)
	{
		fct.PreProcessFormula(x, sizeof(x)/4);
		if ( !fct.fInjected )
		{
			fct.fcnptr = InitializationEval(x, fct.fInjected, fct);
		}
		auto f = (T2 (*)(const T2&))fct.fcnptr;
		return f(x);
	}
	//Evaluates the formula with a T2 type with parameters
	static T2 DoEval(const T2& x, const std::vector<double> p, TFormulaLite<T2>& fct)
	{
		fct.PreProcessFormula(x, sizeof(x)/4);
		if ( !fct.fReadyToEval )
			fct.fReadyToEval = fct.ProcessFormula(p, sizeof(x)/8);
		if( fct.fReadyToEval )
			if ( !fct.fInjected )
			{
				fct.fcnptr = InitializationEval(x, fct.fInjected, fct);
			}
		auto f = (T2 (*)(const T2&, const std::vector<double>&))fct.fcnptr;
		return f(x,p);
	}
	static T2 DoEval(const T2& x, const std::vector<float> p, TFormulaLite<T2>& fct)
	{
		fct.PreProcessFormula(x, sizeof(x)/4);
		if ( !fct.fReadyToEval )
			fct.fReadyToEval = fct.ProcessFormulaFloat(p, sizeof(x)/4);
		if( fct.fReadyToEval )
			if ( !fct.fInjected )
			{
				fct.fcnptr = InitializationEval(x, fct.fInjected, fct);
			}
		auto f = (T2 (*)(const T2&, const std::vector<float>&))fct.fcnptr;
		return f(x,p);
	}
};


void PreProcessFormulaVar(const T *x, Int_t dimension = 1)
{
	fForm.Clear();
	if( dimension == 1 )
	{
		if( !fPreProcessed )
		{
			fForm = "template<class T> T ";
			fForm += GetName();
			fForm += "(const T *x) {return " + fFormula + ";}";
			if( !fFormula.Contains("p[") )
			{
				gInterpreter->ProcessLine(fForm);
				fPreProcessed = true;
			}
		}
	}
	else if( dimension != 1 )
	{
		if( !fPreProcessed )	
		{
			if( !fFormula.Contains("p[") )
			{
				TString xVariable = "][i]";
				fFormula.ReplaceAll("]", xVariable);
				std::vector<TString> result(dimension);  
				fForm = "template<class T> T ";
				fForm += GetName();
				fForm += "(const T *x) { T result = *x; ";
				fForm += "for(Int_t i = 0; i < ";
				fForm += dimension;
				fForm += "; ++i) { result[i] = " ;
				fForm += fFormula;
				fForm += ";} ";
				fForm += "return result;}";
				gInterpreter->ProcessLine(fForm);
				fPreProcessed = true;
			}
		}
	}
}
template <class T3>
class EvaluatorVar
{
	T3 *x;
	static void* fcnptr;
public:
	//Evaluates the formula with multiple variables
	static void* InitializationEvalVar(const T3 *x, Bool_t fInjected, TFormulaLite<T3>& fct)
	{
		TString typeName = fct.GetTypename(*x);
		TString funcName = fct.GetSpecialisedName();
		return fct.GetSpecialisedFcnPtr(funcName, typeName);
	}
	static T3 DoEvalVar(const T3 *x, const Int_t& nvar, TFormulaLite<T3>& fct)
	{
		fct.PreProcessFormulaVar(x, sizeof(x[0])/4);
		if ( !fct.fInjected )
		{
			fct.fcnptr = InitializationEvalVar(x, fct.fInjected, fct);
		}
		auto f = (T3 (*)(const T3*))fct.fcnptr;
		return f(x);
	}
	static T3 DoEvalVar(const T3 *x, const std::vector<double>& p, const Int_t& nvar, TFormulaLite<T3>& fct)
	{
		fct.PreProcessFormulaVar(x);
		if ( !fct.fReadyToEval )
			fct.fReadyToEval = fct.ProcessFormula(p, sizeof(x[0])/8);
		if( fct.fReadyToEval )
			if ( !fct.fInjected )
			{
				fct.fcnptr = InitializationEvalVar(x, fct.fInjected, fct);
			}
		auto f = (T3 (*)(const T3*, const std::vector<double>&))fct.fcnptr;
		return f(x, p);
	}
	static T3 DoEvalVar(const T3 *x, const std::vector<float>& p, const Int_t& nvar, TFormulaLite<T3>& fct)
	{
		fct.PreProcessFormulaVar(x);
		if ( !fct.fReadyToEval )
			fct.fReadyToEval = fct.ProcessFormulaFloat(p, sizeof(x[0])/4);
		if( fct.fReadyToEval )
			if ( !fct.fInjected )
			{
				fct.fcnptr = InitializationEvalVar(x, fct.fInjected, fct);
			}
		auto f = (T3 (*)(const T3*, const std::vector<float>&))fct.fcnptr;
		return f(x, p);
	}
};

	T Eval(const T& x)
	{
		return Evaluator<T>::DoEval(x, *this);	
	}
	T Eval(const T& x, const std::vector<double>& p)
	{
		return Evaluator<T>::DoEval(x, p, *this);	
	}
	T Eval(const T& x, const std::vector<float>& p)
	{
		return Evaluator<T>::DoEval(x, p, *this);	
	}

	T Eval(const T *x, const Int_t& nvar)
	{
		return EvaluatorVar<T>::DoEvalVar(x, nvar, *this);	
	}
	T Eval(const T *x, const std::vector<double>& p, const Int_t& nvar)
	{
		return EvaluatorVar<T>::DoEvalVar(x, p, nvar, *this);	
	}
	T Eval(const T *x, const std::vector<float>& p, const Int_t& nvar)
	{
		return EvaluatorVar<T>::DoEvalVar(x, p, nvar, *this);	
	}
};

template <>
template <class T2> class TFormulaLite<double>::Evaluator
{
	double x;
  	static void* fcnptr;
public:
	static void* InitializationEvalDouble(const double& x, Bool_t fInjected, TFormulaLite<double> &fct)
	{
		Initializer<double> doubleType(x);
		TString typeName = doubleType.GetTypeDouble();
		TString funcName = fct.GetSpecialisedName();
		return fct.GetSpecialisedFcnPtr(funcName, typeName);
	}
	static double DoEval(const double& x, TFormulaLite<double> &fct)
	{
		fct.PreProcessFormula(x);
		if ( !fct.fInjected )
		{
			fct.fcnptr = InitializationEvalDouble(x, fct.fInjected, fct);
		}
		auto f = (double (*)(const double&))fct.fcnptr;
		return f(x);
	}	
	static double DoEval(const double& x, const std::vector<double> p, TFormulaLite<double> &fct)
	{
		fct.PreProcessFormula(x);
		if ( !fct.fReadyToEval )
			fct.fReadyToEval = fct.ProcessFormula(p,1);
		if( fct.fReadyToEval )
			if ( !fct.fInjected )
			{
				fct.fcnptr = InitializationEvalDouble(x, fct.fInjected, fct);
			}
		auto f = (double (*)(const double&, const std::vector<double>&))fct.fcnptr;
		return f(x,p);
	}
};
template <>
template <class T2> class TFormulaLite<float>::Evaluator
{
	float x;
   	static void* fcnptr;
public:
	static void* InitializationEvalFloat(const float& x, Bool_t fInjected, TFormulaLite<float> &fct)
	{
		Initializer<float> floatType(x);
		TString typeName = floatType.GetTypeFloat();
		TString funcName = fct.GetSpecialisedName();
		return fct.GetSpecialisedFcnPtr(funcName, typeName);
	}
	static float DoEval(const float& x, TFormulaLite<float> &fct)
	{
		fct.PreProcessFormula(x);
		if ( !fct.fInjected )
		{
			fct.fcnptr = InitializationEvalFloat(x, fct.fInjected, fct);
		}
		auto f = (float (*)(const float&))fct.fcnptr;
		return f(x);
	}
	static float DoEval(const float& x, const std::vector<double> p, TFormulaLite<float> &fct)
	{
		fct.PreProcessFormula(x);
		if ( !fct.fReadyToEval )
			fct.fReadyToEval = fct.ProcessFormula(p, 1);
		if( fct.fReadyToEval )
			if ( !fct.fInjected )
			{
				fct.fcnptr = InitializationEvalFloat(x, fct.fInjected, fct);
			}
		auto f = (float (*)(const float&, const std::vector<double>&))fct.fcnptr;
		return f(x,p);
	}
};
template <>
template <class T2> class TFormulaLite<Vc::double_v>::Evaluator
{
	Vc::double_v x;
   	static void* fcnptr;
public:
	static void* InitializationEvalDoubleV(const Vc::double_v& x, Bool_t fInjected, TFormulaLite<Vc::double_v> &fct)
	{
		Initializer<Vc::double_v> doubleVType(x);
		TString typeName = doubleVType.GetTypeVDouble();
		TString funcName = fct.GetSpecialisedName();
		return fct.GetSpecialisedFcnPtr(funcName, typeName);
	}
	static Vc::double_v DoEval(const Vc::double_v& x, TFormulaLite<Vc::double_v> &fct)
	{
		fct.PreProcessFormula(x);
		if ( !fct.fInjected )
		{
			fct.fcnptr = InitializationEvalDoubleV(x, fct.fInjected, fct);
		}
		auto f = (Vc::double_v (*)(const Vc::double_v&))fct.fcnptr;
		return f(x);
	}
	static Vc::double_v DoEval(const Vc::double_v& x, const std::vector<double> p, TFormulaLite<Vc::double_v> &fct)
	{
		fct.PreProcessFormula(x);
		if ( !fct.fReadyToEval )
			fct.fReadyToEval = fct.ProcessFormula(p, 1);
		if( fct.fReadyToEval )
			if ( !fct.fInjected )
			{
				fct.fcnptr = InitializationEvalDoubleV(x, fct.fInjected, fct);
			}
		auto f = (Vc::double_v (*)(const Vc::double_v&, const std::vector<double>&))fct.fcnptr;
		return f(x,p);
	}
};
template <>
template <class T2> class TFormulaLite<Vc::float_v>::Evaluator
{
	Vc::float_v x;
   	static void* fcnptr;
public:
	static void* InitializationEvalFloatV(const Vc::float_v& x, Bool_t fInjected, TFormulaLite<Vc::float_v> &fct)
	{
		Initializer<Vc::float_v> floatVType(x);
		TString typeName = floatVType.GetTypeVFloat();
		TString funcName = fct.GetSpecialisedName();
		return fct.GetSpecialisedFcnPtr(funcName, typeName);
	}
	static Vc::float_v DoEval(const Vc::float_v& x, TFormulaLite<Vc::float_v> &fct)
	{
		fct.PreProcessFormula(x);
		if ( !fct.fInjected )
		{
			fct.fcnptr = InitializationEvalFloatV(x, fct.fInjected, fct);
		}
		auto f = (Vc::float_v (*)(const Vc::float_v&))fct.fcnptr;
		return f(x);
	}
	static Vc::float_v DoEval(const Vc::float_v& x, const std::vector<float> p, TFormulaLite<Vc::float_v> &fct)
	{
		fct.PreProcessFormula(x);
		if ( !fct.fReadyToEval )
			fct.fReadyToEval = fct.ProcessFormulaFloat(p, 1);
		if( fct.fReadyToEval )
			if ( !fct.fInjected )
			{
				fct.fcnptr = InitializationEvalFloatV(x, fct.fInjected, fct);
			}
		auto f = (Vc::float_v (*)(const Vc::float_v&, const std::vector<float>&))fct.fcnptr;
		return f(x,p);
	}
};
template <>
template <class T2> class TFormulaLite<std::vector<double>>::Evaluator
{
	std::vector<double> x;
   	static void* fcnptr;
public:
	static void* InitializationEvalDoubleStdV(const std::vector<double>& x, Bool_t fInjected, TFormulaLite<std::vector<double>> &fct)
	{
		Initializer<std::vector<double>> DStdType(x);
		TString typeName = DStdType.GetTypeVD();
		TString funcName = fct.GetSpecialisedName();
		return fct.GetSpecialisedFcnPtr(funcName, typeName);
	}
	static std::vector<double> DoEval(const std::vector<double>& x, TFormulaLite<std::vector<double>> &fct)
	{
		fct.PreProcessFormula(x, x.size());
		if ( !fct.fInjected )
		{
			fct.fcnptr = InitializationEvalDoubleStdV(x, fct.fInjected, fct);
		}
		auto f = (std::vector<double> (*)(const std::vector<double>&))fct.fcnptr;
		return f(x);
	}
	static std::vector<double> DoEval(const std::vector<double>& x, const std::vector<double> p, TFormulaLite<std::vector<double>> &fct)
	{
		fct.PreProcessFormula(x, x.size());
		if ( !fct.fReadyToEval )
			fct.fReadyToEval = fct.ProcessFormula(p, x.size());
		if( fct.fReadyToEval )
			if ( !fct.fInjected )
			{
				fct.fcnptr = InitializationEvalDoubleStdV(x, fct.fInjected, fct);
			}
		auto f = (std::vector<double> (*)(const std::vector<double>&, const std::vector<double>&))fct.fcnptr;
		return f(x,p);
	}
};
template <>
template <class T2> class TFormulaLite<std::vector<float>>::Evaluator
{
	std::vector<float> x;
   	static void* fcnptr;
public:
	static void* InitializationEvalDoubleStdV(const std::vector<float>& x, Bool_t fInjected, TFormulaLite<std::vector<float>> &fct)
	{
		Initializer<std::vector<float>> FStdType(x);
		TString typeName = FStdType.GetTypeVF();
		TString funcName = fct.GetSpecialisedName();
		return fct.GetSpecialisedFcnPtr(funcName, typeName);
	}
	static std::vector<float> DoEval(const std::vector<float>& x, TFormulaLite<std::vector<float>> &fct)
	{
		fct.PreProcessFormula(x, x.size());
		if ( !fct.fInjected )
		{
			fct.fcnptr = InitializationEvalDoubleStdV(x, fct.fInjected, fct);
		}
		auto f = (std::vector<float> (*)(const std::vector<float>&))fct.fcnptr;
		return f(x);
	}
	static std::vector<float> DoEval(const std::vector<float>& x, const std::vector<double> p, TFormulaLite<std::vector<float>> &fct)
	{
		fct.PreProcessFormula(x, x.size());
		if ( !fct.fReadyToEval )
			fct.fReadyToEval = fct.ProcessFormula(p, x.size());
		if( fct.fReadyToEval )
			if ( !fct.fInjected )
			{
				fct.fcnptr = InitializationEvalDoubleStdV(x, fct.fInjected, fct);
			}
		auto f = (std::vector<float> (*)(const std::vector<float>&, const std::vector<double>&))fct.fcnptr;
		return f(x,p);
	}
};
template <>
template <class T3> class TFormulaLite<double>::EvaluatorVar
{
	double *x;
   	static void* fcnptr;
public:
	static void* InitializationEvalVar(const double *x, Bool_t fInjected, TFormulaLite<double>& fct)
	{
		Initializer<double> doubleType(*x);
		TString typeName = doubleType.GetTypeDouble();
		TString funcName = fct.GetSpecialisedName();
		auto fp = fct.GetSpecialisedFcnPtr(funcName, typeName);
		return fp;
	}
	static double DoEvalVar(const double *x, const Int_t& nvar, TFormulaLite<double>& fct)
	{
		fct.PreProcessFormulaVar(x);
		if ( !fct.fInjected )
		{
			fct.fcnptr = InitializationEvalVar(x, fct.fInjected, fct);
		}
		auto f = (double (*)(const double*))fct.fcnptr;
		return f(x);
	}
	static double DoEvalVar(const double *x, const std::vector<double>& p, const Int_t& nvar, TFormulaLite<double>& fct)
	{
		fct.PreProcessFormulaVar(x);
		if ( !fct.fReadyToEval )
			fct.fReadyToEval = fct.ProcessFormula(p, 1, nvar);
		if ( !fct.fInjected )
		{
			fct.fcnptr = InitializationEvalVar(x, fct.fInjected, fct);
		}
		auto f = (double (*)(const double*, const std::vector<double>&))fct.fcnptr;
		return f(x, p);
	}
};
template <>
template <class T3> class TFormulaLite<float>::EvaluatorVar
{
	float *x;

   	static void* fcnptr;
public:
	static void* InitializationEvalVar(const float *x, Bool_t fInjected, TFormulaLite<float>& fct)
	{
		Initializer<float> floatType(*x);
		TString typeName = floatType.GetTypeFloat();
		TString funcName = fct.GetSpecialisedName();
		auto fp = fct.GetSpecialisedFcnPtr(funcName, typeName);
		return fp;
	}
	static float DoEvalVar(const float *x, const Int_t& nvar, TFormulaLite<float>& fct)
	{
		fct.PreProcessFormulaVar(x);
		if ( !fct.fInjected )
		{
			fct.fcnptr = InitializationEvalVar(x, fct.fInjected, fct);
		}
		auto f = (float (*)(const float*))fct.fcnptr;
		return f(x);
	}
	static float DoEvalVar(const float *x, const std::vector<float>& p, const Int_t& nvar, TFormulaLite<float>& fct)
	{
		fct.PreProcessFormulaVar(x);
		if ( !fct.fReadyToEval )
			fct.fReadyToEval = fct.ProcessFormulaFloat(p, 1, nvar);
		if ( !fct.fInjected )
		{
			fct.fcnptr = InitializationEvalVar(x, fct.fInjected, fct);
		}
		auto f = (float (*)(const float*, const std::vector<float>&))fct.fcnptr;
		return f(x, p);
	}
};
template <>
template <class T3> class TFormulaLite<Vc::double_v>::EvaluatorVar
{
	Vc::double_v *x;

   	static void* fcnptr;
public:
	static void* InitializationEvalVar(const Vc::double_v *x, Bool_t fInjected, TFormulaLite<Vc::double_v>& fct)
	{
		Initializer<Vc::double_v> doubleVType(*x);
		TString typeName = doubleVType.GetTypeVDouble();
		TString funcName = fct.GetSpecialisedName();
		auto fp = fct.GetSpecialisedFcnPtr(funcName, typeName);
		return fp;
	}
	static Vc::double_v DoEvalVar(const Vc::double_v *x, const Int_t& nvar, TFormulaLite<Vc::double_v>& fct)
	{
		fct.PreProcessFormulaVar(x);
		if ( !fct.fInjected )
		{
			fct.fcnptr = InitializationEvalVar(x, fct.fInjected, fct);
		}
		auto f = (Vc::double_v (*)(const Vc::double_v*))fct.fcnptr;
		return f(x);
	}
	static Vc::double_v DoEvalVar(const Vc::double_v *x, const std::vector<double>& p, const Int_t& nvar, 
					TFormulaLite<Vc::double_v>& fct)
	{
		fct.PreProcessFormulaVar(x);
		if ( !fct.fReadyToEval )
			fct.fReadyToEval = fct.ProcessFormula(p, 1, nvar);
		if ( !fct.fInjected )
		{
			fct.fcnptr = InitializationEvalVar(x, fct.fInjected, fct);
		}
		auto f = (Vc::double_v (*)(const Vc::double_v*, const std::vector<double>&))fct.fcnptr;
		return f(x, p);
	}
};
template <>
template <class T3> class TFormulaLite<Vc::float_v>::EvaluatorVar
{
	Vc::float_v *x;
   	static void* fcnptr;
public:
	static void* InitializationEvalVar(const Vc::float_v *x, Bool_t fInjected, TFormulaLite<Vc::float_v>& fct)
	{
		Initializer<Vc::float_v> floatVType(*x);
		TString typeName = floatVType.GetTypeVFloat();
		TString funcName = fct.GetSpecialisedName();
		auto fp = fct.GetSpecialisedFcnPtr(funcName, typeName);
		return fp;
	}
	static Vc::float_v DoEvalVar(const Vc::float_v *x, const Int_t& nvar, TFormulaLite<Vc::float_v>& fct)
	{
		fct.PreProcessFormulaVar(x);
		if ( !fct.fInjected )
		{
			fct.fcnptr = InitializationEvalVar(x, fct.fInjected, fct);
		}
		auto f = (Vc::float_v (*)(const Vc::float_v*))fct.fcnptr;
		return f(x);
	}
	static Vc::float_v DoEvalVar(const Vc::float_v *x, const std::vector<float>& p, const Int_t& nvar, 
					TFormulaLite<Vc::float_v>& fct)
	{
		fct.PreProcessFormulaVar(x);
		if ( !fct.fReadyToEval )
			fct.fReadyToEval = fct.ProcessFormulaFloat(p, 1, nvar);
		if ( !fct.fInjected )
		{
			fct.fcnptr = InitializationEvalVar(x, fct.fInjected, fct);
		}
		auto f = (Vc::float_v (*)(const Vc::float_v*, const std::vector<float>&))fct.fcnptr;
		return f(x, p);
	}
};
template <>
template <class T3> class TFormulaLite<std::vector<double>>::EvaluatorVar
{
	std::vector<double> *x;
   	static void* fcnptr;
public:
	static void* InitializationEvalVar(const std::vector<double> *x, Bool_t fInjected, TFormulaLite<std::vector<double>>& fct)
	{
		Initializer<std::vector<double>> doubleStdVType(*x);
		TString typeName = doubleStdVType.GetTypeVD();
		TString funcName = fct.GetSpecialisedName();
		auto fp = fct.GetSpecialisedFcnPtr(funcName, typeName);
		return fp;
	}
	static std::vector<double> DoEvalVar(const std::vector<double> *x, const Int_t& nvar, TFormulaLite<std::vector<double>>& fct)
	{
		fct.PreProcessFormulaVar(x, x[0].size());
		if ( !fct.fInjected )
		{
			fct.fcnptr = InitializationEvalVar(x, fct.fInjected, fct);
		}
		auto f = (std::vector<double> (*)(const std::vector<double> *))fct.fcnptr;
		return f(x);
	}
	static std::vector<double> DoEvalVar(const std::vector<double> *x, const std::vector<double>& p, const Int_t& nvar, 
					TFormulaLite<std::vector<double>>& fct)
	{
		fct.PreProcessFormulaVar(x, x[0].size());
		if ( !fct.fReadyToEval )
			fct.fReadyToEval = fct.ProcessFormula(p,  x[0].size(), nvar);
		if ( !fct.fInjected )
		{
			fct.fcnptr = InitializationEvalVar(x, fct.fInjected, fct);
		}
		auto f = (std::vector<double> (*)(const std::vector<double> *, const std::vector<double>&))fct.fcnptr;
		return f(x, p);
	}
};
template <>
template <class T3> class TFormulaLite<std::vector<float>>::EvaluatorVar
{
	std::vector<float> *x;
   	static void* fcnptr;
public:
	static void* InitializationEvalVar(const std::vector<float> *x, Bool_t fInjected, TFormulaLite<std::vector<float>>& fct)
	{
		Initializer<std::vector<float>> floatStdVType(*x);
		TString typeName = floatStdVType.GetTypeVF();
		TString funcName = fct.GetSpecialisedName();
		auto fp = fct.GetSpecialisedFcnPtr(funcName, typeName);
		return fp;
	}
	static std::vector<float> DoEvalVar(const std::vector<float> *x, const Int_t& nvar, TFormulaLite<std::vector<float>>& fct)
	{
		fct.PreProcessFormulaVar(x, x[0].size());
		if ( !fct.fInjected )
		{
			fct.fcnptr = InitializationEvalVar(x, fct.fInjected, fct);
		}
		auto f = (std::vector<float> (*)(const std::vector<float>*))fct.fcnptr;
		return f(x);
	}
	static std::vector<float> DoEvalVar(const std::vector<float> *x, const std::vector<float>& p, const Int_t& nvar, 
					TFormulaLite<std::vector<float>>& fct)
	{
		fct.PreProcessFormulaVar(x, x[0].size());
		if ( !fct.fReadyToEval )
			fct.fReadyToEval = fct.ProcessFormulaFloat(p,  x[0].size(), nvar);
		if ( !fct.fInjected )
		{
			fct.fcnptr = InitializationEvalVar(x, fct.fInjected, fct);
		}
		auto f = (std::vector<float> (*)(const std::vector<float>*, const std::vector<float>&))fct.fcnptr;
		return f(x, p);
	}
};

#endif
