#include "typeVar.h"

///
/// ctor: init with none
///
TypeVar::TypeVar() : _var() {}

///
/// ctor: init with bool
///
TypeVar::TypeVar(bool n) : _var(n) {}

///
/// ctor: init with int
///
TypeVar::TypeVar(int n) : _var(n) {}

///
/// ctor: init with double
///
TypeVar::TypeVar(double n) : _var(n) {}

///
/// ctor: init with string
///
TypeVar::TypeVar(const std::string& s) : _var(string_t(s)) {}

///
/// ctor: init with string constant
///
TypeVar::TypeVar(const char* s) : _var(string_t(s)) {}

///
/// ctor: init with wide string
///
TypeVar::TypeVar(const std::wstring& s) : _var(wstring_t(s)) {}

///
/// ctor: init with wide string
///
TypeVar::TypeVar(const wchar_t* s) : _var(wstring_t(s)) {}

///
/// ctor: init with TypeVar
///
TypeVar::TypeVar(const TypeVar& v) : _var(v._var) {}

TypeVar::TypeVar(char c):_var(c){}

TypeVar::TypeVar(size_type n):_var(n){}