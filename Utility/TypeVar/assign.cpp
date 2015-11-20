
#include "typeVar.h"

///
/// assign bool to TypeVar
///
TypeVar& TypeVar::operator = (bool n) {
    _var = n;
    return *this;
}

TypeVar& TypeVar::operator = (char c){
	 _var = c;
    return *this;
}
///
/// assign int to TypeVar
///
TypeVar& TypeVar::operator = (int n) {
    _var = n;
    return *this;
}


///
/// assign double to TypeVar
///
TypeVar& TypeVar::operator = (double n) {
    _var = n;
    return *this;
}

///
/// assign string to TypeVar
///
TypeVar& TypeVar::operator = (const std::string& s) {
    _var = string_t(s);
    return *this;
}

///
/// assign string constant to TypeVar
///
TypeVar& TypeVar::operator = (const char* s) {
    _var = string_t(s);
    return *this;
}
    
///
/// assign wide string to TypeVar
///
TypeVar& TypeVar::operator = (const std::wstring& s) {
    _var = wstring_t(s);
    return *this;
}

///
/// assign wide string constant to TypeVar
TypeVar& TypeVar::operator = (const wchar_t* s) {
    _var = wstring_t(s);
    return *this;
}

///
/// assign TypeVar to TypeVar
///
TypeVar& TypeVar::operator = (const TypeVar& v) {
    _var = v._var;
    return *this;
}

TypeVar& TypeVar::operator += (int n){
	if(is_int()) _var=boost::get<int>(_var)+n;
	else throw myException("invalid += to int @TypeVar::operator+=");
	return *this;
}
TypeVar& TypeVar::operator += (double n){
	if(is_double()) _var=boost::get<double>(_var)+n;
	else throw myException("invalid += to double @TypeVar::operator+=");
	return *this;
}
TypeVar& TypeVar::operator -= (int n){
	if(is_int()) _var=boost::get<int>(_var)-n;
	throw myException("invalid -= to int @TypeVar::operator-=");
	return *this;
}
TypeVar& TypeVar::operator -= (double n){
	if(is_double()) _var=boost::get<double>(_var)-n;
	else throw myException("invalid -= to double @TypeVar::operator-=");
	return *this;
}
TypeVar& TypeVar::operator *= (int n){
	if(is_int()) _var=boost::get<int>(_var)*n;
	else throw myException("invalid *= to int @TypeVar::operator*=");
	return *this;
}
TypeVar& TypeVar::operator *= (double n){
	if(is_double()) _var=boost::get<double>(_var)*n;
	else throw myException("invalid *= to double @TypeVar::operator*=");
	return *this;
}
TypeVar& TypeVar::operator /= (int n){
	if(is_int()) _var=boost::get<int>(_var)/n;
	else throw myException("invalid /= to int @TypeVar::operator/=");
	return *this;
}
TypeVar& TypeVar::operator /= (double n){
	if(is_double()) _var=boost::get<double>(_var)/n;
	else throw myException("invalid /= to double @TypeVar::operator/=");
	return *this;
}
TypeVar& TypeVar::operator %= (int n){
	if(is_int()) _var=boost::get<int>(_var)%n;
	else throw myException("invalid %= to int @TypeVar::operator%=");
	return *this;
}

TypeVar& TypeVar::operator = (size_type n){
	 _var = n;
    return *this;
}

TypeVar& TypeVar::operator += (size_type n){
	if(is_size_t()) _var=boost::get<size_type>(_var)+n;
	else throw myException("invalid += to size_type @TypeVar::operator+=");
	return *this;
}
TypeVar& TypeVar::operator -= (size_type n){
	if(is_size_t()) _var=boost::get<size_type>(_var)-n;
	else throw myException("invalid %= to size_type @TypeVar::operator-=");
	return *this;
}
TypeVar& TypeVar::operator *= (size_type n){
	if(is_size_t()) _var=boost::get<size_type>(_var)*n;
	else throw myException("invalid %= to size_type @TypeVar::operator*=");
	return *this;
}
TypeVar& TypeVar::operator /= (size_type n){
	if(is_size_t()) _var=boost::get<size_type>(_var)/n;
	else throw myException("invalid %= to size_type @TypeVar::operator/=");
	return *this;
}
TypeVar& TypeVar::operator %= (size_type n){
	if(is_size_t()) _var=boost::get<size_type>(_var)%n;
	else throw myException("invalid %= to size_type @TypeVar::operator%=");
	return *this;
}