#include "typeVar.h"

int TypeVar::operator+(int n) const{
	if (is_int()) return boost::get<int_t>(_var) + n;
	throw myException("invalid + to int@TypeVar::operator+");
}
double TypeVar::operator+(double n) const{
	if (is_double()) return boost::get<double_t>(_var) + n;
	throw myException("invalid +  to double@TypeVar::operator+");
}

TypeVar& TypeVar::operator++(){
	
	if(is_int()){
		_var=boost::get<int>(_var)+1;
	}else if(is_char()){
		_var=static_cast<char>(boost::get<char>(_var)+1);
	}else if(is_bool()){
		_var=1;
	}else{
		throw myException("invalid ++ @TypeVar::operator++()");
	}
	
	return *this;
}
TypeVar TypeVar::operator++(int){
	TypeVar x(*this);
	++(*this);

	return x;
}
TypeVar TypeVar::operator+(const TypeVar& v)const{
	switch (type()) {
    case type_null :    throw myException("invalid + to none");
    case type_bool :    throw myException("bool + not implemented");	
    case type_int :    if(v.is_int()) return boost::get<int_t>(_var) + boost::get<int_t>(v._var) ;
						throw myException("invalid + to int"); 
    case type_double :  if(v.is_double()) return boost::get<double_t>(_var) + boost::get<double_t>(v._var) ;
						throw myException("invalid + to double"); 
    case type_string :  throw myException("string + not implemented");
    case type_wstring : throw myException("wstring + not implemented");
    default :           throw myException("(unhandled type) + not implemented");
    }
}

int TypeVar::operator-(int n ) const{
	if (is_int()) return boost::get<int_t>(_var) - n;
	throw myException("invalid - to int@TypeVar::operator-");
}
double TypeVar::operator-(double n) const{
	if (is_double()) return boost::get<double_t>(_var) - n;
	throw myException("invalid - to double@TypeVar::operator-");
}
TypeVar TypeVar::operator-(const TypeVar& v)const{
		switch (type()) {
    case type_null :    throw myException("invalid - comparison to none");
    case type_bool :    throw myException("bool - not implemented");	
    case type_int :    if(v.is_int()) return boost::get<int_t>(_var) - boost::get<int_t>(v._var) ;
						throw myException("invalid - to int"); 
    case type_double :  if(v.is_double()) return boost::get<double_t>(_var) - boost::get<double_t>(v._var) ;
						throw myException("invalid - to double"); 
    case type_string :  throw myException("string - not implemented");
    case type_wstring : throw myException("wstring - not implemented");
    default :           throw myException("(unhandled type) - not implemented");
    }
}


int TypeVar::operator*(int n) const{
	if (is_int()) return boost::get<int_t>(_var) * n;
	throw myException("invalid * to int@TypeVar::operator*");
}
double TypeVar::operator*(double n) const{
	if (is_double()) return boost::get<double_t>(_var) * n;
	throw myException("invalid * to double@TypeVar::operator*");
}
TypeVar TypeVar::operator*(const TypeVar& v)const{
		switch (type()) {
    case type_null :    throw myException("invalid * comparison to none");
    case type_bool :    throw myException("bool * not implemented");	
    case type_int :    if(v.is_int()) return boost::get<int_t>(_var) * boost::get<int_t>(v._var) ;
						throw myException("invalid * to int"); 
    case type_double :  if(v.is_double()) return boost::get<double_t>(_var) * boost::get<double_t>(v._var) ;
						throw myException("invalid * to double"); 
    case type_string :  throw myException("string * not implemented");
    case type_wstring : throw myException("wstring * not implemented");
    default :           throw myException("(unhandled type) * not implemented");
    }
}


int TypeVar::operator/(int n) const{
	if (is_int()) return boost::get<int_t>(_var) / n;
	throw myException("invalid /  to int@TypeVar::operator/");
}
double TypeVar::operator/(double n) const{
	if (is_double()) return boost::get<double_t>(_var) / n;
	throw myException("invalid /  to double@TypeVar::operator/");
}
TypeVar TypeVar::operator/(const TypeVar& v)const{
		switch (type()) {
    case type_null :    throw myException("invalid /  to none");
    case type_bool :    throw myException("bool / not implemented");	
    case type_int :    if(v.is_int()) return boost::get<int_t>(_var) / boost::get<int_t>(v._var) ;
						throw myException("invalid / to int"); 
    case type_double :  if(v.is_double()) return boost::get<double_t>(_var) / boost::get<double_t>(v._var) ;
						throw myException("invalid / to double"); 
    case type_string :  throw myException("string / not implemented");
    case type_wstring : throw myException("wstring / not implemented");
    default :           throw myException("(unhandled type) / not implemented");
    }
}

	
int TypeVar::operator%(int n) const{
	if (is_int()) return boost::get<int_t>(_var) % n;
	throw myException("invalid % to int@TypeVar::operator%");
}
TypeVar TypeVar::operator%(const TypeVar& v)const{
		switch (type()) {
    case type_null :    throw myException("invalid % to none");
    case type_bool :    throw myException("invalid % not implemented");	
    case type_int :    if(v.is_int()) return boost::get<int_t>(_var) % boost::get<int_t>(v._var) ;
						throw myException("invalid % to int"); 
    case type_double :  throw myException("invalid % to double"); 
    case type_string :  throw myException("invalid % to string");
    case type_wstring : throw myException("invalid % to wstring");
	case type_size_t:  if(v.is_size_t()) return boost::get<size_type>(_var) % boost::get<size_type>(v._var) ;
						throw myException("invalid % to size_t");
    default :           throw myException("(unhandled type) % not implemented");
    }
}

TypeVar::size_type TypeVar::operator-(size_type n) const{
	if (is_size_t()) return boost::get<size_type>(_var) - n;
	throw myException("invalid -  to int@TypeVar::operator-");
}
TypeVar::size_type TypeVar::operator+(size_type n) const{
	if (is_size_t()) return boost::get<size_type>(_var) + n;
	throw myException("invalid +  to int@TypeVar::operator+");
}
TypeVar::size_type TypeVar::operator*(size_type n) const{
	if (is_size_t()) return boost::get<size_type>(_var) * n;
	throw myException("invalid *  to int@TypeVar::operator*");
}
TypeVar::size_type TypeVar::operator/(size_type n) const{
	if (is_size_t()) return boost::get<size_type>(_var) / n;
	throw myException("invalid /  to int@TypeVar::operator/");
}
TypeVar::size_type TypeVar::operator%(size_type n) const{
	if (is_size_t()) return boost::get<size_type>(_var) % n;
	throw myException("invalid %  to int@TypeVar::operator%");
}