
#include "typeVar.h"

///
/// TypeVar == bool
///
bool TypeVar::operator == (bool n) const {
    return *this == TypeVar(n);
}

bool TypeVar::operator == (char c) const {
    return *this == TypeVar(c);
}

bool TypeVar::operator == (size_type n) const{
	return *this == TypeVar(n);
}
///
/// TypeVar == int
///
bool TypeVar::operator == (int n) const {
    return *this == TypeVar(n);
}

///
/// TypeVar == double
///
bool TypeVar::operator == (double n) const {
    return *this == TypeVar(n);
}

///
/// TypeVar == string
///
bool TypeVar::operator == (const std::string& s) const {
    return *this == TypeVar(s);
}

///
/// TypeVar == string constant
///
bool TypeVar::operator == (const char* s) const {
    return *this == TypeVar(s);
}

///
/// TypeVar == wide string
///
bool TypeVar::operator == (const std::wstring& s) const {
    return *this == TypeVar(s);
}

///
/// TypeVar == wide string constant
///
bool TypeVar::operator == (const wchar_t* s) const {
    return *this == TypeVar(s);
}

///
/// TypeVar == TypeVar
///
struct TypeVar::equal_visitor : public boost::static_visitor<bool>
{
    // Different types
    template <typename T, typename U>
    result_type operator () (const T& lhs, const U& rhs) const
    {
        return false;
    }
    // Same types
    template <typename T>
    result_type operator () (const T& lhs, const T& rhs) const
    {
        // Use overloaded functions to handle explicit specialization
        return equal(lhs, rhs);
    }

private:
    bool equal(const null_t&, const null_t&) const
    {
        return true;
    }
    bool equal(const bool_t& lhs, const bool_t& rhs) const
    {
        return lhs == rhs;
    }
	bool equal(const char_t& lhs, const char_t& rhs) const
    {
        return lhs == rhs;
    }
    bool equal(const int_t& lhs, const int_t& rhs) const
    {
        return lhs == rhs;
    }
    bool equal(const double_t& lhs, const double_t& rhs) const
    {
        return lhs == rhs;
    }
    bool equal(const string_t& lhs, const string_t& rhs) const
    {
        return *(lhs.ps) == *(rhs.ps);
    }
    bool equal(const wstring_t& lhs, const wstring_t& rhs) const
    {
        return *(lhs.ps) == *(rhs.ps);
    }
	bool equal(const size_type& lhs, const size_type& rhs) const
    {
        return lhs == rhs;
    }
};
bool TypeVar::operator == (const TypeVar& v) const {
    return boost::apply_visitor(equal_visitor(), _var, v._var);
}

///
/// TypeVar != bool
///
bool TypeVar::operator != (bool n) const {
    return *this != TypeVar(n);
}

bool TypeVar::operator != (char n) const {
    return *this != TypeVar(n);
}
///
/// TypeVar != size_type
///
bool TypeVar::operator != (size_type n) const {
    return *this != TypeVar(n);
}

///
/// TypeVar != int
///
bool TypeVar::operator != (int n) const {
    return *this != TypeVar(n);
}

///
/// TypeVar != double
///
bool TypeVar::operator != (double n) const {
    return *this != TypeVar(n);
}

///
/// TypeVar != string
///
bool TypeVar::operator != (const std::string& s) const {
    return *this != TypeVar(s);
}

///
/// TypeVar != string constant
///
bool TypeVar::operator != (const char* s) const {
    return *this != TypeVar(s);
}

///
/// TypeVar != wide string
///
bool TypeVar::operator != (const std::wstring& s) const {
    return *this != TypeVar(s);
}

///
/// TypeVar != wide string constant
///
bool TypeVar::operator != (const wchar_t* s) const {
    return *this != TypeVar(s);
}

///
/// TypeVar != TypeVar
///
bool TypeVar::operator != (const TypeVar& v) const {
    return !boost::apply_visitor(equal_visitor(), _var, v._var);
}

///
/// TypeVar < bool
///
bool TypeVar::operator < (bool n) const {
    const TypeVar v(n);
    return less_var()(*this, v);
}

bool TypeVar::operator < (char n) const {
    const TypeVar v(n);
    return less_var()(*this, v);
}
///
/// TypeVar < int
///
bool TypeVar::operator < (int n) const {
    const TypeVar v(n);
    return less_var()(*this, v);
}

///
/// TypeVar < size_t
///
bool TypeVar::operator < (size_type n) const {
    const TypeVar v(n);
    return less_var()(*this, v);
}

///
/// TypeVar < double
///
bool TypeVar::operator < (double n) const {
    const TypeVar v(n);
    return less_var()(*this, v);
}

///
/// TypeVar < string
///
bool TypeVar::operator < (const std::string& s) const {
    const TypeVar v(s);
    return less_var()(*this, v);
}

///
/// TypeVar < string constant
///
bool TypeVar::operator < (const char* s) const {
    const TypeVar v(s);
    return less_var()(*this, v);
}

///
/// TypeVar < wide string
///
bool TypeVar::operator < (const std::wstring& s) const {
    const TypeVar v(s);
    return less_var()(*this, v);
}

///
/// TypeVar < wide string constant
///
bool TypeVar::operator < (const wchar_t* s) const {
    const TypeVar v(s);
    return less_var()(*this, v);
}

///
/// TypeVar < TypeVar
///
bool TypeVar::operator < (const TypeVar& v) const {
    return less_var()(*this, v);
}

///
/// TypeVar <= bool
///
bool TypeVar::operator <= (bool n) const {
    if (is_bool()) return boost::get<bool_t>(_var) <= n;
    throw myException("invalid <= comparison to bool");
}

bool TypeVar::operator <= (char n) const {
    if (is_char()) return boost::get<char_t>(_var) <= n;
    throw myException("invalid <= comparison to char");
}

///
/// TypeVar <= int
///
bool TypeVar::operator <= (int n) const {
    if (is_int()) return boost::get<int_t>(_var) <= n;
    throw myException("invalid <= comparison to int");
}

///
/// TypeVar <= double
///
bool TypeVar::operator <= (double n) const {
    if (is_double()) return boost::get<double_t>(_var) <= n;
    throw myException("invalid <= comparison to double");
}

///
/// TypeVar <= string
///
bool TypeVar::operator <= (const std::string& s) const {
    if (is_string()) return *boost::get<string_t>(_var).ps <= s;
    throw myException("invalid <= comparison to string");
}

///
/// TypeVar <= string constant
///
bool TypeVar::operator <= (const char* s) const {
    if (is_string()) return *boost::get<string_t>(_var).ps <= s;
    throw myException("invalid <= comparison to char*");
}
    
///
/// TypeVar <= wide string
///
bool TypeVar::operator <= (const std::wstring& s) const {
    if (is_wstring()) return *boost::get<wstring_t>(_var).ps <= s;
    throw myException("invalid <= comparison to wstring");
}

///
/// TypeVar <= wide string constant
///
bool TypeVar::operator <= (const wchar_t* s) const {
    if (is_wstring()) return *boost::get<wstring_t>(_var).ps <= s;
    throw myException("invalid <= comparison to wchar_t*");
}

///
/// TypeVar <= TypeVar
///
bool TypeVar::operator <= (const TypeVar& v) const {
    switch (type()) {
    case type_null :    throw myException("invalid <= comparison to none");
	case type_char:		return v.is_char() && boost::get<char_t>(_var) <= boost::get<char_t>(v._var);
    case type_bool :    return v.is_bool() && boost::get<bool_t>(_var) <= boost::get<bool_t>(v._var);
    case type_int :     return v.is_int() && boost::get<int_t>(_var) <= boost::get<int_t>(v._var);
    case type_double :  return v.is_double() && boost::get<double_t>(_var) <= boost::get<double_t>(v._var);
    case type_string :  return v.is_string() && *boost::get<string_t>(_var).ps <= *boost::get<string_t>(v._var).ps;
    case type_wstring : return v.is_wstring() && *boost::get<wstring_t>(_var).ps <= *boost::get<wstring_t>(v._var).ps;
	case type_size_t:   return v.is_size_t() && boost::get<size_type>(_var) <= boost::get<size_type>(v._var);
    default :           throw myException("(unhandled type) <= not implemented");
    }
}

///
/// TypeVar > bool
///
bool TypeVar::operator > (bool n) const {
    if (is_bool()) return boost::get<bool_t>(_var) > n;
    throw myException("invalid > comparison to bool");
}

bool TypeVar::operator > (char n) const {
    if (is_char()) return boost::get<char_t>(_var) > n;
    throw myException("invalid > comparison to bool");
}
///
/// TypeVar > int
///
bool TypeVar::operator > (int n) const {
    if (is_int()) return boost::get<int_t>(_var) > n;
    throw myException("invalid > comparison to int");
}

///
/// TypeVar > double
///
bool TypeVar::operator > (double n) const {
    if (is_double()) return boost::get<double_t>(_var) > n;
    throw myException("invalid > comparison to double");
}

///
/// TypeVar > string
///
bool TypeVar::operator > (const std::string& s) const {
    if (is_string()) return *boost::get<string_t>(_var).ps > s;
    throw myException("invalid > comparison to string");
}

///
/// TypeVar > string constant
bool TypeVar::operator > (const char* s) const {
    if (is_string()) return *boost::get<string_t>(_var).ps > s;
    throw myException("invalid > comparison to char*");
}

///
/// TypeVar > wide string
///
bool TypeVar::operator > (const std::wstring& s) const {
    if (is_wstring()) return *boost::get<wstring_t>(_var).ps > s;
    throw myException("invalid > comparison to wstring");
}

///
/// TypeVar > wide string constant
///
bool TypeVar::operator > (const wchar_t* s) const {
    if (is_wstring()) return *boost::get<wstring_t>(_var).ps > s;
    throw myException("invalid > comparison to wchar_t*");
}

///
/// TypeVar > TypeVar
///
bool TypeVar::operator > (const TypeVar& v) const {
    switch (type()) {
    case type_null :    throw myException("invalid > comparison to none");
	case type_char :    return v.is_char() && boost::get<char_t>(_var) > boost::get<char_t>(v._var);
    case type_bool :    return v.is_bool() && boost::get<bool_t>(_var) > boost::get<bool_t>(v._var);
    case type_int :     return v.is_int() && boost::get<int_t>(_var) > boost::get<int_t>(v._var);
    case type_double :  return v.is_double() && boost::get<double_t>(_var) > boost::get<double_t>(v._var);
    case type_string :  return v.is_string() && *boost::get<string_t>(_var).ps > *boost::get<string_t>(v._var).ps;
    case type_wstring : return v.is_wstring() && *boost::get<wstring_t>(_var).ps > *boost::get<wstring_t>(v._var).ps;
	case type_size_t:   return v.is_size_t() && boost::get<size_type>(_var) > boost::get<size_type>(v._var);
    default :           throw myException("(unhandled type) > not implemented");
    }
}

///
/// TypeVar >= bool
///
bool TypeVar::operator >= (bool n) const {
    if (is_bool()) return boost::get<bool_t>(_var) >= n;
    throw myException("invalid >= comparison to bool");
}

bool TypeVar::operator >= (char n) const {
    if (is_char()) return boost::get<char_t>(_var) >= n;
    throw myException("invalid >= comparison to char");
}


///
/// TypeVar >= int
///
bool TypeVar::operator >= (int n) const {
    if (is_int()) return boost::get<int_t>(_var) >= n;
    throw myException("invalid >= comparison to int");
}

///
/// TypeVar >= double
///
bool TypeVar::operator >= (double n) const {
    if (is_double()) return boost::get<double_t>(_var) >= n;
    throw myException("invalid >= comparison to double");
}

///
/// TypeVar >= string
///
bool TypeVar::operator >= (const std::string& s) const {
    if (is_string()) return *boost::get<string_t>(_var).ps >= s;
    throw myException("invalid >= comparison to string");
}

///
/// TypeVar >= string constant
///
bool TypeVar::operator >= (const char* s) const {
    if (is_string()) return *boost::get<string_t>(_var).ps >= s;
    throw myException("invalid >= comparison to char*");
}
    
///
/// TypeVar >= wide string
///
bool TypeVar::operator >= (const std::wstring& s) const {
    if (is_wstring()) return *boost::get<wstring_t>(_var).ps >= s;
    throw myException("invalid >= comparison to wstring");
}

///
/// TypeVar >= wide string constant
///
bool TypeVar::operator >= (const wchar_t* s) const {
    if (is_wstring()) return *boost::get<wstring_t>(_var).ps >= s;
    throw myException("invalid >= comparison to wchar_t*");
}

///
/// TypeVar >= TypeVar
///
bool TypeVar::operator >= (const TypeVar& v) const {
    switch (type()) {
    case type_null :    throw myException("invalid >= comparison to none");
    case type_bool :    return v.is_bool() && boost::get<bool_t>(_var) >= boost::get<bool_t>(v._var);
	case type_char :     return v.is_int() && boost::get<char_t>(_var) >= boost::get<char_t>(v._var);
    case type_int :     return v.is_int() && boost::get<int_t>(_var) >= boost::get<int_t>(v._var);
    case type_double :  return v.is_double() && boost::get<double_t>(_var) >= boost::get<double_t>(v._var);
    case type_string :  return v.is_string() && *boost::get<string_t>(_var).ps >= *boost::get<string_t>(v._var).ps;
    case type_wstring : return v.is_wstring() && *boost::get<wstring_t>(_var).ps >= *boost::get<wstring_t>(v._var).ps;
	case type_size_t:   return v.is_size_t() && boost::get<size_type>(_var) >= boost::get<size_type>(v._var);
    default :           throw myException("(unhandled type) >= not implemented");
    }
}

bool TypeVar::operator >= (size_type n) const{
	if (is_size_t()) return boost::get<size_type>(_var) >= n;
    throw myException("invalid >= comparison to size_t");
} 
  
bool TypeVar::operator <= (size_type n) const{
	if (is_size_t()) return boost::get<size_type>(_var) <= n;
    throw myException("invalid <= comparison to size_t");
} 
bool TypeVar::operator > (size_type n) const{
	if (is_size_t()) return boost::get<size_type>(_var) > n;
    throw myException("invalid > comparison to size_t");
} 
