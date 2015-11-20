
#include "typeVar.h"
#include "../myexcept.h"
///
/// cast to bool
///
TypeVar::operator bool() const {
    try {
        return boost::get<bool_t>(_var);
    } catch (const boost::bad_get&) {
        throw myException("cannot convert to bool");
    }
}

TypeVar::operator char() const {
    try {
        return boost::get<char_t>(_var);
    } catch (const boost::bad_get&) {
        throw myException("cannot convert to char");
    }
}

///
/// cast to int
///
TypeVar::operator int() const {
    try {
        return boost::get<int_t>(_var);
    } catch (const boost::bad_get&) {
        throw myException("cannot convert to int");
    }
}

///
/// cast to double
///
TypeVar::operator double() const {
    try {
        return boost::get<double_t>(_var);
    } catch (const boost::bad_get&) {
        throw myException("cannot convert to double");
    }
}

///
/// cast to string
///
TypeVar::operator std::string() const {
    try {
        return *boost::get<string_t>(_var).ps;
    } catch (const boost::bad_get&) {
        throw myException("cannot convert to string");
    }
}

///
/// cast to wide string
///
TypeVar::operator std::wstring() const {
    try {
        return *boost::get<wstring_t>(_var).ps;
    } catch (const boost::bad_get&) {
        throw myException("cannot convert to wstring");
    }
}

TypeVar::operator size_type() const {
    try {
		if(is_size_t())      return boost::get<size_type>(_var);
		else if(is_int() ) return static_cast<size_type>(boost::get<int>(_var));
    } catch (const boost::bad_get&) {
        throw myException("cannot convert to size_type");
    }
}
///
/// @return type name
///
struct TypeVar::name_visitor : public boost::static_visitor<std::string>
{
    result_type operator () (const null_t&) const { return "null"; }
    result_type operator () (const bool_t&) const { return "bool"; }
	result_type operator () (const char_t&) const { return "char"; }
    result_type operator () (const int_t&) const { return "int"; }
    result_type operator () (const double_t&) const { return "double"; }
    result_type operator () (const string_t& value) const { return "string"; }
    result_type operator () (const wstring_t& value) const { return "wstring"; }
	result_type operator () (const size_type& ptr) const { return "size_type"; }
};
std::string TypeVar::name() const {
    return boost::apply_visitor(name_visitor(), _var);
}

///
/// @return type identifier
///
TypeVar::code TypeVar::type() const {
    return code(_var.which());
}

