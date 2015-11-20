
#include "typeVar.h"

namespace TypeVar_None{
const TypeVar typeNone;
};
bool TypeVar::less_var::operator () (const TypeVar& lhs, const TypeVar& rhs) {
    // if the two vars are of different types, order by type
    code lht = lhs.type(), rht = rhs.type();
    if (lht != rht) return lht < rht;

    // they are of the same type, order by value
    switch (lht) {
    case type_null : return false;
	case type_char : return boost::get<char_t>(lhs._var) < boost::get<char_t>(rhs._var);
    case type_bool : return boost::get<bool_t>(lhs._var) < boost::get<bool_t>(rhs._var);
    case type_int : return boost::get<int_t>(lhs._var) < boost::get<int_t>(rhs._var);
    case type_double : return boost::get<double_t>(lhs._var) < boost::get<double_t>(rhs._var);
    case type_string : return *(boost::get<string_t>(lhs._var).ps) < *(boost::get<string_t>(rhs._var).ps);
    case type_wstring : return *(boost::get<wstring_t>(lhs._var).ps) < *(boost::get<wstring_t>(rhs._var).ps);
    default : throw myException("unhandled type");
    }
}
/*
///
/// append a bool to a collection
///
TypeVar& TypeVar::operator () (bool n) { return operator() (TypeVar(n)); }

TypeVar& TypeVar::operator () (char n) { return operator() (TypeVar(n)); }

///
/// append an int to a collection
///
TypeVar& TypeVar::operator () (int n) { return operator() (TypeVar(n)); }

///
/// append a double to a collection
///
TypeVar& TypeVar::operator () (double n) { return operator() (TypeVar(n)); }

///
/// append a string to a collection
///
TypeVar& TypeVar::operator () (const std::string& s) { return operator() (TypeVar(s)); }

///
/// append a string constant to a collection
///
TypeVar& TypeVar::operator () (const char* s) { return operator() (TypeVar(s)); }

///
/// append a wide string to a collection
///
TypeVar& TypeVar::operator () (const std::wstring& s) { return operator() (TypeVar(s)); }

///
/// append a wide string constant to a collection
///
TypeVar& TypeVar::operator () (const wchar_t* s) { return operator() (TypeVar(s)); }

*/
///
/// count of objects in a collection or characters in a string
///
struct TypeVar::count_visitor : public boost::static_visitor<size_type>
{
    result_type operator () (const null_t&) const { return 0; }
	result_type operator () (const char_t&) const { return 1; }
    result_type operator () (const bool_t&) const { return 1; }
    result_type operator () (const int_t&) const { return 1; }
    result_type operator () (const double_t&) const { return 1; }
	result_type operator () (const size_type&) const { return 1; }
    result_type operator () (const string_t& value) const { return value.ps->length(); }
    result_type operator () (const wstring_t& value) const { return value.ps->length(); }
};
TypeVar::size_type TypeVar::count() const {
    return boost::apply_visitor(count_visitor(), _var);
}
/*
///
/// index a collection with a double
///
TypeVar& TypeVar::operator [] (double n) { return operator[] (TypeVar(n)); }

///
/// index a collection with a string
///
TypeVar& TypeVar::operator [] (const std::string& s) { return operator[] (TypeVar(s)); }

///
/// index a collection with a string constant
///
TypeVar& TypeVar::operator [] (const char* s) { return operator[] (TypeVar(s)); }

///
/// index a collection with a wide string
///
TypeVar& TypeVar::operator [] (const std::wstring& s) { return operator[] (TypeVar(s)); }

///
/// index a collection with a wide string constant
///
TypeVar& TypeVar::operator [] (const wchar_t* s) { return operator[] (TypeVar(s)); }
    

	*/
///
/// write a TypeVar to an ostream
///
std::ostream& TypeVar::_write_var(std::ostream& os) const {
    switch (type()) {
    case type_null :    os << "null"; return os;
    case type_bool:     os << (boost::get<bool_t>(_var) ? "true" : "false"); return os;
	case type_char:     os << boost::get<char_t>(_var); return os;
    case type_int :     os << boost::get<int_t>(_var); return os;
    case type_double :  os << boost::get<double_t>(_var); return os;
    case type_string :  return _write_string(os);
    case type_wstring : return _write_wstring(os);
    default :           throw myException("TypeVar::_write_var(ostream) unhandled type");
    }
}

///
/// write a string to an ostream
///
std::ostream& TypeVar::_write_string(std::ostream& os) const {
    assert(is_string());
    //os << '"';
    for (const char* s = (*boost::get<string_t>(_var).ps).c_str(); *s; ++s)
        switch (*s) {
        case '\b' : os << "\\b"; break;
        case '\r' : os << "\\r"; break;
        case '\n' : os << "\\n"; break;
        case '\f' : os << "\\f"; break;
        case '\t' : os << "\\t"; break;
        case '\\' : os << "\\\\"; break;
        case '\"' : os << "\\\""; break;
        //case '/' : os << "\\/"; break;
        default :
            if (std::iscntrl(*s)) os << "0" << std::oct << std::setw(3) << std::setfill('0') << int(*s);
            else os << *s;
        }
    //os << '"';
    return os;
}

///
/// write a wide string to an ostream
///
std::ostream& TypeVar::_write_wstring(std::ostream& os) const {
    assert(is_wstring());
    //os << '"';
    for (const wchar_t* s = (*boost::get<wstring_t>(_var).ps).c_str(); *s; ++s)
        switch (*s) {
        case L'\b' : os << L"\\b"; break;
        case L'\r' : os << L"\\r"; break;
        case L'\n' : os << L"\\n"; break;
        case L'\f' : os << L"\\f"; break;
        case L'\t' : os << L"\\t"; break;
        case L'\\' : os << L"\\\\"; break;
        case L'\"' : os << L"\\\""; break;
        //case L'/' : os << L"\\/"; break;
        default :
            if (std::iswcntrl(*s)) os << L"0" << std::oct << std::setw(3) << std::setfill('0') << int(*s);
            else os << *s;
        }
    //os << '"';
    return os;
}


///
/// write a TypeVar to a wostream
///
std::wostream& TypeVar::_write_var(std::wostream& os) const {
    switch (type()) {
    case type_null :    os << "null"; return os;
    case type_bool:     os << (boost::get<bool_t>(_var) ? "true" : "false"); return os;
	case type_char :     os << boost::get<char_t>(_var); return os;
    case type_int :     os << boost::get<int_t>(_var); return os;
    case type_double :  os << boost::get<double_t>(_var); return os;
    case type_string :  return _write_string(os);
    case type_wstring : return _write_wstring(os);
    default :           throw myException("TypeVar::_write_var(wostream) unhandled type");
    }
}

///
/// write a string to a wostream
///
std::wostream& TypeVar::_write_string(std::wostream& os) const {
    assert(is_string());
    os << '\'';
    for (const char* s = (*boost::get<string_t>(_var).ps).c_str(); *s; ++s)
        switch (*s) {
        case '\b' : os << "\\b"; break;
        case '\r' : os << "\\r"; break;
        case '\n' : os << "\\n"; break;
        case '\f' : os << "\\f"; break;
        case '\t' : os << "\\t"; break;
        case '\\' : os << "\\\\"; break;
        case '\'' : os << "\\'"; break;
        default :
            if (*s < ' ') os << "0" << std::oct << std::setw(3) << std::setfill(L'0') << int(*s);
            else os << *s;
        }
    os << '\'';
    return os;
}

///
/// write a wide string to a wostream
///
std::wostream& TypeVar::_write_wstring(std::wostream& os) const {
    assert(is_wstring());
    os << '\'';
    for (const wchar_t* s = (*boost::get<wstring_t>(_var).ps).c_str(); *s; ++s)
        switch (*s) {
        case '\b' : os << L"\\b"; break;
        case '\r' : os << L"\\r"; break;
        case '\n' : os << L"\\n"; break;
        case '\f' : os << L"\\f"; break;
        case '\t' : os << L"\\t"; break;
        case '\\' : os << L"\\\\"; break;
        case '\'' : os << L"\\'"; break;
        default :
            if (*s < ' ') os << "0" << std::oct << std::setw(3) << std::setfill(L'0') << int(*s);
            else os << *s;
        }
    os << '\'';
    return os;
}
    