#ifndef DYNAMIC_VAR_HPP
#define DYNAMIC_VAR_HPP

#include "../myexcept.h"
#include "../definition.h"

class TypeVar {
public :
    typedef std::size_t size_type;
    // Note to dynamic developer: Make sure that code and the variant list for var_t always match
    enum code { type_null = 0, type_bool, type_char, type_int, type_double, type_string, type_wstring,type_size_t};

    TypeVar();
    TypeVar(bool);
	TypeVar(char c);
    TypeVar(int n);
    TypeVar(double n);
    TypeVar(const std::string& s);
    TypeVar(const char* s);
    TypeVar(const std::wstring& s);
    TypeVar(const wchar_t* s);
    TypeVar(const TypeVar& v);
	TypeVar(size_type n);

    TypeVar& operator = (bool);
    TypeVar& operator = (char c);
	TypeVar& operator = (int n);
    TypeVar& operator = (double n);
    TypeVar& operator = (const std::string& s);
    TypeVar& operator = (const char* s);
    TypeVar& operator = (const std::wstring& s);
    TypeVar& operator = (const wchar_t* s);
    TypeVar& operator = (const TypeVar& v);
	TypeVar& operator = (size_type n);

	TypeVar& operator += (int n);
    TypeVar& operator += (double n);
	TypeVar& operator -= (int n);
    TypeVar& operator -= (double n);
	TypeVar& operator *= (int n);
    TypeVar& operator *= (double n);
	TypeVar& operator /= (int n);
    TypeVar& operator /= (double n);
	TypeVar& operator %= (int n);

	TypeVar& operator += (size_type n);
	TypeVar& operator -= (size_type n);
	TypeVar& operator *= (size_type n);
	TypeVar& operator /= (size_type n); 
	TypeVar& operator %= (size_type n);


    operator bool() const;
	operator char() const;
    operator int() const;
    operator double() const;
    operator std::string() const;
    operator std::wstring() const;
	operator size_type () const;

    enum code type() const;
    std::string name() const;

	int operator+(int) const;
	double operator+(double) const;
	TypeVar operator+(const TypeVar&)const;

	//waring: there is no range management 
	TypeVar& operator++();
	TypeVar operator++(int);

	int operator-(int) const;
	double operator-(double) const;
	TypeVar operator-(const TypeVar&)const;

	int operator*(int) const;
	double operator*(double) const;
	TypeVar operator*(const TypeVar&)const;

	int operator/(int) const;
	double operator/(double) const;
	TypeVar operator/(const TypeVar&)const;

	
	int operator%(int) const;
	TypeVar operator%(const TypeVar&)const;

	size_type operator-(size_type n) const;
	size_type operator+(size_type n) const;
	size_type operator*(size_type n) const;
	size_type operator/(size_type n) const;
	size_type operator%(size_type n) const;

    bool operator == (bool) const;
	bool operator == (char c) const;
    bool operator == (int n) const;
    bool operator == (double n) const;
    bool operator == (const std::string& s) const;
    bool operator == (const char* s) const;
    bool operator == (const std::wstring& s) const;
    bool operator == (const wchar_t* s) const;
    bool operator == (const TypeVar& v) const;
	bool operator == (size_type n) const;

    bool operator != (bool) const;
	bool operator != (char) const;
    bool operator != (int n) const;
    bool operator != (double n) const;
    bool operator != (const std::string& s) const;
    bool operator != (const char* s) const;
    bool operator != (const std::wstring& s) const;
    bool operator != (const wchar_t* s) const;
    bool operator != (const TypeVar& v) const;
	bool operator != (size_type n) const;


    bool operator < (bool) const;
	bool operator < (char) const;
    bool operator < (int n) const;
    bool operator < (double n) const;
    bool operator < (const std::string& s) const;
    bool operator < (const char* s) const;
    bool operator < (const std::wstring& s) const;
    bool operator < (const wchar_t* s) const;
    bool operator < (const TypeVar& v) const;
	bool operator < (size_type n) const;

    bool operator <= (bool) const;
	 bool operator <= (char) const;
    bool operator <= (int n) const;
    bool operator <= (double n) const;
    bool operator <= (const std::string& s) const;
    bool operator <= (const char* s) const;
    bool operator <= (const std::wstring& s) const;
    bool operator <= (const wchar_t* s) const;
    bool operator <= (const TypeVar& v) const;
	bool operator <= (size_type n) const;

    bool operator > (bool) const;
	bool operator > (char) const;
    bool operator > (int n) const;
    bool operator > (double n) const;
    bool operator > (const std::string& s) const;
    bool operator > (const char* s) const;
    bool operator > (const std::wstring& s) const;
    bool operator > (const wchar_t* s) const;
    bool operator > (const TypeVar& v) const;
	bool operator > (size_type n) const;

    bool operator >= (bool) const;
	bool operator >= (char n) const;
    bool operator >= (int n) const;
    bool operator >= (double n) const;
    bool operator >= (const std::string& s) const;
    bool operator >= (const char* s) const;
    bool operator >= (const std::wstring& s) const;
    bool operator >= (const wchar_t* s) const;
    bool operator >= (const TypeVar& v) const;
    bool operator >= (size_type n) const;    
		
    /// is TypeVar a null?
    bool is_null() const { return type() == type_null; }
    /// is TypeVar a bool?
    bool is_bool() const { return type() == type_bool; }

	 bool is_char() const { return type() == type_char; }
    /// is TypeVar an int?
    bool is_int() const { return type() == type_int; }
    /// is TypeVar a double?
    bool is_double() const { return type() == type_double; }
    /// is TypeVar a numeric type?
    bool is_numeric() const { return is_int() || is_double(); }
    /// is TypeVar a string?
    bool is_string() const { return type() == type_string; }
    /// is TypeVar a wide string?
    bool is_wstring() const { return type() == type_wstring; }
    /// is TypeVar a string type?
    bool is_string_type() const { return is_string() || is_wstring(); }
	
	bool is_size_t() const{ return type()==type_size_t; }

     
    std::ostream& _write_var(std::ostream& os) const;
    std::ostream& _write_string(std::ostream& os) const;
    std::ostream& _write_wstring(std::ostream& os) const;
   
    std::wostream& _write_var(std::wostream& os) const;
    std::wostream& _write_string(std::wostream& os) const;
    std::wostream& _write_wstring(std::wostream& os) const;
        
    size_type count() const;

    /// TypeVar comparison functor
    struct less_var {
        /// TypeVar comparison function
        bool operator () (const TypeVar& lhs, const TypeVar& rhs);
    };


private :
	friend	int operator-(int n, const TypeVar& v);
	friend	double operator-(double n, const TypeVar& v);
	friend	int operator/(int n, const TypeVar& v);
	friend	double operator/(double n, const TypeVar& v);
	friend int operator%(int n, const TypeVar& v);

	friend	int operator+(int n, const TypeVar& v);
	friend	double operator+(double n, const TypeVar& v);
	friend	int operator*(int n, const TypeVar& v);
	friend	double operator*(double n, const TypeVar& v);

	friend	void operator-=(int& n, const TypeVar& v);
	friend	void operator-=(double &n, const TypeVar& v);
	friend	void operator/=(int &n, const TypeVar& v);
	friend	void operator/=(double &n, const TypeVar& v);
	friend void operator%=(int& n, const TypeVar& v);
	friend	void operator+=(int& n, const TypeVar& v);
	friend	void operator+=(double &n, const TypeVar& v);
	friend	void operator*=(int & n, const TypeVar& v);
	friend	void operator*=(double & n, const TypeVar& v);

    typedef boost::blank null_t;

    struct string_t {
        string_t() : ps(boost::make_shared<std::string>()) {}
        string_t(const std::string& s) : ps(boost::make_shared<std::string>(s)) {}
        string_t(const char* s) : ps(boost::make_shared<std::string>(s)) {}

        boost::shared_ptr<std::string>  ps;
    };
        
    struct wstring_t {
        wstring_t() : ps(boost::make_shared<std::wstring>()) {}
        wstring_t(const std::wstring& s) : ps(boost::make_shared<std::wstring>(s)) {}
        wstring_t(const wchar_t* s) : ps(boost::make_shared<std::wstring>(s)) {}

        boost::shared_ptr<std::wstring>  ps;
    };

    typedef bool bool_t;
    typedef int int_t;
    typedef double double_t;
	typedef char char_t;
	


    typedef boost::variant<null_t, bool_t, char_t,int_t, double_t, string_t, wstring_t,size_type> var_t;

    var_t _var;

    struct name_visitor;
    struct count_visitor;
    struct equal_visitor;
	
};

/// ostream << TypeVar
inline std::ostream& operator << (std::ostream& os, const TypeVar& v) { return v._write_var(os); }
/// wostream << TypeVar
inline std::wostream& operator << (std::wostream& os, const TypeVar& v) { return v._write_var(os); }


inline	int operator-(int n, const TypeVar& v){
	if(v.is_int()) return n-boost::get<int>(v._var);
	throw myException("invalid - from int to var @ int operator-");
}
inline	double operator-(double n, const TypeVar& v){
	if(v.is_double()) return n-boost::get<double>(v._var);
	throw myException("invalid - from double to var @ double operator-");
}

inline	int operator/(int n, const TypeVar& v){
	if(v.is_int()) return n/boost::get<int>(v._var);
	throw myException("invalid / from int to var @ int operator /");
}
inline	double operator/(double n, const TypeVar& v){
	if(v.is_double()) return n/boost::get<double>(v._var);
	throw myException("invalid / from double to var @ double operator /");
}

inline	int operator%(int n, const TypeVar& v){
	if(v.is_int()) return n%boost::get<int>(v._var);
	throw myException("invalid % from int to var @ int operator %");
}

inline int operator+(int n, const TypeVar& v){
	if(v.is_int()) return n+boost::get<int>(v._var);
	throw myException("invalid + from int to var @ int operator +");
}
inline double operator+(double n, const TypeVar& v){
	if(v.is_double()) return n+boost::get<double>(v._var);
	throw myException("invalid + from double to var @ double operator +");
}
inline int operator*(int n, const TypeVar& v){
	if(v.is_int()) return n*boost::get<int>(v._var);
	throw myException("invalid * from int to var @ int operator *");
}
inline double operator*(double n, const TypeVar& v){
	if(v.is_double()) return n*boost::get<double>(v._var);
	throw myException("invalid * from double to var @ double operator *");
}

inline void operator-=(int& n, const TypeVar& v){
	if(v.is_int())  n-=boost::get<int>(v._var);
	else throw myException("invalid -= from int to var @ int operator -=");
}
inline void operator-=(double &n, const TypeVar& v){
	if(v.is_double()) n-=boost::get<double>(v._var);
	else throw myException("invalid -= from double to var @ double operator -=");
}
inline void operator/=(int &n, const TypeVar& v){
	if(v.is_int())  n/=boost::get<int>(v._var);
	else throw myException("invalid /= from int to var @ int operator /=");
}
inline void operator/=(double &n, const TypeVar& v){
	if(v.is_double()) n/=boost::get<double>(v._var);
	else throw myException("invalid /= from double to var @ double operator /=");
}
inline void operator%=(int& n, const TypeVar& v){
	if(v.is_int())  n%=boost::get<int>(v._var);
	else throw myException("invalid %= from int to var @ int operator %=");
}
inline void operator+=(int& n, const TypeVar& v){
	if(v.is_int())  n+=boost::get<int>(v._var);
	else throw myException("invalid += from int to var @ int operator +=");
}
inline void operator+=(double &n, const TypeVar& v){
	if(v.is_double()) n+=boost::get<double>(v._var);
	else throw myException("invalid += from double to var @ double operator +=");
}
inline void operator*=(int & n, const TypeVar& v){
	if(v.is_int())  n*=boost::get<int>(v._var);
	else throw myException("invalid *= from int to var @ int operator *=");
}
inline void operator*=(double & n, const TypeVar& v){
	if(v.is_double()) n*=boost::get<double>(v._var);
	else throw myException("invalid *= from double to var @ double operator *=");
}

typedef map<Param,TypeVar> ParamMap;

#endif /* DYNAMIC_VAR_HPP */
