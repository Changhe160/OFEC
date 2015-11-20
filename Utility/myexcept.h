
/// \file myexcept.h
/// Exception handler.

#ifndef EXCEPTION_LIB
#define EXCEPTION_LIB

#include "include.h"

///
/// exception class thrown by OFEC
///
class myException : public std::exception
{
public :
    ///
    /// construct exception with message
    ///
    /// @param message The exception message.
    ///
    /// The message pointer *must* remain valid during the entire lifetime of the exception.
    ///
    myException(const char* message) : std::exception() { _message = message; }
    ~myException() throw() {}

    ///
    /// @return exception message
    ///
    const char* what() const throw() { return _message; }

private :
    const char* _message;
};


#endif                            // end of EXCEPTION_LIB


// body file: myexcept.cpp


///@}

