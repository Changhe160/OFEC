#include <stdio.h>
#include <stdarg.h>
#include "LKH.h"

/* 
 * The printff function prints a message and flushes stdout.
 */

void LKH::LKHAlg::printff(char *fmt, ...)
{
    va_list args;

    va_start(args, fmt);
    vprintf(fmt, args);
    va_end(args);
    fflush(stdout);
}
