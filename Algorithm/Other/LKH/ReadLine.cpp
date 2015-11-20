#include "LKH.h"

/*      
 * The ReadLine function reads the next input line from a file. The function
 * handles the problem that an input line may be terminated by a carriage
 * return, a newline, both, or EOF.
 */

static boost::thread_specific_ptr<vector<char> > Buffer;
static boost::thread_specific_ptr<int> MaxBuffer;

static int EndOfLine(FILE * InputFile, int c)
{
    int EOL = (c == '\r' || c == '\n');
    if (c == '\r') {
        c = fgetc(InputFile);
        if (c != '\n' && c != EOF)
            ungetc(c, InputFile);
    }
    return EOL;
}

char *LKH::LKHAlg::ReadLine(FILE * InputFile)
{
    int i, c;

    if (!Buffer.get())
	{
		MaxBuffer.reset(new int(80));
		Buffer.reset(new vector<char>(*MaxBuffer));
	}
    for (i = 0; (c = fgetc(InputFile)) != EOF && !EndOfLine(InputFile, c);
         i++) {
        if (i >= *MaxBuffer - 1) {
            *MaxBuffer *= 2;
			Buffer.reset(new vector<char>(*MaxBuffer));
        }
        (*Buffer.get())[i] = (char) c;
    }
    (*Buffer.get())[i] = '\0';
    if (!LastLine || (int) strlen(LastLine) < i) {
        free(LastLine);
        assert(LastLine = (char *) malloc((i + 1) * sizeof(char)));
    }
    strcpy(LastLine, &(*Buffer.get())[0]);
    return c == EOF && i == 0 ? 0 : &(*Buffer.get())[0];
}
