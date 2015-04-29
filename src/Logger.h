using namespace std;

#include <iostream>
#include <streambuf>
#include <fstream>
#include <sstream>


class Logger
{
public:

	ofstream logFile;

	Logger();
	~Logger();

	//Logger& operator<< (std::ostream& (*pfun)(std::ostream&));
	/*Logger& operator<<(const char* sb);
	friend Logger operator<<(Logger &log, int i);*/
	
	//template <typename ANY>
	//friend Logger& operator<<(Logger& logClass, ANY val);
	
	
	friend Logger& operator<<(Logger &log, int i);
	friend Logger& operator<<(Logger &log, long i);
	friend Logger& operator<<(Logger &log, float i);
	friend Logger& operator<<(Logger &log, double i);
	friend Logger& operator<<(Logger &log, const char* sb);
	friend Logger& operator<<(Logger &log, string sb);

	template <class T>
	friend Logger& operator<<(Logger &log, T i);
	Logger& operator << (std::ostream& (*pfun)(std::ostream&));
};

