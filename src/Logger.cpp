#include "Logger.h"


Logger::Logger()
{
	// Open file for log
	if(!logFile.is_open())
		logFile.open("log.txt", std::ios_base::out | std::ios_base::trunc);

};

Logger::~Logger()
{
	// Close file bofore die.
	logFile.flush();
	logFile.close();
};


Logger& Logger::operator << (std::ostream& (*pfun)(std::ostream&))
{
	pfun(logFile);
	pfun(std::cout);
	return (*this);
}

// friend functions

template <class T>
Logger& operator<<(Logger &log, T i)
{
	log.logFile << i;
	std::cout << i;
	return log;
};

Logger& operator<<(Logger &log, const char* sb)
{
	log.logFile << sb;
	std::cout << sb;
	return log;
};

Logger& operator<<(Logger &log, string sb)
{
	log.logFile << sb;
	std::cout << sb;
	return log;
};



Logger& operator<<(Logger &log, int i)
{
	log.logFile << i;
	std::cout << i;
	return log;
};

Logger& operator<<(Logger &log, long i)
{
	log.logFile << i;
	std::cout << i;
	return log;
};

Logger& operator<<(Logger &log, float i)
{
	log.logFile << i;
	std::cout << i;
	return log;
};

Logger& operator<<(Logger &log, double i)
{
	log.logFile << i;
	std::cout << i;
	return log;
};
