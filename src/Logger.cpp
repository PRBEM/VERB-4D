/**
 * \file Logger.cpp
 * \brief Logs the messages sent, including information, warnings, and errors
 */


#include "Logger.h"

/// creates logger that will write to log.txt
void Logger::createInstance()
{	
	delete sInstance;
	sInstance = new Logger;
	mDebugLevel = DEBUG_LEVEL_ERROR;
	logFile.open("log.txt");
}

/// deletes instance of logger and closes log file ('log.txt')
void Logger::deleteInstance()
{
	delete sInstance;	
	logFile.flush();
	logFile.close();
}

/// Can set the Debug Level
void Logger::setDebugLevel(Logger::DebugLevel level)
{
	mDebugLevel = level;
}

/// can signify errors while writing message
void Logger::writeError(const std::string& message)
{
	// if level of error is greator than set debug level
	if (mDebugLevel >= DEBUG_LEVEL_ERROR) {
		// type out [Error]
		const char *type = "[Error] ";
		string str = type + message;
		// write out message
		writeMessage(str, MessageType::MESSAGE_ERROR); // Add [Error] //message = message.append("[Error]: ");....
	}
#ifdef _DEBUG
	assert(0);
#endif	
	// exits/fails upon error
	exit(EXIT_FAILURE);
}

/// can signify warning while writing message
void Logger::writeWarning(const std::string& message)
{
	// if level of warning is greator than set debug level
	if (mDebugLevel >= DEBUG_LEVEL_WARNING) {
		// type out [Error]
		const char *type = "[Warning] ";
		string str = type + message;
		// write out message
		writeMessage(str, MessageType::MESSAGE_WARNING);
	}
}

/// Write out message to logger of type MESSAGE_INFO if debug level allows it
void Logger::writeMessage(const std::string& message)
{
	if (mDebugLevel >= DEBUG_LEVEL_MESSAGE) {

		writeMessage(message, MessageType::MESSAGE_INFO);
	}
}

/// Add message directly to the end of a file
void Logger::appendToFile(const std::string& message){
	logFile << message << flush;
}

/// Write out message to logger with the inputted message type 
void Logger::writeMessage(const std::string& message, MessageType type)
{
	std::string text(message);
			//std::replace(text.begin(), text.end(), '\n', ' '); // Fix it later
	// Send message to standard console out
	cout << message;
	// add message to end of file
	appendToFile(text);
#ifdef _DEBUG
	writeIDEDebugString(text, type);
#endif
}

//void Log_Windows::writeIDEDebugString(const std::string& message, MessageType type)
void Logger::writeIDEDebugString(const std::string& message, MessageType type) 
{
	//cerr << message.c_str() << endl;
}

Logger::Streamer::Streamer(Logger::MessageType messageType)
: std::ostream(new StringBuffer(messageType))
{
}

/// Deconstructor - deletes stream buffer
Logger::Streamer::~Streamer()
{
	delete rdbuf();
}

Logger::Streamer::StringBuffer::StringBuffer(Logger::MessageType messageType)
: mMessageType(messageType)
{
}

/// Writes out any unwritten characters to output if string buffers not synchronized, else does nothing
Logger::Streamer::StringBuffer::~StringBuffer()
{
	pubsync();
}

// Creates message, warning and error variables for Logger
Logger::Streamer Logger::message(Logger::MESSAGE_INFO);
Logger::Streamer Logger::warning(Logger::MESSAGE_WARNING);
Logger::Streamer Logger::error(Logger::MESSAGE_ERROR);

/// Checks to make sure Logger instance is in sync with correct message type
int Logger::Streamer::StringBuffer::sync()
{
	if (Logger::sInstance == NULL) {
		return 0;
	}
	std::string text(str());
	if (text.empty()) {
		return 0;
	}
	str("");
	switch (mMessageType) {
	case MESSAGE_INFO:
		Logger::sInstance->writeMessage(text);
		break;

	case MESSAGE_WARNING:
		Logger::sInstance->writeWarning(text);
		break;

	case MESSAGE_ERROR:
		Logger::sInstance->writeError(text);
		break;
	}
	return 0;
};

Logger* Logger::sInstance = 0;
Logger::DebugLevel Logger::mDebugLevel = Logger::DebugLevel::DEBUG_LEVEL_DISABLED;
ofstream Logger::logFile;



/*
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
*/