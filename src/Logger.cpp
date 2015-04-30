#include "Logger.h"

void Logger::createInstance()
{	
	delete sInstance;
	sInstance = new Logger;
	mDebugLevel = DEBUG_LEVEL_ERROR;
	logFile.open("log.txt");
}


void Logger::deleteInstance()
{
	delete sInstance;	
	logFile.flush();
	logFile.close();
}

void Logger::setDebugLevel(Logger::DebugLevel level)
{
	mDebugLevel = level;
}

void Logger::writeError(const std::string& message)
{
	if (mDebugLevel >= DEBUG_LEVEL_ERROR) {

		const char *type = "[Error] ";
		string str = type + message;

		writeMessage(str, MessageType::MESSAGE_ERROR); // Add [Error] //message = message.append("[Error]: ");....
	}
#ifdef _DEBUG
	assert(0);
#endif	
	exit(EXIT_FAILURE);
}

void Logger::writeWarning(const std::string& message)
{
	if (mDebugLevel >= DEBUG_LEVEL_WARNING) {
		const char *type = "[Warning] ";
		string str = type + message;

		writeMessage(str, MessageType::MESSAGE_WARNING);
	}
}

void Logger::writeMessage(const std::string& message)
{
	if (mDebugLevel >= DEBUG_LEVEL_MESSAGE) {

		writeMessage(message, MessageType::MESSAGE_INFO);
	}
}

void Logger::appendToFile(const std::string& message){
	logFile << message << flush;
}

void Logger::writeMessage(const std::string& message, MessageType type)
{
	std::string text(message);
	//std::replace(text.begin(), text.end(), '\n', ' '); // Fix it later
	cout << message;
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

Logger::Streamer::~Streamer()
{
	delete rdbuf();
}

Logger::Streamer::StringBuffer::StringBuffer(Logger::MessageType messageType)
: mMessageType(messageType)
{
}

Logger::Streamer::StringBuffer::~StringBuffer()
{
	pubsync();
}

Logger::Streamer Logger::message(Logger::MESSAGE_INFO);
Logger::Streamer Logger::warning(Logger::MESSAGE_WARNING);
Logger::Streamer Logger::error(Logger::MESSAGE_ERROR);

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