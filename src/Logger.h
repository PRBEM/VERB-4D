#pragma once

#include <iostream>
#include <streambuf>
#include <fstream>
#include <sstream>
#include <stdio.h> 
#include <assert.h>
using namespace std;

class Logger
{
public:

	

	enum DebugLevel
	{
		DEBUG_LEVEL_DISABLED = 0,
		DEBUG_LEVEL_MESSAGE = 1,
		DEBUG_LEVEL_WARNING = 2,
		DEBUG_LEVEL_ERROR = 3,		
	};

	enum MessageType
	{
		MESSAGE_INFO = 1,
		MESSAGE_WARNING = 2,
		MESSAGE_ERROR = 3,
	};

	static void createInstance();
	static void deleteInstance();
	static void setDebugLevel(DebugLevel level);



	class Streamer : public std::ostream {
	public:

		Streamer(Logger::MessageType messageType);
		~Streamer();

	private:
		class StringBuffer : public std::stringbuf {
		public:
			StringBuffer(Logger::MessageType messageType);
			~StringBuffer();
			virtual int sync();

		private:
			Logger::MessageType mMessageType;
		};
	};

	static Streamer message;
	static Streamer warning;
	static Streamer error;

	

	protected:
	Logger() {};
	Logger(const Logger& thisobject);
	Logger& operator=(const Logger&);
	~Logger() {};

	void writeIDEDebugString(const std::string& message, MessageType type);

	

	private:
		void writeMessage(const std::string& message);
		void writeWarning(const std::string& message);
		void writeError(const std::string& message);

		void appendToFile(const std::string& message);

		void writeMessage(const std::string& message, MessageType type);

		static DebugLevel mDebugLevel;
		static Logger* sInstance; // ??
		static ofstream logFile;
	};


	

	//Logger& operator<< (std::ostream& (*pfun)(std::ostream&));
	/*Logger& operator<<(const char* sb);
	friend Logger operator<<(Logger &log, int i);*/
	
	//template <typename ANY>
	//friend Logger& operator<<(Logger& logClass, ANY val);
	
	
	/*friend Logger& operator<<(Logger &log, int i);
	friend Logger& operator<<(Logger &log, long i);
	friend Logger& operator<<(Logger &log, float i);
	friend Logger& operator<<(Logger &log, double i);
	friend Logger& operator<<(Logger &log, const char* sb);
	friend Logger& operator<<(Logger &log, string sb);

	template <class T>
	friend Logger& operator<<(Logger &log, T i);
	Logger& operator << (std::ostream& (*pfun)(std::ostream&));*/
//};

